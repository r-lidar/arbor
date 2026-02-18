#' @export
#' @rdname arbor-studio
arbor_studio_seeds <- function(las, seeds, ...)
{
  # ============================================================================
  # SETUP FUNCTIONS
  # ============================================================================

  setup_palette <- function(ids)
  {
    unique_ids <- unique(ids)
    colors_palette <- lidR::pastel.colors(length(unique_ids))
    names(colors_palette) <- as.character(unique_ids)
    return(colors_palette)
  }

  get_colors <- function(ids)
  {
    cols <- colors_palette[as.character(ids)]
    cols[is.na(cols)] <- "white"
    return(cols)
  }

  colors_palette <- setup_palette(seeds$treeID)

  # ============================================================================
  # RENDER 3D SCENE
  # ============================================================================

  passage <- hag <- NULL

  pc <- lidR::filter_poi(las, hag < 2)
  long_passages <- lidR::filter_poi(las, hag < 4, passage > 150)

  x <- plot_semantic(pc)
  lidR::plot(long_passages, add = x, pal = "gray", size = 2)
  rgl::axes3d()
  rgl::title3d(main = "TreeID Editor", xlab = "X", ylab = "Y", zlab = "Z")
  dim = rgl::par3d("windowRect")
  dim[3] = dim[3] + 256
  dim[4] = dim[4] + 256
  rgl::par3d(windowRect = dim)

  seeds$X <- seeds$X - x[1]
  seeds$Y <- seeds$Y - x[2]

  # Render editable seeds
  seeds_obj_id <- rgl::points3d(seeds$X, seeds$Y, seeds$Z, col = get_colors(seeds$treeID), size = 7)

  # ============================================================================
  # STATE VARIABLES
  # ============================================================================

  state <- "IDLE"  # States: IDLE, SELECT_KEEP, SELECT_CHANGE
  keep_seed <- NULL
  change_seeds <- list()

  # Visual tracking
  keep_marker <- NULL
  keep_text <- NULL
  change_markers <- list()
  change_texts <- list()
  change_lines <- list()

  # Control variables
  selection_active <- FALSE
  app_running <- TRUE

  # ============================================================================
  # VISUAL UPDATE FUNCTIONS
  # ============================================================================

  refresh_seeds_display <- function()
  {
    if (!is.null(seeds_obj_id))
    {
      tryCatch(rgl::pop3d(id = seeds_obj_id), error = function(e) invisible(NULL))
    }
    colors_palette <<- setup_palette(seeds$treeID)
    seeds_obj_id <<- rgl::points3d(seeds$X, seeds$Y, seeds$Z, col = get_colors(seeds$treeID), size = 7)
  }

  clear_keep_markers <- function()
  {
    if (!is.null(keep_marker)) rgl::pop3d(id = keep_marker)
    if (!is.null(keep_text)) rgl::pop3d(id = keep_text)
    keep_marker <<- NULL
    keep_text <<- NULL
  }

  clear_change_markers <- function()
  {
    for (m in change_markers) rgl::pop3d(id = m)
    for (t in change_texts) rgl::pop3d(id = t)
    for (l in change_lines) rgl::pop3d(id = l)
    change_markers <<- list()
    change_texts <<- list()
    change_lines <<- list()
  }

  add_keep_marker <- function(seed_data)
  {
    keep_marker <<- rgl::points3d(seed_data$X, seed_data$Y, seed_data$Z, col = "blue", size = 15)
    keep_text <<- rgl::text3d(seed_data$X, seed_data$Y, seed_data$Z, texts = paste("REF:", seed_data$treeID), col = "blue", cex = 1.5, adj = c(0, -1.5))
  }

  add_change_marker <- function(seed_data)
  {
    marker <- rgl::points3d(seed_data$X, seed_data$Y, seed_data$Z, col = "red", size = 15)
    text_obj <- rgl::text3d(seed_data$X, seed_data$Y, seed_data$Z, texts = paste("CHANGE:", seed_data$treeID), col = "red", cex = 1.5, adj = c(0, -1.5))
    line_obj <- rgl::segments3d(c(keep_seed$X, seed_data$X), c(keep_seed$Y, seed_data$Y), c(keep_seed$Z, seed_data$Z), col = "red", lwd = 3)

    change_markers[[length(change_markers) + 1]] <<- marker
    change_texts[[length(change_texts) + 1]] <<- text_obj
    change_lines[[length(change_lines) + 1]] <<- line_obj
  }

  # ============================================================================
  # UI UPDATE FUNCTIONS
  # ============================================================================

  update_status <- function(text, color = "black", bg = "gray95")
  {
    tcltk::tkconfigure(status_label, text = paste("Status:", text), foreground = color, background = bg)
    tcltk::tcl("update", "idletasks")
  }

  update_selection_info <- function(text, bg = "lightyellow")
  {
    tcltk::tkconfigure(selection_info_label, text = text, background = bg)
    tcltk::tcl("update", "idletasks")
  }

  update_button_states <- function()
  {
    if (state == "IDLE")
    {
      tcltk::tkconfigure(btn_start_selection, state = "normal")
      tcltk::tkconfigure(btn_exit, state = "normal")
    }
    else  # SELECT_KEEP or SELECT_CHANGE
    {
      tcltk::tkconfigure(btn_start_selection, state = "disabled")
      tcltk::tkconfigure(btn_exit, state = "disabled")
    }
    tcltk::tcl("update", "idletasks")
  }

  # ============================================================================
  # SELECTION LOGIC (BLOCKING LOOP)
  # ============================================================================

  selection_loop <- function()
  {
    state <<- "SELECT_KEEP"
    update_status("Selection mode ACTIVE")
    update_selection_info("Click a reference seed in the 3D viewer")

    repeat
    {
      # Check if RGL window was closed
      if (length(rgl::rgl.dev.list()) == 0 || rgl::rgl.cur() == 0)
      {
        #cat("\nRGL window closed - exiting\n")
        selection_active <<- FALSE
        app_running <<- FALSE
        break
      }

      # Check if app was closed
      if (!app_running)
      {
        break
      }

      # Update tcltk events
      tcltk::tcl("update", "idletasks")

      # Wait for user click (BLOCKING)
      sel <- tryCatch(
      {
        rgl::identify3d(seeds$X, seeds$Y, seeds$Z, n = 1, plot = FALSE, buttons = "middle", tolerance = 40)
      },
      error = function(e) NULL)

      if (is.null(sel) || length(sel) == 0)
      {
        next
      }

      # Validate selection
      if (sel < 1 || sel > nrow(seeds))
      {
        update_status("Invalid selection, try again", "orange")
        next
      }

      selected_seed <- seeds@data[sel, ]
      selected_id <- selected_seed$treeID

      #cat(sprintf("\nClicked seed index: %d, treeID: %s\n", sel, selected_id))
      #cat(sprintf("Current state: %s\n", state))

      # ======================================================================
      # STATE MACHINE LOGIC
      # ======================================================================

      if (state == "SELECT_KEEP")
      {
        keep_seed <<- selected_seed
        state <<- "SELECT_CHANGE"

        add_keep_marker(keep_seed)

        update_status(sprintf("KEEP seed selected: %s", keep_seed$treeID))
        update_selection_info("Click CHANGE seeds")
      }
      else if (state == "SELECT_CHANGE")
      {
        # Check if clicked KEEP seed again -> VALIDATE MERGE
        if (selected_id == keep_seed$treeID)
        {
          if (length(change_seeds) == 0)
          {
            # No CHANGE seeds selected -> cancel and return to IDLE
            update_status("No CHANGE seeds - canceling selection")

            clear_keep_markers()
            clear_change_markers()

            keep_seed <<- NULL
            change_seeds <<- list()
            state <<- "IDLE"
            selection_active <<- FALSE

            update_selection_info("Ready to start new selection")
            update_button_states()

            Sys.sleep(1)
            break
          }

          # PERFORM MERGE
          change_ids <- sapply(change_seeds, function(s) s$treeID)

          update_status(sprintf("Merging %d seeds into %s...",  length(change_ids), keep_seed$treeID))

          # Update data
          for (change_id in change_ids)
          {
            seeds$treeID[seeds$treeID == change_id] <<- keep_seed$treeID
          }

          # Refresh visuals
          refresh_seeds_display()

          update_status(sprintf("Merged %d seeds successfully!", length(change_ids)), "darkgreen", "lightgreen")
          update_selection_info(sprintf("Merged %d IDs into %s", length(change_ids), keep_seed$treeID), "lightgreen")

          #cat(sprintf("Successfully merged %d IDs\n\n", length(change_ids)))

          # CLEANUP & RESET STATE
          clear_keep_markers()
          clear_change_markers()

          keep_seed <<- NULL
          change_seeds <<- list()
          state <<- "IDLE"
          selection_active <<- FALSE

          Sys.sleep(1.5)

          update_status("Ready", "black")
          update_selection_info("Ready to start new selection")
          update_button_states()

          break
        }
        else
        {
          # Check if this CHANGE seed is already selected
          already_selected <- FALSE
          existing_index <- NULL

          if (length(change_seeds) > 0)
          {
            for (i in seq_along(change_seeds))
            {
              if (change_seeds[[i]]$treeID == selected_id)
              {
                already_selected <- TRUE
                existing_index <- i
                break
              }
            }
          }

          if (already_selected)
          {
            # UNSELECT THIS CHANGE SEED
            update_status(sprintf("Unselecting CHANGE seed: %s", selected_id))

            # Remove visual elements
            rgl::pop3d(id = change_markers[[existing_index]])
            rgl::pop3d(id = change_texts[[existing_index]])
            rgl::pop3d(id = change_lines[[existing_index]])

            # Remove from lists
            change_seeds <<- change_seeds[-existing_index]
            change_markers <<- change_markers[-existing_index]
            change_texts <<- change_texts[-existing_index]
            change_lines <<- change_lines[-existing_index]

            update_selection_info(sprintf("KEEP: %s | CHANGE seeds: %d (click KEEP to validate)", keep_seed$treeID, length(change_seeds)))
          }
          else
          {
            # ADD NEW CHANGE SEED
            change_seeds[[length(change_seeds) + 1]] <<- selected_seed
            add_change_marker(selected_seed)

            update_status(sprintf("CHANGE seed added: %s", selected_id))
            update_selection_info(sprintf("KEEP: %s | CHANGE seeds: %d (click KEEP to validate)", keep_seed$treeID, length(change_seeds)))
          }
        }
      }
    }

    selection_active <<- FALSE
  }

  # ============================================================================
  # BUTTON ACTIONS
  # ============================================================================

  start_selection <- function()
  {
    if (selection_active) return()

    selection_active <<- TRUE
    update_button_states()

    # Start blocking selection loop
    selection_loop()

    # When loop exits, update UI
    if (app_running)
    {
      update_status("Ready", "black")
      update_selection_info("Ready to start new selection")
      update_button_states()
    }
  }

  exit_application <- function()
  {
    window_open <<- FALSE
    selection_active <<- FALSE
    app_running <<- FALSE

    # Close tcltk window
    tryCatch(tcltk::tkdestroy(tt), error = function(e) invisible(NULL))
  }

  # ============================================================================
  # BUILD TCLTK UI
  # ============================================================================

  window_open <- TRUE

  tt <- tcltk::tktoplevel()
  tcltk::tkwm.title(tt, "Arbor Studio - Seed Editor")

  # Bind window close event
  #tcltk::tkwm.protocol(tt, "WM_DELETE_WINDOW", prevent_close)

  # Main frame
  main_frame <- tcltk::tkframe(tt, borderwidth = 2, relief = "groove")
  tcltk::tkpack(main_frame, fill = "both", expand = TRUE, padx = 10, pady = 10)

  # Title
  title_label <- tcltk::tklabel(main_frame, text = "Tree Seed Editor", font = tcltk::tkfont.create(size = 12, weight = "bold"))
  tcltk::tkpack(title_label, pady = c(5, 10))

  # Status frame
  status_frame <- tcltk::tkframe(main_frame, borderwidth = 1, relief = "sunken")
  tcltk::tkpack(status_frame, fill = "x", padx = 5, pady = 5)

  status_label <- tcltk::tklabel(status_frame, text = "Status: Ready", anchor = "w", font = tcltk::tkfont.create(size = 10))
  tcltk::tkpack(status_label, fill = "x", padx = 5, pady = 5)

  # Selection info frame
  info_frame <- tcltk::tkframe(main_frame, borderwidth = 1, relief = "sunken")
  tcltk::tkpack(info_frame, fill = "x", padx = 5, pady = 5)

  selection_info_label <- tcltk::tklabel(info_frame, text = "Ready to start new selection", anchor = "w", font = tcltk::tkfont.create(size = 9))
  tcltk::tkpack(selection_info_label, fill = "x", padx = 5, pady = 5)

  # Separator
  tcltk::tkpack(tcltk::tkframe(main_frame, height = 2, borderwidth = 1, relief = "sunken"), fill = "x", pady = 10)

  # Button section title
  button_section_label <- tcltk::tklabel(main_frame, text = "Actions", font = tcltk::tkfont.create(size = 10, weight = "bold"))
  tcltk::tkpack(button_section_label, pady = c(0, 5))

  # Button 1: Start Selection
  btn_start_selection <- tcltk::tkbutton(main_frame, text = "Start selection", command = start_selection, pady = 10)
  tcltk::tkpack(btn_start_selection, fill = "x", padx = 10, pady = 8)

  # Button 2: Exit
  btn_exit <- tcltk::tkbutton(main_frame, text = "Exit and return", command = exit_application, pady = 8, fg = "red")
  tcltk::tkpack(btn_exit, fill = "x", padx = 10, pady = 5)

  # Separator
  tcltk::tkpack(tcltk::tkframe(main_frame, height = 2, borderwidth = 1, relief = "sunken"), fill = "x", pady = 10)


  # Help text
  help_text <- "
Workflow:
1. Click 'Start selection' button
2. Click (middle) to select one ref. seed in 3D
3. Click (middle) to select several CHANGE seeds
4. Click (middle) the ref. seed again to merge
5. Redo 1. for more seeds to merge
7. Click 'Exit' to close editor

Note:
* App. buttons are disabled during selection mode
* Validation happens by clicking REFERENCE seed again
* Clicking on an already selected seed = unselect
* DO NOT close the rgl windows while in selection mode
"
  help_label <- tcltk::tklabel(main_frame, text = help_text, justify = "left", font = tcltk::tkfont.create(size = 8), fg = "gray40")
  tcltk::tkpack(help_label, pady = c(5, 10))

  # Set initial button states
  update_button_states()

  # ============================================================================
  # MAIN LOOP
  # ============================================================================

  # Monitor RGL window
  check_rgl_window <- function()
  {
    if (!app_running) return()

    if (length(rgl::rgl.dev.list()) == 0 || rgl::rgl.cur() == 0)
    {
      #cat("\nRGL window closed - exiting\n")
      exit_application()
      return()
    }

    tcltk::tcl("after", 500, check_rgl_window)
  }

  check_rgl_window()

  # Wait for window to close
  while (window_open)
  {
    tryCatch(
    {
      tcltk::tcl("update")  # Process tcltk events
      Sys.sleep(0.05)       # Small delay to prevent CPU spinning
    },
    error = function(e)
    {
      window_open <<- FALSE  # Exit loop on any error
    })
  }

  # Close RGL window
  if (length(rgl::rgl.dev.list()) > 0)
  {
    tryCatch(rgl::close3d(), error = function(e) invisible(NULL))
  }

  # Restore original coordinates
  seeds$X <- seeds$X + x[1]
  seeds$Y <- seeds$Y + x[2]

  return(seeds)
}
