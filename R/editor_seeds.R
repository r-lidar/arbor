#' Interactive Tree Seed ID Editor (RGL)
#'
#' Launches an interactive \pkg{rgl} 3D editor to manually merge tree seed
#' identifiers by clicking on points in a point cloud scene.
#'
#' The workflow is based on a simple state machine:
#' \enumerate{
#'   \item Select a \strong{KEEP} seed (middle mouse button)
#'   \item Select one or more \strong{CHANGE} seeds (middle mouse button)
#'   \item Click the KEEP seed again to validate the merge (middle mouse button)
#' }
#' Multiple red \strong{EXIT} markers are displayed around the scene to allow
#' easy termination of the editor (middle mouse button). Closing the RGL window or
#' pressing ESC during selection \strong{DOES NOT} exit safely.
#'
#' @param las A \code{LAS} used as
#'   background point cloud for visual context.
#' @param seeds A LAS containing seeds.
#' @param ... for debugging use only
#'
#' @return
#' The updated \code{seeds} object with modified \code{treeID} values.
#' @export
editor_seeds <- function(las, seeds, ...)
{
  # ============================================================================
  # SETUP FUNCTIONS
  # ============================================================================

  setup_palette <- function(ids)
  {
    unique_ids <- unique(ids)
    colors_palette <- pastel.colors(length(unique_ids))
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

  cat("Rendering background point cloud...\n")

  pc <- lidR::filter_poi(las, hag < 2)
  long_passages <- lidR::filter_poi(las, hag < 4, passage > 150)

  x <- plot_semantic(pc)
  lidR::plot(long_passages, add = x, pal = "gray", size = 2)

  seeds$X <- seeds$X - x[1]
  seeds$Y <- seeds$Y - x[2]

  # Render editable seeds
  seeds_obj_id <- rgl::points3d(seeds$X, seeds$Y, seeds$Z,
                                col = get_colors(seeds$treeID), size = 7)

  rgl::axes3d()
  rgl::title3d(main = "TreeID Editor", xlab = "X", ylab = "Y", zlab = "Z")

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
    seeds_obj_id <<- rgl::points3d(seeds$X, seeds$Y, seeds$Z,
                                   col = get_colors(seeds$treeID), size = 7)
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
    keep_marker <<- rgl::points3d(seed_data$X, seed_data$Y, seed_data$Z,
                                  col = "blue", size = 15)
    keep_text <<- rgl::text3d(seed_data$X, seed_data$Y, seed_data$Z,
                              texts = paste("KEEP:", seed_data$treeID),
                              col = "blue", cex = 1.5, adj = c(0, -1.5))
  }

  add_change_marker <- function(seed_data)
  {
    marker <- rgl::points3d(seed_data$X, seed_data$Y, seed_data$Z,
                            col = "red", size = 15)
    text_obj <- rgl::text3d(seed_data$X, seed_data$Y, seed_data$Z,
                            texts = paste("CHANGE:", seed_data$treeID),
                            col = "red", cex = 1.5, adj = c(0, -1.5))
    line_obj <- rgl::segments3d(c(keep_seed$X, seed_data$X),
                                c(keep_seed$Y, seed_data$Y),
                                c(keep_seed$Z, seed_data$Z),
                                col = "red", lwd = 3)

    change_markers[[length(change_markers) + 1]] <<- marker
    change_texts[[length(change_texts) + 1]] <<- text_obj
    change_lines[[length(change_lines) + 1]] <<- line_obj
  }

  # ============================================================================
  # UI UPDATE FUNCTIONS
  # ============================================================================

  update_status <- function(text, color = "black", bg = "gray95")
  {
    tkconfigure(status_label, text = paste("Status:", text),
                foreground = color, background = bg)
    tcltk::tcl("update", "idletasks")
  }

  update_selection_info <- function(text, bg = "lightyellow")
  {
    tkconfigure(selection_info_label, text = text, background = bg)
    tcltk::tcl("update", "idletasks")
  }

  update_button_states <- function()
  {
    if (state == "IDLE")
    {
      tkconfigure(btn_start_selection, state = "normal")
      tkconfigure(btn_exit, state = "normal")
    }
    else  # SELECT_KEEP or SELECT_CHANGE
    {
      tkconfigure(btn_start_selection, state = "disabled")
      tkconfigure(btn_exit, state = "disabled")
    }
    tcltk::tcl("update", "idletasks")
  }

  # ============================================================================
  # SELECTION LOGIC (BLOCKING LOOP)
  # ============================================================================

  selection_loop <- function()
  {
    cat("\n======================================================\n")
    cat("SELECTION MODE ACTIVATED\n")
    cat("Step 1: Middle-click KEEP seed (will turn blue)\n")
    cat("Step 2: Middle-click CHANGE seeds (will turn red)\n")
    cat("Step 3: Click CHANGE seed again to unselect it\n")
    cat("Step 4: Click KEEP seed again to VALIDATE merge\n")
    cat("======================================================\n\n")

    state <<- "SELECT_KEEP"
    update_status("Selection mode ACTIVE - Click KEEP seed in 3D viewer", "blue", "lightblue")
    update_selection_info("Step 1: Click a KEEP seed in the 3D viewer", "lightblue")

    repeat
    {
      # Check if RGL window was closed
      if (length(rgl::rgl.dev.list()) == 0 || rgl::rgl.cur() == 0)
      {
        cat("\n✓ RGL window closed - exiting\n")
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
          rgl::identify3d(seeds$X, seeds$Y, seeds$Z, n = 1,
                          plot = FALSE, buttons = "middle", tolerance = 40)
        },
        error = function(e) NULL
      )

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

      cat(sprintf("\nClicked seed index: %d, treeID: %s\n", sel, selected_id))
      cat(sprintf("Current state: %s\n", state))

      # ======================================================================
      # STATE MACHINE LOGIC
      # ======================================================================

      if (state == "SELECT_KEEP")
      {
        keep_seed <<- selected_seed
        state <<- "SELECT_CHANGE"

        add_keep_marker(keep_seed)

        update_status(sprintf("KEEP seed selected: %s", keep_seed$treeID), "blue", "lightblue")
        update_selection_info(sprintf("KEEP: %s | Click CHANGE seeds (or KEEP again to finish)",
                                      keep_seed$treeID), "lightblue")

        cat(sprintf(">>> KEEP seed selected: ID %s\n", keep_seed$treeID))
      }
      else if (state == "SELECT_CHANGE")
      {
        # Check if clicked KEEP seed again → VALIDATE MERGE
        if (selected_id == keep_seed$treeID)
        {
          if (length(change_seeds) == 0)
          {
            # No CHANGE seeds selected → cancel and return to IDLE
            update_status("No CHANGE seeds - canceling selection", "purple", "lavender")
            cat(">>> No CHANGE seeds to merge. Returning to IDLE.\n")

            clear_keep_markers()
            clear_change_markers()

            keep_seed <<- NULL
            change_seeds <<- list()
            state <<- "IDLE"
            selection_active <<- FALSE

            update_selection_info("Ready to start new selection", "lightyellow")
            update_button_states()

            Sys.sleep(1)
            break
          }

          # PERFORM MERGE
          change_ids <- sapply(change_seeds, function(s) s$treeID)

          update_status(sprintf("Merging %d seeds into %s...",
                                length(change_ids), keep_seed$treeID),
                        "yellow", "lightyellow")

          cat("\n=======================================================\n")
          cat(sprintf("MERGING: %s → %s\n",
                      paste(change_ids, collapse = ", "), keep_seed$treeID))
          cat("=======================================================\n")

          # Update data
          for (change_id in change_ids)
          {
            seeds$treeID[seeds$treeID == change_id] <<- keep_seed$treeID
          }

          # Refresh visuals
          refresh_seeds_display()

          update_status(sprintf("✓ Merged %d seeds successfully!", length(change_ids)),
                        "dark green", "lightgreen")
          update_selection_info(sprintf("✓ Merged %d IDs into %s",
                                        length(change_ids), keep_seed$treeID),
                                "lightgreen")

          cat(sprintf("✓ Successfully merged %d IDs\n\n", length(change_ids)))

          # CLEANUP & RESET STATE
          clear_keep_markers()
          clear_change_markers()

          keep_seed <<- NULL
          change_seeds <<- list()
          state <<- "IDLE"
          selection_active <<- FALSE

          Sys.sleep(1.5)

          update_status("Ready", "black")
          update_selection_info("Ready to start new selection", "lightyellow")
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
            update_status(sprintf("Unselecting CHANGE seed: %s", selected_id), "purple")
            cat(sprintf(">>> CHANGE seed unselected: ID %s\n", selected_id))

            # Remove visual elements
            rgl::pop3d(id = change_markers[[existing_index]])
            rgl::pop3d(id = change_texts[[existing_index]])
            rgl::pop3d(id = change_lines[[existing_index]])

            # Remove from lists
            change_seeds <<- change_seeds[-existing_index]
            change_markers <<- change_markers[-existing_index]
            change_texts <<- change_texts[-existing_index]
            change_lines <<- change_lines[-existing_index]

            update_selection_info(sprintf("KEEP: %s | CHANGE seeds: %d (click KEEP to validate)",
                                          keep_seed$treeID, length(change_seeds)),
                                  "lightblue")
          }
          else
          {
            # ADD NEW CHANGE SEED
            change_seeds[[length(change_seeds) + 1]] <<- selected_seed
            add_change_marker(selected_seed)

            cat(sprintf(">>> CHANGE seed #%d added: ID %s\n",
                        length(change_seeds), selected_id))

            update_status(sprintf("CHANGE seed added: %s", selected_id), "orange", "lightyellow")
            update_selection_info(sprintf("KEEP: %s | CHANGE seeds: %d (click KEEP to validate)",
                                          keep_seed$treeID, length(change_seeds)),
                                  "lightyellow")
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
      update_selection_info("Ready to start new selection", "lightyellow")
      update_button_states()
    }
  }

  exit_application <- function()
  {
    cat("\n✓ Exiting editor\n")

    selection_active <<- FALSE
    app_running <<- FALSE

    # Close RGL window
    if (length(rgl::rgl.dev.list()) > 0)
    {
      tryCatch(rgl::close3d(), error = function(e) invisible(NULL))
    }

    # Close tcltk window
    tryCatch(tkdestroy(tt), error = function(e) invisible(NULL))
  }

  # ============================================================================
  # BUILD TCLTK UI
  # ============================================================================

  tt <- tcltk::tktoplevel()
  tcltk::tkwm.title(tt, "Arbor Studio - Seed Editor")

  # Bind window close event
  tcltk::tkwm.protocol(tt, "WM_DELETE_WINDOW", exit_application)

  # Main frame
  main_frame <- tcltk::tkframe(tt, borderwidth = 2, relief = "groove")
  tcltk::tkpack(main_frame, fill = "both", expand = TRUE, padx = 10, pady = 10)

  # Title
  title_label <- tcltk::tklabel(main_frame, text = "Tree Seed Editor",
                                font = tcltk::tkfont.create(size = 12, weight = "bold"))
  tcltk::tkpack(title_label, pady = c(5, 10))

  # Status frame
  status_frame <- tcltk::tkframe(main_frame, borderwidth = 1, relief = "sunken",
                                 background = "gray95")
  tcltk::tkpack(status_frame, fill = "x", padx = 5, pady = 5)

  status_label <- tcltk::tklabel(status_frame, text = "Status: Ready",
                                 anchor = "w", font = tcltk::tkfont.create(size = 10),
                                 background = "gray95")
  tcltk::tkpack(status_label, fill = "x", padx = 5, pady = 5)

  # Selection info frame
  info_frame <- tcltk::tkframe(main_frame, borderwidth = 1, relief = "sunken",
                               background = "lightyellow")
  tcltk::tkpack(info_frame, fill = "x", padx = 5, pady = 5)

  selection_info_label <- tcltk::tklabel(info_frame,
                                         text = "Ready to start new selection",
                                         anchor = "w",
                                         font = tcltk::tkfont.create(size = 9),
                                         background = "lightyellow")
  tcltk::tkpack(selection_info_label, fill = "x", padx = 5, pady = 5)

  # Separator
  tcltk::tkpack(tcltk::tkframe(main_frame, height = 2, borderwidth = 1,
                               relief = "sunken"), fill = "x", pady = 10)

  # Button section title
  button_section_label <- tcltk::tklabel(main_frame, text = "Actions",
                                         font = tcltk::tkfont.create(size = 10,
                                                                     weight = "bold"))
  tcltk::tkpack(button_section_label, pady = c(0, 5))

  # Button 1: Start Selection
  btn_start_selection <- tcltk::tkbutton(main_frame,
                                         text = "🎯 Start Selection",
                                         command = start_selection,
                                         width = 25, pady = 10,
                                         fg = "dark blue",
                                         font = tcltk::tkfont.create(size = 10, weight = "bold"))
  tcltk::tkpack(btn_start_selection, fill = "x", padx = 10, pady = 8)

  # Separator
  tcltk::tkpack(tcltk::tkframe(main_frame, height = 2, borderwidth = 1,
                               relief = "sunken"), fill = "x", pady = 10)

  # Button 2: Exit
  btn_exit <- tcltk::tkbutton(main_frame,
                              text = "❌ Exit",
                              command = exit_application,
                              width = 25, pady = 8,
                              fg = "red")
  tcltk::tkpack(btn_exit, fill = "x", padx = 10, pady = 5)

  # Help text
  help_text <- "
📋 Workflow:
1. Click 'Start Selection' button to begin
2. Middle-click KEEP seed in 3D viewer (turns blue)
3. Middle-click CHANGE seeds to merge (turn red)
4. Click CHANGE seeds again to unselect them
5. Click KEEP seed again to VALIDATE and merge
6. Click 'Exit' to close editor

⚠️ Note:
• Buttons are disabled during selection mode
• Validation happens by clicking KEEP seed again
• Selection mode exits after validation
"
  help_label <- tcltk::tklabel(main_frame, text = help_text,
                               justify = "left",
                               font = tcltk::tkfont.create(size = 8),
                               fg = "gray40")
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
      cat("\n✓ RGL window closed - exiting\n")
      exit_application()
      return()
    }

    tcltk::tcl("after", 500, check_rgl_window)
  }

  check_rgl_window()

  # Wait for window to close
  tcltk::tkwait.window(tt)

  # Restore original coordinates
  seeds$X <- seeds$X + x[1]
  seeds$Y <- seeds$Y + x[2]

  return(seeds)
}
