#' @export
#' @rdname arbor-studio
arbor_studio_instance <- function(las, ...)
{
  attribute = "treeID"

  # ---- Input validation ----
  if (!inherits(las, "LAS")) {
    stop("'las' must be a LAS object", call. = FALSE)
  }

  if (!attribute %in% names(las@data)) {
    stop(sprintf("Attribute '%s' not found in LAS object", attribute), call. = FALSE)
  }

  if (!interactive()) {
    stop("arbor_studio_instance() requires an interactive R session", call. = FALSE)
  }

  # ---- Setup environment ----
  # Use an environment to avoid global assignment and properly scope variables
  env <- new.env(parent = emptyenv())

  # Initialize state
  env$las_orig <- las
  env$las_orig$pointID <- seq_len(lidR::npoints(las))
  env$current_las <- las
  env$mouse_active <- FALSE
  env$current_state <- "Ready"
  env$selected_treeID <- NULL
  env$modified <- FALSE
  env$return_value <- NULL
  env$attribute <- attribute
  env$app_running <- TRUE
  env$center_x <- mean(las$X)
  env$center_y <- mean(las$Y)

  # ---- Helper functions ----

  redraw_scene <- function() {
    rgl::clear3d()
    plot_instance(env$current_las, add = c(env$center_x, env$center_y),  bg = "black")
  }

  nearest_point_treeID <- function(x, y, z, las_data) {
    xyz <- cbind(las_data$X, las_data$Y, las_data$Z)
    d2 <- rowSums((xyz - c(x, y, z))^2)
    las_data[[env$attribute]][which.min(d2)]
  }

  update_status <- function(status_text, color = "black") {
    env$current_state <- status_text
    if (!is.null(env$status_label)) {
      tryCatch({
        tcltk::tkconfigure(env$status_label,
                           text = paste("Status:", status_text),
                           fg = color)
        tcltk::tcl("update", "idletasks")
      }, error = function(e) {
        # Window might be closed
      })
    }
  }

  # ---- Exit application ----

  exit_application <- function() {
    update_status("Closing application...", "blue")

    env$app_running <- FALSE

    # Close RGL window
    tryCatch({
      rgl::close3d()
    }, error = function(e) {
      # Already closed
    })

    # Set return value
    if (env$modified) {
      env$return_value <- env$las_orig
    } else {
      env$return_value <- NULL
    }

    # Close Tcl/Tk window
    tryCatch({
      if (!is.null(env$tt)) {
        tcltk::tkdestroy(env$tt)
      }
    }, error = function(e) {
      # Already closed
    })
  }

  # ---- Extract context ----

  extract_and_display_context <- function() {
    if (env$mouse_active) {
      update_status("Please wait for current operation to finish", "red")
      return()
    }

    if (is.null(env$selected_treeID)) {
      update_status("No tree selected. Please select a tree first.", "red")
      return()
    }

    update_status("Extracting tree context...", "blue")

    tryCatch({
      env$current_las <- extract_tree_context(env$las_orig, env$selected_treeID)

      tree_info <- sprintf("Tree ID: %s | Height: %.2fm [CONTEXT EXTRACTED]",
                           env$selected_treeID,
                           max(env$current_las$Z) - min(env$current_las$Z))
      tcltk::tkconfigure(env$tree_info_label, text = tree_info)

      redraw_scene()

      # Update button states
      tcltk::tkconfigure(env$btn_pick, state = "disabled")
      tcltk::tkconfigure(env$btn_extract, state = "disabled")
      tcltk::tkconfigure(env$btn_export, state = "normal")
      tcltk::tkconfigure(env$btn_restore, state = "normal")
      tcltk::tkconfigure(env$btn_reassign_to_tree, state = "normal")
      tcltk::tkconfigure(env$btn_reassign_to_na, state = "normal")

      update_status(paste("Context extracted for Tree ID", env$selected_treeID), "green")

      tcltk::tcl("after", 2000, function() {
        if (grepl("Context extracted", env$current_state)) {
          update_status("Ready - Context extracted", "green")
        }
      })

    }, error = function(e) {
      update_status(paste("Error:", e$message), "red")
    })
  }

  # ---- Restore original ----

  restore_original <- function() {
    if (env$mouse_active) {
      update_status("Cannot restore while operation is active", "red")
      return()
    }

    update_status("Restoring original view...", "blue")

    env$current_las <- env$las_orig
    env$selected_treeID <- NULL
    redraw_scene()
    env$mouse_active <- FALSE

    tcltk::tkconfigure(env$tree_info_label, text = "All trees visible")

    # Update button states
    tcltk::tkconfigure(env$btn_pick, state = "normal")
    tcltk::tkconfigure(env$btn_extract, state = "disabled")
    tcltk::tkconfigure(env$btn_export, state = "disabled")
    tcltk::tkconfigure(env$btn_restore, state = "disabled")
    tcltk::tkconfigure(env$btn_reassign_to_tree, state = "disabled")
    tcltk::tkconfigure(env$btn_reassign_to_na, state = "disabled")

    update_status("Original view restored", "green")

    tcltk::tcl("after", 1500, function() {
      if (env$current_state == "Original view restored") {
        update_status("Ready", "black")
      }
    })
  }

  # ---- Point picking ----

  enable_tree_picking <- function() {
    if (env$mouse_active) {
      update_status("Tree picking already active", "orange")
      return()
    }

    env$mouse_active <- TRUE
    update_status("Click a point in the 3D viewer to select a tree...", "blue")

    # Disable buttons during selection
    tcltk::tkconfigure(env$btn_pick, state = "disabled")
    tcltk::tkconfigure(env$btn_extract, state = "disabled")
    tcltk::tkconfigure(env$btn_export, state = "disabled")
    tcltk::tkconfigure(env$btn_restore, state = "disabled")
    tcltk::tkconfigure(env$btn_reassign_to_tree, state = "disabled")
    tcltk::tkconfigure(env$btn_reassign_to_na, state = "disabled")
    tcltk::tkconfigure(env$btn_exit, state = "disabled")

    # Force update
    tcltk::tcl("update", "idletasks")

    ids <- rgl::rgl.ids(type = "shapes")

    if (nrow(ids) > 0) {
      vertices <- rgl::rgl.attrib(ids$id[1], "vertices")
      selected_idx <- rgl::identify3d(x = vertices[, 1],
                                      y = vertices[, 2],
                                      z = vertices[, 3],
                                      n = 1)

      if (length(selected_idx) > 0 && selected_idx > 0) {
        pt <- vertices[selected_idx, ]
        treeID <- nearest_point_treeID(pt[1] + env$center_x,
                                       pt[2] + env$center_y,
                                       pt[3],
                                       env$las_orig)
        env$selected_treeID <- treeID

        update_status(paste("Tree ID", treeID, "selected - Click 'Extract Context' to view"), "green")

        tree_info <- sprintf("Tree ID: %s selected [Context not extracted yet]", treeID)
        tcltk::tkconfigure(env$tree_info_label, text = tree_info)

        tcltk::tkconfigure(env$btn_pick, state = "normal")
        tcltk::tkconfigure(env$btn_extract, state = "normal")
        tcltk::tkconfigure(env$btn_exit, state = "normal")

        tcltk::tcl("after", 100, function() {
          tcltk::tkfocus(env$btn_extract)
        })

      } else {
        update_status("No point selected", "orange")
        env$selected_treeID <- NULL
        tcltk::tkconfigure(env$btn_pick, state = "normal")
        tcltk::tkconfigure(env$btn_exit, state = "normal")
      }
    } else {
      update_status("No shapes found in scene", "red")
      tcltk::tkconfigure(env$btn_pick, state = "normal")
      tcltk::tkconfigure(env$btn_exit, state = "normal")
    }

    env$mouse_active <- FALSE
  }

  # ---- Reassign points to tree ----

  reassign_points_to_current_tree <- function() {
    if (env$mouse_active) {
      update_status("Please wait for current operation to finish", "red")
      return()
    }

    if (is.null(env$selected_treeID)) {
      update_status("No tree selected. Cannot reassign.", "red")
      return()
    }

    if (is.null(env$current_las) || nrow(env$current_las@data) == 0 ||
        identical(env$current_las, env$las_orig)) {
      update_status("No context extracted. Please extract context first.", "red")
      return()
    }

    env$mouse_active <- TRUE
    update_status(paste("Draw a polygon in 3D viewer to select points for Tree ID",
                        env$selected_treeID), "blue")

    # Disable all buttons
    tcltk::tkconfigure(env$btn_pick, state = "disabled")
    tcltk::tkconfigure(env$btn_extract, state = "disabled")
    tcltk::tkconfigure(env$btn_export, state = "disabled")
    tcltk::tkconfigure(env$btn_restore, state = "disabled")
    tcltk::tkconfigure(env$btn_reassign_to_tree, state = "disabled")
    tcltk::tkconfigure(env$btn_reassign_to_na, state = "disabled")
    tcltk::tkconfigure(env$btn_exit, state = "disabled")

    # Force update
    tcltk::tcl("update", "idletasks")

    tryCatch({
      selection_fn <- rgl::select3d()
      selected <- selection_fn(env$current_las$X - env$center_x,
                               env$current_las$Y - env$center_y,
                               env$current_las$Z)

      if (sum(selected) > 0) {
        env$current_las@data[[env$attribute]][selected] <- env$selected_treeID

        ids <- env$current_las$pointID[selected]
        env$las_orig@data[[env$attribute]][ids] <- env$selected_treeID
        env$modified <- TRUE

        redraw_scene()

        update_status(paste(sum(selected), "points reassigned to Tree ID",
                            env$selected_treeID), "green")

        tree_info <- sprintf("Tree ID: %s | Points: %d | Height: %.2fm [MODIFIED]",
                             env$selected_treeID,
                             sum(env$current_las[[env$attribute]] == env$selected_treeID, na.rm = TRUE),
                             max(env$current_las$Z, na.rm = TRUE) - min(env$current_las$Z, na.rm = TRUE))
        tcltk::tkconfigure(env$tree_info_label, text = tree_info)

      } else {
        update_status("No points selected", "orange")
      }

    }, error = function(e) {
      update_status(paste("Selection error:", e$message), "red")
    })

    # Re-enable buttons
    tcltk::tkconfigure(env$btn_pick, state = "disabled")
    tcltk::tkconfigure(env$btn_extract, state = "disabled")
    tcltk::tkconfigure(env$btn_export, state = "normal")
    tcltk::tkconfigure(env$btn_restore, state = "normal")
    tcltk::tkconfigure(env$btn_reassign_to_tree, state = "normal")
    tcltk::tkconfigure(env$btn_reassign_to_na, state = "normal")
    tcltk::tkconfigure(env$btn_exit, state = "normal")

    env$mouse_active <- FALSE

    tcltk::tcl("after", 2000, function() {
      if (grepl("reassigned", env$current_state)) {
        update_status("Ready - Context extracted", "green")
      }
    })
  }

  # ---- Reassign points to NA ----

  reassign_points_to_na <- function() {
    if (env$mouse_active) {
      update_status("Please wait for current operation to finish", "red")
      return()
    }

    if (is.null(env$current_las) || nrow(env$current_las@data) == 0 ||
        identical(env$current_las, env$las_orig)) {
      update_status("No context extracted. Please extract context first.", "red")
      return()
    }

    env$mouse_active <- TRUE
    update_status("Draw a polygon in 3D viewer to select points for NA assignment", "blue")

    # Disable all buttons
    tcltk::tkconfigure(env$btn_pick, state = "disabled")
    tcltk::tkconfigure(env$btn_extract, state = "disabled")
    tcltk::tkconfigure(env$btn_export, state = "disabled")
    tcltk::tkconfigure(env$btn_restore, state = "disabled")
    tcltk::tkconfigure(env$btn_reassign_to_tree, state = "disabled")
    tcltk::tkconfigure(env$btn_reassign_to_na, state = "disabled")
    tcltk::tkconfigure(env$btn_exit, state = "disabled")

    # Force update
    tcltk::tcl("update", "idletasks")

    tryCatch({
      selection_fn <- rgl::select3d()
      selected <- selection_fn(env$current_las$X - env$center_x,
                               env$current_las$Y - env$center_y,
                               env$current_las$Z)

      if (sum(selected) > 0) {
        env$current_las@data[[env$attribute]][selected] <- NA

        ids <- env$current_las$pointID[selected]
        env$las_orig@data[[env$attribute]][ids] <- NA
        env$modified <- TRUE

        redraw_scene()

        update_status(paste(sum(selected), "points reassigned to NA"), "green")

        tree_info <- sprintf("Tree ID: %s | Points: %d | Height: %.2fm [MODIFIED]",
                             ifelse(is.null(env$selected_treeID), "NA", as.character(env$selected_treeID)),
                             nrow(env$current_las@data),
                             max(env$current_las$Z, na.rm = TRUE) - min(env$current_las$Z, na.rm = TRUE))
        tcltk::tkconfigure(env$tree_info_label, text = tree_info)

      } else {
        update_status("No points selected", "orange")
      }

    }, error = function(e) {
      update_status(paste("Selection error:", e$message), "red")
    })

    # Re-enable buttons
    tcltk::tkconfigure(env$btn_pick, state = "disabled")
    tcltk::tkconfigure(env$btn_extract, state = "disabled")
    tcltk::tkconfigure(env$btn_export, state = "normal")
    tcltk::tkconfigure(env$btn_restore, state = "normal")
    tcltk::tkconfigure(env$btn_reassign_to_tree, state = "normal")
    tcltk::tkconfigure(env$btn_reassign_to_na, state = "normal")
    tcltk::tkconfigure(env$btn_exit, state = "normal")

    env$mouse_active <- FALSE

    tcltk::tcl("after", 2000, function() {
      if (grepl("reassigned", env$current_state)) {
        update_status("Ready - Context extracted", "green")
      }
    })
  }

  # ---- Export tree ----

  export_selected_tree <- function() {
    if (env$mouse_active) {
      update_status("Please wait for current operation to finish", "red")
      return()
    }

    if (is.null(env$current_las) || nrow(env$current_las@data) == 0 ||
        identical(env$current_las, env$las_orig)) {
      update_status("No context extracted. Please extract context first.", "red")
      return()
    }

    update_status("Exporting tree...", "blue")

    tryCatch({
      filename <- tcltk::tclvalue(tcltk::tkgetSaveFile(
        title = "Export Tree",
        defaultextension = ".las",
        filetypes = "{{LAS files} {.las}} {{LAZ files} {.laz}} {{All files} *}"
      ))

      if (nchar(filename) > 0) {
        lidR::writeLAS(env$current_las, filename)
        update_status(paste("Tree exported to:", basename(filename)), "green")

        tcltk::tcl("after", 3000, function() {
          if (grepl("exported", env$current_state)) {
            update_status("Ready - Context extracted", "green")
          }
        })
      } else {
        update_status("Export cancelled", "orange")
        tcltk::tcl("after", 1500, function() {
          update_status("Ready - Context extracted", "green")
        })
      }

    }, error = function(e) {
      update_status(paste("Export error:", e$message), "red")
    })
  }

  # ---- Build UI ----

  # Initial plot
  rgl::open3d()
  rgl::bg3d("black")
  redraw_scene()

  # Create Tcl/Tk window
  env$tt <- tcltk::tktoplevel()
  tcltk::tkwm.title(env$tt, "Arbor Studio - Instance Module")
  tcltk::tkwm.protocol(env$tt, "WM_DELETE_WINDOW", exit_application)

  main_frame <- tcltk::tkframe(env$tt, borderwidth = 2, relief = "groove")
  tcltk::tkpack(main_frame, fill = "both", expand = TRUE, padx = 10, pady = 10)

  title_label <- tcltk::tklabel(main_frame, text = "Tree Instance Viewer",
                                font = tcltk::tkfont.create(size = 12, weight = "bold"))
  tcltk::tkpack(title_label, pady = c(5, 10))

  # Status frame
  status_frame <- tcltk::tkframe(main_frame, borderwidth = 1, relief = "sunken", background = "gray95")
  tcltk::tkpack(status_frame, fill = "x", padx = 5, pady = 5)

  env$status_label <- tcltk::tklabel(status_frame, text = "Status: Ready",
                                     anchor = "w", font = tcltk::tkfont.create(size = 10),
                                     background = "gray95")
  tcltk::tkpack(env$status_label, fill = "x", padx = 5, pady = 5)

  # Tree info frame
  info_frame <- tcltk::tkframe(main_frame, borderwidth = 1, relief = "sunken", background = "lightyellow")
  tcltk::tkpack(info_frame, fill = "x", padx = 5, pady = 5)

  env$tree_info_label <- tcltk::tklabel(info_frame, text = "All trees visible",
                                        anchor = "w", font = tcltk::tkfont.create(size = 9),
                                        background = "lightyellow")
  tcltk::tkpack(env$tree_info_label, fill = "x", padx = 5, pady = 5)

  # Separator
  tcltk::tkpack(tcltk::tkframe(main_frame, height = 2, borderwidth = 1, relief = "sunken"),
                fill = "x", pady = 10)

  # Buttons
  button_section_label <- tcltk::tklabel(main_frame, text = "Actions",
                                         font = tcltk::tkfont.create(size = 10, weight = "bold"))
  tcltk::tkpack(button_section_label, pady = c(0, 5))

  env$btn_pick <- tcltk::tkbutton(main_frame, text = "Select Tree",
                                  command = enable_tree_picking,
                                  width = 25, pady = 5)
  tcltk::tkpack(env$btn_pick, fill = "x", padx = 10, pady = 3)

  env$btn_extract <- tcltk::tkbutton(main_frame, text = "Extract Context",
                                     command = extract_and_display_context,
                                     width = 25, pady = 5, fg = "dark green")
  tcltk::tkpack(env$btn_extract, fill = "x", padx = 10, pady = 3)

  env$btn_export <- tcltk::tkbutton(main_frame, text = "Export Tree",
                                    command = export_selected_tree,
                                    width = 25, pady = 5, fg = "dark blue")
  tcltk::tkpack(env$btn_export, fill = "x", padx = 10, pady = 3)

  env$btn_restore <- tcltk::tkbutton(main_frame, text = "Restore View",
                                     command = restore_original,
                                     width = 25, pady = 5)
  tcltk::tkpack(env$btn_restore, fill = "x", padx = 10, pady = 3)

  # Separator
  tcltk::tkpack(tcltk::tkframe(main_frame, height = 2, borderwidth = 1, relief = "sunken"),
                fill = "x", pady = 10)

  reassign_section_label <- tcltk::tklabel(main_frame, text = "Point Reassignment",
                                           font = tcltk::tkfont.create(size = 10, weight = "bold"))
  tcltk::tkpack(reassign_section_label, pady = c(0, 5))

  env$btn_reassign_to_tree <- tcltk::tkbutton(main_frame, text = "Reassign Points to Tree",
                                              command = reassign_points_to_current_tree,
                                              width = 25, pady = 5, fg = "dark orange")
  tcltk::tkpack(env$btn_reassign_to_tree, fill = "x", padx = 10, pady = 3)

  env$btn_reassign_to_na <- tcltk::tkbutton(main_frame, text = "Reassign Points to NA",
                                            command = reassign_points_to_na,
                                            width = 25, pady = 5, fg = "dark red")
  tcltk::tkpack(env$btn_reassign_to_na, fill = "x", padx = 10, pady = 3)

  # Separator
  tcltk::tkpack(tcltk::tkframe(main_frame, height = 2, borderwidth = 1, relief = "sunken"),
                fill = "x", pady = 10)

  env$btn_exit <- tcltk::tkbutton(main_frame, text = "Exit",
                                  command = exit_application,
                                  width = 25, pady = 5, fg = "red")
  tcltk::tkpack(env$btn_exit, fill = "x", padx = 10, pady = 3)

  # Set initial button states
  tcltk::tkconfigure(env$btn_extract, state = "disabled")
  tcltk::tkconfigure(env$btn_export, state = "disabled")
  tcltk::tkconfigure(env$btn_restore, state = "disabled")
  tcltk::tkconfigure(env$btn_reassign_to_tree, state = "disabled")
  tcltk::tkconfigure(env$btn_reassign_to_na, state = "disabled")

  # Help text
  help_text <- "
Workflow:
1. 'Select Tree' -> Click point in 3D viewer
2. 'Extract Context' -> View tree & surroundings
3. 'Reassign Points' -> Draw polygon to modify
4. 'Export Tree' -> Save to file
5. 'Restore View' -> Show all trees
6. 'Exit' -> Close & return modified data
"
  help_label <- tcltk::tklabel(main_frame, text = help_text,
                               justify = "left", font = tcltk::tkfont.create(size = 8),
                               fg = "gray40")
  tcltk::tkpack(help_label, pady = c(5, 10))

  update_status("Ready", "black")

  # Force window to appear
  tcltk::tcl("update", "idletasks")

  # ---- Event loop (non-blocking) ----
  message("Arbor Studio Instance Module launched. Close the window to return to R console.")

  tcltk::tkwait.window(env$tt)

  # Use a simple check loop instead of tkwait.window
  # This returns control to R immediately
  return(env$las_orig)
}
