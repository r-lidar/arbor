# @file studio_instance.R
# Project: Arbor
#
# Copyright (C) 2026 Jean-Romain Roussel (r-lidar) <info @ r-lidar.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

#' @export
#' @rdname arbor-studio
#' @noRd
arbor_studio_instance <- function(las, ...)
{
  . <- R <- G <- B <- treeID <- NULL

  attribute = "treeID"

  las = colorize_trees(las, FALSE)
  col_table = las@data[, list(R = R[1], G = G[1], B = B[1]), by = treeID]
  data.table::setkey(col_table, treeID)

  # ---- Input validation ----
  if (!inherits(las, "LAS"))  stop("'las' must be a LAS object", call. = FALSE)
  if (!attribute %in% names(las@data)) stop(sprintf("Attribute '%s' not found in LAS object", attribute), call. = FALSE)
  if (!interactive()) stop("arbor_studio_instance() requires an interactive R session", call. = FALSE)

  # ---- Setup environment ----
  # Use an environment to avoid global assignment and properly scope variables
  env <- new.env(parent = emptyenv())

  # Initialize state
  env$las_orig <- las
  env$las_orig@data$pointID <- seq_len(lidR::npoints(las))
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
  env$col_table <- col_table

  # ---- Helper functions ----

  redraw_scene <- function()
  {
    rgl::clear3d()
    lidR::plot(env$current_las, add = c(env$center_x, env$center_y),  bg = "black", color = "RGB")
  }

  update_status <- function(status_text, color = "black")
  {
    env$current_state <- status_text
    if (!is.null(env$status_label))
    {
        tcltk::tkconfigure(env$status_label, text = paste("Status:", status_text), fg = color)
        tcltk::tcl("update", "idletasks")
    }
  }

  # ---- Exit application ----

  exit_application <- function()
  {
    # Atomic check-and-set to prevent re-entry
    if (!isTRUE(env$app_running)) return()

    env$app_running <- FALSE  # Set IMMEDIATELY before any other operations

    update_status("Closing application...")

    # Cancel any pending tcl callbacks
    tryCatch({
      tcltk::tcl("after", "cancel", check_rgl_window)
    }, error = function(e) invisible(NULL))

    # Set return value
    if (isTRUE(env$modified)) {
      env$return_value <- env$las_orig
    } else {
      env$return_value <- NULL
    }

    # Destroy window
    tryCatch({
      if (!is.null(env$tt) && tcltk::tclvalue(tcltk::tkwinfo("exists", env$tt)) == "1") {
        tcltk::tkdestroy(env$tt)
        env$tt <- NULL
      }
    }, error = function(e) invisible(NULL))
  }

  # ---- Extract context ----

  extract_and_display_context <- function()
  {
    if (env$mouse_active) {
      update_status("Please wait for current operation to finish", "red")
      return()
    }

    if (is.null(env$selected_treeID)) {
      update_status("No tree selected. Please select a tree first.", "red")
      return()
    }

    update_status("Extracting tree context...")

    tryCatch({
      env$current_las <- extract_tree_context(env$las_orig, env$selected_treeID)

      i = which(env$current_las[[env$attribute]] ==  env$selected_treeID)
      env$current_las@data[["R"]][i] = 185*255
      env$current_las@data[["G"]][i] = 0
      env$current_las@data[["B"]][i] = 0

      tree_info <- sprintf("Tree ID: %s | Height: %.2fm", env$selected_treeID, max(env$current_las$Z) - min(env$current_las$Z))
      tcltk::tkconfigure(env$tree_info_label, text = tree_info)

      redraw_scene()

      # Update button states
      tcltk::tkconfigure(env$btn_pick, state = "disabled")
      tcltk::tkconfigure(env$btn_extract, state = "disabled")
      tcltk::tkconfigure(env$btn_export, state = "normal")
      tcltk::tkconfigure(env$btn_restore, state = "normal")
      tcltk::tkconfigure(env$btn_reassign_to_tree, state = "normal")
      tcltk::tkconfigure(env$btn_reassign_to_tree_click, state = "normal")
      tcltk::tkconfigure(env$btn_reassign_to_na, state = "normal")

      update_status(paste("Context extracted for Tree ID", env$selected_treeID))

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

  restore_original <- function()
  {
    if (env$mouse_active) {
      update_status("Cannot restore while operation is active", "red")
      return()
    }

    update_status("Restoring original view...")

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
    tcltk::tkconfigure(env$btn_reassign_to_tree_click, state = "disabled")
    tcltk::tkconfigure(env$btn_reassign_to_na, state = "disabled")

    update_status("Original view restored")

    tcltk::tcl("after", 1500, function() {
      if (env$current_state == "Original view restored") {
        update_status("Ready", "black")
      }
    })
  }

  # ---- Point picking ----

  enable_tree_picking <- function()
  {
    if (env$mouse_active) {
      update_status("Tree picking already active")
      return()
    }

    env$mouse_active <- TRUE
    update_status("Click a point in the 3D viewer to select a tree...")

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

    if (nrow(ids) > 0)
    {
      selected_idx <- rgl::identify3d(x = env$las_orig$X-env$center_x, y = env$las_orig$Y-env$center_y, z = env$las_orig$Z, n = 1, plot = FALSE, buttons = "middle", tolerance = 40)

      if (length(selected_idx) > 0 && selected_idx > 0)
      {
        treeID <- env$las_orig$treeID[selected_idx]
        env$selected_treeID <- treeID

        rgl::texts3d(
          env$las_orig$X[selected_idx] - env$center_x,
          env$las_orig$Y[selected_idx] - env$center_y,
          env$las_orig$Z[selected_idx],
          texts = treeID,
          color = "white",
          cex = 1.2
        )

        update_status(paste("Tree ID", treeID, "selected - Click 'Extract Context' to view"))

        tree_info <- sprintf("Tree ID: %s selected [Context not extracted yet]", treeID)
        tcltk::tkconfigure(env$tree_info_label, text = tree_info)

        tcltk::tkconfigure(env$btn_pick, state = "normal")
        tcltk::tkconfigure(env$btn_extract, state = "normal")
        tcltk::tkconfigure(env$btn_exit, state = "normal")

        tcltk::tcl("after", 100, function()
        {
          tcltk::tkfocus(env$btn_extract)
        })

      }
      else
      {
        update_status("No point selected", "orange")
        env$selected_treeID <- NULL
        tcltk::tkconfigure(env$btn_pick, state = "normal")
        tcltk::tkconfigure(env$btn_exit, state = "normal")
      }
    }
    else
    {
      update_status("No shapes found in scene", "red")
      tcltk::tkconfigure(env$btn_pick, state = "normal")
      tcltk::tkconfigure(env$btn_exit, state = "normal")
    }

    env$mouse_active <- FALSE
  }

  # ---- Reassign points to tree ----

  reassign_points_to_current_tree <- function(mode = "rectangle")
  {
    mode = match.arg(mode, c("click", "rectangle"))

    if (env$mouse_active) {
      update_status("Please wait for current operation to finish", "blue")
      return()
    }

    if (is.null(env$selected_treeID)) {
      update_status("No tree selected. Cannot reassign.", "red")
      return()
    }

    if (is.null(env$current_las) || nrow(env$current_las@data) == 0 || identical(env$current_las, env$las_orig)) {
      update_status("No context extracted. Please extract context first.", "red")
      return()
    }

    env$mouse_active <- TRUE

    if (mode == "rectangle") {
      update_status(paste("Draw a rectangle in 3D viewer to select points for Tree ID", env$selected_treeID))
    } else {
      update_status(paste("Click a point to reassign all points with its ID to Tree ID", env$selected_treeID))
    }

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

    tryCatch(
    {
      if (mode == "rectangle")
      {
        update_status("Draw a rectancle to select points", "blue")

        selection_fn <- rgl::select3d()
        selected <- selection_fn(env$current_las$X - env$center_x, env$current_las$Y - env$center_y, env$current_las$Z)

        if (sum(selected) > 0)
        {
          col = env$col_table[.(env$selected_treeID)]
          R = col$R
          G = col$G
          B = col$B
          env$current_las@data[[env$attribute]][selected] <- env$selected_treeID
          env$current_las@data[["R"]][selected] <- 185*255
          env$current_las@data[["G"]][selected] <- 0
          env$current_las@data[["B"]][selected] <- 0

          ids <- env$current_las$pointID[selected]
          env$las_orig@data[[env$attribute]][ids] <- env$selected_treeID
          env$las_orig@data[["R"]][ids] <- R
          env$las_orig@data[["G"]][ids] <- G
          env$las_orig@data[["B"]][ids] <- B
          env$modified <- TRUE

          redraw_scene()

          update_status(paste(sum(selected), "points reassigned to Tree ID",  env$selected_treeID), "green")

          tree_info <- sprintf("Tree ID: %s | Points: %d | Height: %.2fm",
                               env$selected_treeID,
                               sum(env$current_las[[env$attribute]] == env$selected_treeID, na.rm = TRUE),
                               max(env$current_las$Z, na.rm = TRUE) - min(env$current_las$Z, na.rm = TRUE))
          tcltk::tkconfigure(env$tree_info_label, text = tree_info)
        }
        else
        {
          update_status("No points selected", "orange")
        }

      }
      else if (mode == "click")
      {
        update_status("Click on a point to select its ID...", "blue")

        clicked <- rgl::identify3d(x = env$current_las$X - env$center_x, y = env$current_las$Y - env$center_y, z = env$current_las$Z, n = 1, plot = FALSE, buttons = "middle", tolerance = 40)

        if (length(clicked) > 0 && clicked[1] > 0)
        {
          # Get the ID of the clicked point
          clicked_id <- env$current_las@data[[env$attribute]][clicked[1]]

          if (is.na(clicked_id))
          {
            update_status("Clicked point has no ID (NA). Cannot reassign.", "orange")
          }
          else
          {
            # Find all points with the same ID in current_las
            points_to_reassign <- env$current_las@data[[env$attribute]] == clicked_id & !is.na(env$current_las@data[[env$attribute]])
            n_points <- sum(points_to_reassign)

            if (n_points > 0)
            {
              # Reassign in current_las
              col = env$col_table[.(env$selected_treeID)]
              R = col$R
              G = col$G
              B = col$B
              env$current_las@data[[env$attribute]][points_to_reassign] <- env$selected_treeID
              env$current_las@data[["R"]][points_to_reassign] <- 185*255
              env$current_las@data[["G"]][points_to_reassign] <- 0
              env$current_las@data[["B"]][points_to_reassign] <- 0

              # Reassign in original las
              ids <- env$current_las$pointID[points_to_reassign]
              env$las_orig@data[[env$attribute]][ids] <- env$selected_treeID
              env$las_orig@data[["R"]][ids] <- R
              env$las_orig@data[["G"]][ids] <- G
              env$las_orig@data[["B"]][ids] <- B
              env$modified <- TRUE

              redraw_scene()

              update_status(paste(n_points, "points with ID", clicked_id, "reassigned to Tree ID", env$selected_treeID), "green")

              tree_info <- sprintf("Tree ID: %s | Points: %d | Height: %.2fm",
                                   env$selected_treeID,
                                   sum(env$current_las[[env$attribute]] == env$selected_treeID, na.rm = TRUE),
                                   max(env$current_las$Z, na.rm = TRUE) - min(env$current_las$Z, na.rm = TRUE))
              tcltk::tkconfigure(env$tree_info_label, text = tree_info)
            }
            else
            {
              update_status("No points found with the selected ID", "orange")
            }
          }
        }
        else
        {
          update_status("No point clicked", "orange")
        }
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

    tcltk::tcl("after", 2000, function()
    {
      if (grepl("reassigned", env$current_state))
      {
        update_status("Ready - Context extracted", "green")
      }
    })
  }


  # ---- Reassign points to NA ----

  reassign_points_to_na <- function()
  {
    if (env$mouse_active)
    {
      update_status("Please wait for current operation to finish")
      return()
    }

    if (is.null(env$current_las) || nrow(env$current_las@data) == 0 || identical(env$current_las, env$las_orig))
    {
      update_status("No context extracted. Please extract context first.")
      return()
    }

    env$mouse_active <- TRUE
    update_status("Draw a polygon in 3D viewer to select points for NA assignment")

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

    tryCatch(
    {
      selection_fn <- rgl::select3d()
      selected <- selection_fn(env$current_las$X - env$center_x, env$current_las$Y - env$center_y, env$current_las$Z)
      selected[env$current_las[[env$attribute]] != env$selected_treeID] = FALSE

      if (sum(selected) > 0)
      {
        env$current_las@data[[env$attribute]][selected] <- NA
        env$current_las@data[["R"]][selected] <- 38250
        env$current_las@data[["G"]][selected] <- 38250
        env$current_las@data[["B"]][selected] <- 38250

        ids <- env$current_las$pointID[selected]
        env$las_orig@data[[env$attribute]][ids] <- NA
        env$las_orig@data[["R"]][ids] <- 38250
        env$las_orig@data[["G"]][ids] <- 38250
        env$las_orig@data[["B"]][ids] <- 38250
        env$modified <- TRUE

        redraw_scene()

        update_status(paste(sum(selected), "points reassigned to NA"), "green")

        tree_info <- sprintf("Tree ID: %s | Points: %d | Height: %.2fm",
                             ifelse(is.null(env$selected_treeID), "NA", as.character(env$selected_treeID)),
                             nrow(env$current_las@data),
                             max(env$current_las$Z, na.rm = TRUE) - min(env$current_las$Z, na.rm = TRUE))
        tcltk::tkconfigure(env$tree_info_label, text = tree_info)

      }
      else
      {
        update_status("No points selected", "orange")
      }

    },
    error = function(e)
    {
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

    tcltk::tcl("after", 2000, function()
    {
      if (grepl("reassigned", env$current_state))
      {
        update_status("Ready - Context extracted", "green")
      }
    })
  }

  # ---- Export tree ----

  export_selected_tree <- function()
  {
    if (env$mouse_active)
    {
      update_status("Please wait for current operation to finish")
      return()
    }

    if (is.null(env$current_las) || nrow(env$current_las@data) == 0 || identical(env$current_las, env$las_orig))
    {
      update_status("No context extracted. Please extract context first.")
      return()
    }

    update_status("Exporting tree...")

    tryCatch(
    {
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

  title_label <- tcltk::tklabel(main_frame, text = "Tree Instance Viewer", font = tcltk::tkfont.create(size = 12, weight = "bold"))
  tcltk::tkpack(title_label, pady = c(5, 10))

  # Status frame
  status_frame <- tcltk::tkframe(main_frame, borderwidth = 1, relief = "sunken")
  tcltk::tkpack(status_frame, fill = "x", padx = 5, pady = 5)

  env$status_label <- tcltk::tklabel(status_frame, text = "Status: Ready", anchor = "w", font = tcltk::tkfont.create(size = 10))
  tcltk::tkpack(env$status_label, fill = "x", padx = 5, pady = 5)

  # Tree info frame
  info_frame <- tcltk::tkframe(main_frame, borderwidth = 1, relief = "sunken")
  tcltk::tkpack(info_frame, fill = "x", padx = 5, pady = 5)

  env$tree_info_label <- tcltk::tklabel(info_frame, text = "All trees visible", anchor = "w", font = tcltk::tkfont.create(size = 9))
  tcltk::tkpack(env$tree_info_label, fill = "x", padx = 5, pady = 5)

  # Separator
  tcltk::tkpack(tcltk::tkframe(main_frame, height = 2, borderwidth = 1, relief = "sunken"), fill = "x", pady = 10)

  # Buttons
  button_section_label <- tcltk::tklabel(main_frame, text = "Actions", font = tcltk::tkfont.create(size = 10, weight = "bold"))
  tcltk::tkpack(button_section_label, pady = c(0, 5))

  env$btn_pick <- tcltk::tkbutton(main_frame, text = "Select Tree", command = enable_tree_picking, pady = 5)
  tcltk::tkpack(env$btn_pick, fill = "x", padx = 10, pady = 3)

  env$btn_extract <- tcltk::tkbutton(main_frame, text = "Extract Context", command = extract_and_display_context, pady = 5)
  tcltk::tkpack(env$btn_extract, fill = "x", padx = 10, pady = 3)

  env$btn_export <- tcltk::tkbutton(main_frame, text = "Export Tree", command = export_selected_tree, pady = 5)
  tcltk::tkpack(env$btn_export, fill = "x", padx = 10, pady = 3)

  env$btn_restore <- tcltk::tkbutton(main_frame, text = "Restore View", command = restore_original, pady = 5)
  tcltk::tkpack(env$btn_restore, fill = "x", padx = 10, pady = 3)

  # Separator
  tcltk::tkpack(tcltk::tkframe(main_frame, height = 2, borderwidth = 1, relief = "sunken"), fill = "x", pady = 10)

  reassign_section_label <- tcltk::tklabel(main_frame, text = "Point Reassignment", font = tcltk::tkfont.create(size = 10, weight = "bold"))
  tcltk::tkpack(reassign_section_label, pady = c(0, 5))

  env$btn_reassign_to_tree <- tcltk::tkbutton(main_frame, text = "Reassign Points to Tree", command = function() reassign_points_to_current_tree(mode = "rectangle"), pady = 5)
  tcltk::tkpack(env$btn_reassign_to_tree, fill = "x", padx = 10, pady = 3)

  env$btn_reassign_to_tree_click <- tcltk::tkbutton(main_frame, text = "Reassign by Click", command = function() reassign_points_to_current_tree(mode = "click"))
  tcltk::tkpack(env$btn_reassign_to_tree_click, fill = "x", padx = 10, pady = 3)

  env$btn_reassign_to_na <- tcltk::tkbutton(main_frame, text = "Reassign Points to NA", command = reassign_points_to_na, pady = 5)
  tcltk::tkpack(env$btn_reassign_to_na, fill = "x", padx = 10, pady = 3)

  # Separator
  tcltk::tkpack(tcltk::tkframe(main_frame, height = 2, borderwidth = 1, relief = "sunken"), fill = "x", pady = 10)

  env$btn_exit <- tcltk::tkbutton(main_frame, text = "Exit", command = exit_application, pady = 5, fg = "red")
  tcltk::tkpack(env$btn_exit, fill = "x", padx = 10, pady = 3)

  # Set initial button states
  tcltk::tkconfigure(env$btn_extract, state = "disabled")
  tcltk::tkconfigure(env$btn_export, state = "disabled")
  tcltk::tkconfigure(env$btn_restore, state = "disabled")
  tcltk::tkconfigure(env$btn_reassign_to_tree, state = "disabled")
  tcltk::tkconfigure(env$btn_reassign_to_tree_click, state = "disabled")
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
  help_label <- tcltk::tklabel(main_frame, text = help_text, justify = "left", font = tcltk::tkfont.create(size = 8), fg = "gray40")
  tcltk::tkpack(help_label, pady = c(5, 10))

  update_status("Ready", "black")

  # Monitor RGL window
  check_rgl_window <- function()
  {
    # Check flag FIRST - don't do anything if exiting
    if (!isTRUE(env$app_running)) return()

    # Check if RGL window still exists
    rgl_exists <- tryCatch({
      length(rgl::rgl.dev.list()) > 0 && rgl::rgl.cur() > 0
    }, error = function(e) FALSE)

    if (!rgl_exists && isTRUE(env$app_running)) {
      exit_application()
      return()
    }

    # Schedule next check only if still running
    if (isTRUE(env$app_running)) {
      tcltk::tcl("after", 500, check_rgl_window)
    }
  }
  check_rgl_window()

  # Force window to appear
  tcltk::tcl("update", "idletasks")

  # Wait for window to close
  while (env$app_running)
  {
    result <- tryCatch(
      {
        tcltk::tcl("update")  # Process tcltk events
        TRUE
      },
      error = function(e)
      {
        # Any error means we should exit
        FALSE
      })

    if (!result) {
      env$app_running <- FALSE
      break
    }

    Sys.sleep(0.05)  # Small delay to prevent CPU spinning
  }

  # Close RGL window
  if (length(rgl::rgl.dev.list()) > 0)
  {
    tryCatch(rgl::close3d(), error = function(e) invisible(NULL))
  }

  # Use a simple check loop instead of tkwait.window
  # This returns control to R immediately
  return(env$las_orig)
}
