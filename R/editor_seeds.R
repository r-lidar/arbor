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
  setup_palette <- function(ids)
  {
    unique_ids <- unique(ids)
    colors_palette <- pastel.colors(length(unique_ids))
    names(colors_palette) <- as.character(unique_ids)
    return(colors_palette)
  }

  hud_id <- NULL
  update_hud <- function(text, color = "black")
  {
    if (!is.null(hud_id)) rgl::pop3d(id = hud_id)
    ids <- rgl::mtext3d(text, edge = "y+-", line = 2, at = NULL, pos = NA, col = color, cex = 1.5, font = 2)
    assign("hud_id", ids, envir = parent.frame())
  }

  colors_palette <- setup_palette(seeds$treeID)

  get_colors <- function(ids)
  {
    cols <- colors_palette[as.character(ids)]
    cols[is.na(cols)] <- "white"
    return(cols)
  }

  cat("Rendering background point cloud...\n")

  pc <- lidR::filter_poi(las, hag < 2)
  long_passages <- lidR::filter_poi(las, hag < 4, passage > 150)

  x <- plot_semantic(pc)
  lidR::plot(long_passages, add = x, pal = "gray", size = 2)

  seeds$X <- seeds$X - x[1]
  seeds$Y <- seeds$Y - x[2]

  # C. ADD LARGE, EASY-TO-FIND EXIT ZONE
  # Create a prominent EXIT plane at a fixed screen position
  bbox_x <- range(seeds$X)
  bbox_y <- range(seeds$Y)
  bbox_z <- range(seeds$Z)

  # Create multiple EXIT markers for visibility from any angle
  exit_size <- diff(bbox_x) * 0.05
  exit_positions <- data.frame(
    X = c(bbox_x[2] + exit_size, bbox_x[1] - exit_size,  mean(bbox_x), mean(bbox_x)),
    Y = c(mean(bbox_y), mean(bbox_y), bbox_y[2] + exit_size, bbox_y[1] - exit_size),
    Z = rep(bbox_z[2] + diff(bbox_z) * 0.1, 4)
  )

  # Large visible EXIT markers
  exit_marker_ids <- rgl::points3d(exit_positions$X, exit_positions$Y, exit_positions$Z, col = "red", size = 30)
  exit_text_ids <- rgl::text3d(exit_positions$X, exit_positions$Y, exit_positions$Z, texts = rep("EXIT", 4), col = "red", cex = 3, adj = c(0.5, -1))

  # Create EXIT connection plane for easier clicking
  exit_plane_x <- c(exit_positions$X, exit_positions$X[1])
  exit_plane_y <- c(exit_positions$Y, exit_positions$Y[1])
  exit_plane_z <- c(exit_positions$Z, exit_positions$Z[1])
  exit_plane_id <- rgl::points3d(exit_plane_x, exit_plane_y, exit_plane_z, col = "red", lwd = 5)

  # D. RENDER EDITABLE SEEDS
  seeds_obj_id <-rgl::points3d(seeds$X, seeds$Y, seeds$Z, col = get_colors(seeds$treeID), size = 7)

  rgl::axes3d()
  rgl::title3d(main = "TreeID Editor - Middle-Click to Edit or EXIT (red markers)", xlab = "X", ylab = "Y", zlab = "Z")

  # E. INSTRUCTIONS
  cat("\n======================================================\n")
  cat("WORKFLOW:\n")
  cat("  1. Middle-click KEEP seed (blue highlight)\n")
  cat("  2. Middle-click CHANGE seeds (red highlights)\n")
  cat("  3. Middle-click CHANGE seed AGAIN to unselect it\n")
  cat("  4. Middle-click KEEP seed AGAIN to validate merge\n")
  cat("\n")
  cat("EXIT OPTIONS:\n")
  cat("  • Middle-click any RED EXIT marker (4 around scene)\n")
  cat("  • Close the RGL window\n")
  cat("  • Press ESC during identify3d() to cancel and exit\n")
  cat("======================================================\n\n")

  # State machine
  state <- "SELECT_KEEP"
  keep_seed <- NULL
  change_seeds <- list()

  # Visual tracking
  keep_marker <- NULL
  keep_text <- NULL
  change_markers <- list()
  change_texts <- list()
  change_lines <- list()

  # Combine seeds + exit points
  all_points_x <- c(seeds$X, exit_positions$X)
  all_points_y <- c(seeds$Y, exit_positions$Y)
  all_points_z <- c(seeds$Z, exit_positions$Z)

  n_seeds <- nrow(seeds)
  exit_indices <- (n_seeds + 1):(n_seeds + nrow(exit_positions))

  # Interactive loop
  repeat
  {
    # === EXIT CONDITION 1: Window closed ===
    # (not actually working)
    if (length(rgl.dev.list()) == 0 || rgl.cur() == 0)
    {
      cat("\n✓ Exiting editor (window closed)\n")
      return(seeds)
    }

    # Update HUD based on state
    if (state == "SELECT_KEEP")
    {
      update_hud("Step 1: Click KEEP seed | Click EXIT to quit", "blue")
    }
    else if (state == "SELECT_CHANGE")
    {
      n_changes <- length(change_seeds)
      update_hud(sprintf("KEEP: %s | Click CHANGE seeds or KEEP to merge (%d) | EXIT to quit", keep_seed$treeID, n_changes), "orange")
    }

    # Wait for user click
    sel <- tryCatch(
    {
      rgl::identify3d(all_points_x, all_points_y, all_points_z, n = 1, plot = FALSE, buttons = "middle", tolerance = 40)
    },
    error = function(e)
    {
      # User pressed ESC or error occurred
      # (not actually working)
      cat("\n✓ Exiting editor (ESC or error)\n")
      return(NULL)
    })

    # === EXIT CONDITION 2: ESC pressed or error ===
    if (is.null(sel))
    {
      break
    }

    if (length(sel) == 0)
    {
      next  # No selection, continue
    }

    # === EXIT CONDITION 3: EXIT marker clicked ===
    if (sel %in% exit_indices)
    {
      update_hud("Exiting editor...", "green")
      cat("\n✓ Exiting editor (EXIT marker clicked)\n")
      Sys.sleep(0.5)
      break
    }

    # Validate selection is a seed
    if (sel < 1 || sel > n_seeds)
    {
      cat("[WARN] Invalid selection, try again\n")
      next
    }

    # Get selected seed
    selected_seed <- seeds@data[sel, ]
    selected_id <- selected_seed$treeID

    cat(sprintf("\nClicked seed index: %d, treeID: %s\n", sel, selected_id))
    cat(sprintf("Current state: %s\n", state))

    # STATE MACHINE LOGIC
    if (state == "SELECT_KEEP")
    {
      keep_seed <- selected_seed
      state <- "SELECT_CHANGE"

      # Visual feedback
      keep_marker <- rgl::points3d(keep_seed$X, keep_seed$Y, keep_seed$Z, col = "blue", size = 15)
      keep_text <- rgl::text3d(keep_seed$X, keep_seed$Y, keep_seed$Z, texts = paste("KEEP:", keep_seed$treeID), col = "blue", cex = 1.5, adj = c(0, -1.5))
      cat(sprintf("\n>>> KEEP seed selected: ID %s\n", keep_seed$treeID))
    }
    else if (state == "SELECT_CHANGE")
    {
      # Check if clicked KEEP seed again → VALIDATE
      if (selected_id == keep_seed$treeID)
      {
        if (length(change_seeds) == 0)
        {
          update_hud("No CHANGE seeds selected. Select some first!", "red")
          cat("⚠ No CHANGE seeds to merge\n")
          Sys.sleep(1.5)
          next
        }

        # PERFORM MERGE
        change_ids <- sapply(change_seeds, function(s) s$treeID)

        update_hud(sprintf("Merging %d seeds into %s...", length(change_ids), keep_seed$treeID), "yellow")

        cat("\n=======================================================\n")
        cat(sprintf("MERGING: %s → %s\n",  paste(change_ids, collapse = ", "), keep_seed$treeID))
        cat("=======================================================\n")

        # Update data
        for (change_id in change_ids)
        {
          seeds$treeID[seeds$treeID == change_id] <- keep_seed$treeID
        }

        # Refresh visuals
        pop3d(id = seeds_obj_id)
        colors_palette <- setup_palette(seeds$treeID)
        seeds_obj_id <- rgl::points3d(seeds$X, seeds$Y, seeds$Z, col = get_colors(seeds$treeID), size = 7)

        update_hud(sprintf("✓ Merged %d seeds successfully!", length(change_ids)), "green")
        cat(sprintf("✓ Successfully merged %d IDs\n\n", length(change_ids)))
        Sys.sleep(1.5)

        # CLEANUP & RESET STATE
        pop3d(id = keep_marker)
        pop3d(id = keep_text)
        for (m in change_markers) rgl::pop3d(id = m)
        for (t in change_texts)   rgl::pop3d(id = t)
        for (l in change_lines)   rgl::pop3d(id = l)

        keep_seed <- NULL
        change_seeds <- list()
        change_markers <- list()
        change_texts <- list()
        change_lines <- list()
        state <- "SELECT_KEEP"

      }
      else
      {
        # CHECK IF THIS CHANGE SEED IS ALREADY SELECTED → UNSELECT IT
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
          update_hud(sprintf("Unselecting CHANGE seed: %s", selected_id), "purple")
          cat(sprintf(">>> CHANGE seed unselected: ID %s\n", selected_id))

          # Remove visual elements
          rgl::pop3d(id = change_markers[[existing_index]])
          rgl::pop3d(id = change_texts[[existing_index]])
          rgl::pop3d(id = change_lines[[existing_index]])

          # Remove from lists
          change_seeds   <- change_seeds[-existing_index]
          change_markers <- change_markers[-existing_index]
          change_texts   <- change_texts[-existing_index]
          change_lines   <- change_lines[-existing_index]

          Sys.sleep(0.5)
        }
        else
        {
          # ADD NEW CHANGE SEED
          change_seeds[[length(change_seeds) + 1]] <- selected_seed

          # Visual feedback
          marker   <- rgl::points3d(selected_seed$X, selected_seed$Y, selected_seed$Z, col = "red", size = 15)
          text_obj <- rgl::text3d(selected_seed$X, selected_seed$Y, selected_seed$Z, texts = paste("CHANGE:", selected_id), col = "red", cex = 1.5, adj = c(0, -1.5))
          line_obj <- rgl::segments3d(c(keep_seed$X, selected_seed$X), c(keep_seed$Y, selected_seed$Y), c(keep_seed$Z, selected_seed$Z), col = "red", lwd = 3)

          change_markers[[length(change_markers) + 1]] <- marker
          change_texts[[length(change_texts) + 1]] <- text_obj
          change_lines[[length(change_lines) + 1]] <- line_obj

          cat(sprintf(">>> CHANGE seed #%d added: ID %s\n", length(change_seeds), selected_id))
        }
      }
    }
  }

  # Ensure window closes gracefully
  if (length(rgl.dev.list()) > 0)
  {
    tryCatch(rgl::close3d(), error = function(e) invisible(NULL))
  }

  seeds$X <- seeds$X + x[1]
  seeds$Y <- seeds$Y + x[2]
  return(seeds)
}
