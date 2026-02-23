generate_cage = function(circles, params, header)
{
  ans = generate_cage_cpp(circles,  params$path_finder$decimation)
  cage$passage = 1000
  cage$hag = 0
  cage = suppressWarnings(lidR::LAS(cage, header))
  return(cage)
}

# r_generate_cage = function(circles, params, header)
# {
#   # Convert circle to points
#   res = params$path_finder$decimation*0.75
#   circle_points_list <- lapply(1:nrow(circles), function(i)
#   {
#     generate_circle_points(circles$X[i], circles$Y[i], circles$Z[i], circles$R[i], step = res)
#   })
#   circle_points = do.call(rbind, circle_points_list)
#
#   if (FALSE) qpoints3d(circle_points)
#
#   connectors <- generate_all_connectors(circles, step_z = res)
#
#   if (FALSE) qpoints3d(connectors)
#
#   cage = rbind(circle_points, connectors)
#   cage$passage = 1000
#   cage$hag = 0
#   cage = suppressWarnings(lidR::LAS(cage, header))
#
#   if (FALSE) qpoints3d(cage@data)
#
#   return(cage)
# }
#
# # Convert to points
# generate_circle_points <- function(x, y, z, r, step = 0.1)
# {
#   # Circumference
#   circumference <- 2 * pi * r
#   # Number of points for approximately every `step` meters
#   n_points <- ceiling(circumference / step)
#   # Angles for points
#   theta <- seq(0, 2*pi, length.out = n_points + 1)[-1]  # remove last point to avoid duplicate
#   # Generate points
#   data.frame(
#     X = x + r * cos(theta),
#     Y = y + r * sin(theta),
#     Z = rep(z, n_points)
#   )
# }
#
# # --- Generate random points on a disk (XY plane) ---
# generate_disk_points <- function(x, y, z, r, n = 1000)
# {
#   theta  <- stats::runif(n, 0, 2*pi)
#   radius <- sqrt(stats::runif(n)) * r
#   data.frame(
#     X = x + radius * cos(theta),
#     Y = y + radius * sin(theta),
#     Z = rep(z, n)
#   )
# }
#
# generate_disk_radii <- function(x, y, z, r, n = 8)
# {
#   angles <- seq(0, 2*pi, length.out = n + 1)[- (n + 1)]  # remove last to avoid duplicating 0
#
#   data.frame(
#     X = x + r * cos(angles),
#     Y = y + r * sin(angles),
#     Z = rep(z, n)
#   )
# }
#
#
# # --- Generate points for one group ---
# generate_connectors <- function(df_group, step_z = 0.05)
# {
#   df_group <- df_group[order(df_group$Z), ]
#   connectors <- list()
#
#   for (i in seq_len(nrow(df_group) - 1))
#   {
#     c1 <- df_group[i, ]
#     c2 <- df_group[i + 1, ]
#
#     # vertical interpolation sequence (every step_z)
#     z_seq <- seq(c1$Z, c2$Z, by = step_z)
#     if (utils::tail(z_seq, 1) != c2$Z)
#       z_seq <- c(z_seq, c2$Z)  # include top
#
#     t_seq <- (z_seq - c1$Z) / (c2$Z - c1$Z)
#
#     # interpolate center and radius
#     x_seq <- c1$X + t_seq * (c2$X - c1$X)
#     y_seq <- c1$Y + t_seq * (c2$Y - c1$Y)
#     r_seq <- c1$R + t_seq * (c2$R - c1$R)
#
#     # add 4-radii disks
#     for (j in seq_along(z_seq))
#     {
#       connectors[[length(connectors) + 1]] <- generate_disk_radii(x_seq[j], y_seq[j], z_seq[j], r_seq[j])
#     }
#   }
#
#   connectors = do.call(rbind, connectors)
#
#   return(connectors)
# }
#
# # Wrapper for all groups
# generate_all_connectors <- function(df, step_z = 0.05)
# {
#   groups <- unique(df$id)
#   connectors <- list()
#
#   for (g in seq_along(groups))
#   {
#     df_g <- df[df$id == groups[g], ]
#     connectors[[g]] <- generate_connectors(df_g, step_z)
#   }
#
#   connectors = do.call(rbind, connectors)
#
#   if (FALSE) qpoints3d(connectors)
#
#   return(connectors)
# }
