library(lidR)
library(arbor)

files = list.files("/media/jr/easystore/r-lidar/Projets/2025/MNRF-MLS/wd/las/Sta/output/its/", full.names = TRUE)
files = list.files("/home/jr/Documents/Bastien/Individual trees/SainteVeronique/las", full.names = TRUE)
files = list.files("/home/jr/Documents/r-lidar/clients/fsinvestor/Zambia/JasonFarm/its/", full.names = TRUE, pattern = ".las")
files = list.files("/home/jr/Documents/r-lidar/clients/fsinvestor/Indonesia/Walk1Area2/ITS/", full.names = TRUE, pattern = ".las")
files = list.files("~/Téléchargements/QSM St Vero/Pall_registered/", full.names = TRUE, pattern = ".las")

file = "~/Téléchargements/QSM St Vero/Pall_registered/33BOJ.las"
file = "~/Téléchargements/QSM St Vero/Pall_registered/84ERS.las" # Vstem 2.3 m³
file = "~/Téléchargements/QSM St Vero/Pall_registered/34BOJ.las" # Error branches plus grosses que parents
file = "~/Téléchargements/QSM St Vero/Pall_registered/36BOJ.las" # Vstem 2.6 m³  Error branches plus grosses que parents
file = "~/Téléchargements/QSM St Vero/Pall_registered/45BOJ.las" # Vstem 2.5 m³

file = "/home/jr/Documents/r-lidar/clients/fsinvestor/Rwanda/Eucalyptus/Stoneaucd6_output/ITS/tree_65.las"
file = "/home/jr/Documents/r-lidar/clients/fsinvestor/Rwanda/Eucalyptus/Stoneaucd6_output/ITS/tree_130.las"
file = "/home/jr/Documents/r-lidar/clients/fsinvestor/Rwanda/Eucalyptus/Stoneaucd6_output/ITS/tree_177.las"

# Buttress
file = "/home/jr/Documents/r-lidar inc/arbor/Tree bank/bad-dbh/DUC0001-02_44.las"

# Big tree non circular
file = "/home/jr/Documents/r-lidar inc/arbor/Tree bank/non-circular/tree_1114.las"
file = "/home/jr/Documents/r-lidar inc/arbor/Tree bank/non-circular/tree_54.las"
file = "/home/jr/Documents/r-lidar inc/arbor/Tree bank/non-circular/tree_681.las"
file = "/home/jr/Documents/r-lidar inc/arbor/Tree bank/non-circular/tree_705.las"
file = "/home/jr/Documents/r-lidar inc/arbor/Tree bank/non-circular/tree_848.las"
file = "/home/jr/Documents/r-lidar inc/arbor/Tree bank/non-circular/tree_939.las"

file = "/home/jr/Documents/r-lidar inc/arbor/Tree bank/non-circular/wytham_15240.las"
file = "/home/jr/Documents/r-lidar inc/arbor/Tree bank/non-circular/wytham_15408.las"

# ==== Single QSM ======

# QSM for a random tree
i = sample(seq_along(files), 1)
file <- files[i]
files = files[-i]
file
tree <- lidR::readLAS(file)
tree@data$treeID = as.integer(tree@data$treeID)

params= default_arbor_parameters
qsm = qsm(tree, params)
stem = qsm_stem(qsm)
merch = qsm_merchantable(qsm, 0.045)
merch_stem  = qsm_stem(qsm) |> qsm_merchantable(0.045)
plot_qsm(qsm)
plot_qsm(stem)
plot_qsm(merch)
plot_qsm(merch_stem)

u = microbenchmark::microbenchmark(qsm_dbh(qsm), qsm_dbh_cpp(qsm), times = 50)
ggplot2::autoplot(u)

x = plot_semantic(tree)
plot_qsm(qsm, add = x, color = "branch_order", cylinder = T)

dbh = qsm_dbh_cpp(qsm)
x = plot_semantic(tree) |> add_dbh3d(dbh, lwd = 3)
plot(qsm, add = x, pal = "chocolate4")


tr = qsm[axis_ID == 1]
m = lm(radius ~ 0+subtree_length, data = tr[tr$measure == TRUE,])
plot(tr$subtree_length, tr$theoric_radius, col = "gray50")
points(tr$subtree_length, tr$radius, col = "blue")
abline(m)

distance_to_root = function(dt)
{
  # Build adjacency (children of each parent)
  children_list <- split(dt$cyl_ID, dt$parent_ID)

  # Initialize vector
  dist_root <- numeric(nrow(dt))
  names(dist_root) <- dt$cyl_ID

  # Recursive DFS function
  accumulate_dist <- function(node, current_dist) {
    dist_root[[as.character(node)]] <<- current_dist
    children <- children_list[[as.character(node)]]
    if (!is.null(children)) {
      for (child in children) {
        cyl_len <- dt[cyl_ID == child, cyl_length]
        accumulate_dist(child, current_dist + cyl_len)
      }
    }
  }

  # Run from the root (parent_ID == 0)
  root_id <- dt[parent_ID == 0, cyl_ID]
  accumulate_dist(root_id, 0)

  # Attach to table
  dt[, distance_to_root := dist_root[as.character(cyl_ID)]]
}



x = plot_semantic(tree)
passage <- lidR::filter_poi(tree, passage > 1)
sampled = lidR::decimate_points(tree, random_per_voxel(0.04))

data = data.table::copy(tree)

z = data$Z
z = z - min(z)
lay = lidR:::round_any(z, 0.05)
data@data$lay = lay
data@data[, Z := Z + lay*2, by = lay]
plot(data)
data = connected_components(data, 0.05, 1, connectivity = 26, "treeID")
plot(data, color = "treeID")

skel = data@data[, .(X = median(X), Y = median(Y), Z = median(Z), lay = min(lay)), by = treeID]
skel$treeID = NULL
skel[, Z := Z - lay*2, by = lay]
skel = LAS(skel)

x= plot(tree)
plot(skel, pal = "red", size = 4, add=x)

gr = compute_point_network(skel, k = 25, max_gap = 2, power = 1, downward = T)


library(igraph)
g <- graph_from_data_frame(gr, directed = FALSE)
mst <- mst(g, weights = E(g)$cost)
mst = as_data_frame(mst, what = "edges")
mst$from = as.integer(mst$from)
mst$to = as.integer(mst$to)

startX = skel$X[mst$from]
startY = skel$Y[mst$from]
startZ = skel$Z[mst$from]
endX = skel$X[mst$to]
endY = skel$Y[mst$to]
endZ = skel$Z[mst$to]
pid1 = mst$from
pid2 = mst$to
mst = data.table::data.table(startX, startY, startZ,  endX,  endY,  endZ, pid1, pid2)
plot_qsm(mst)
mst = qsm_topology(mst)
mst = qsm_architecture(mst)
plot_qsm(mst, color = "subtree_length")
lidR:::.pan3d(2)

len = mst[, max(subtree_length), by = axis_ID]
valid_axis = len[V1 > 2]$axis_ID

u = mst[axis_ID %in% valid_axis]

u = distance_to_root(u)

x = plot_semantic(tree)
plot_qsm(u, add = x, color = "distance_to_root")

u$radius = 0
v = qsm_simplify(u, 0.1)
v$radius = NULL

x = plot_semantic(tree)
plot_qsm(v, add = x, color = "subtree_length")

lidR:::.pan3d(2)
