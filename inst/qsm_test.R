library(lidR)
library(arbor)

files = list.files("/media/jr/easystore/r-lidar/Projets/2025/MNRF-MLS/wd/las/Sta/output/its/", full.names = TRUE)
files = list.files("/home/jr/Documents/Bastien/Individual trees/SainteVeronique/las", full.names = TRUE)
files = list.files("/home/jr/Documents/r-lidar/clients/fsinvestor/Zambia/JasonFarm/its/", full.names = TRUE, pattern = ".las")
files = list.files("/home/jr/Documents/r-lidar/clients/fsinvestor/Indonesia/Walk1Area2/ITS/", full.names = TRUE, pattern = ".las")
files = list.files("~/Téléchargements/QSM St Vero/Pall_registered/", full.names = TRUE, pattern = ".las")


file = "~/Téléchargements/QSM St Vero/Pall_registered/33BOJ.las"
file = "~/Téléchargements/QSM St Vero/Pall_registered/84ERS.las" # Vstem 2.3 m³
file = "~/Téléchargements/QSM St Vero/Pall_registered/34BOJ.las"
file = "~/Téléchargements/QSM St Vero/Pall_registered/36BOJ.las" # Vstem 2.6 m³
file = "~/Téléchargements/QSM St Vero/Pall_registered/45BOJ.las" # Vstem 2.5 m³ !!

file = "/home/jr/Documents/r-lidar/clients/fsinvestor/Rwanda/Eucalyptus/Stoneaucd6_output/ITS/tree_65.las"
file = "/home/jr/Documents/r-lidar/clients/fsinvestor/Rwanda/Eucalyptus/Stoneaucd6_output/ITS/tree_130.las"
file = "/home/jr/Documents/r-lidar/clients/fsinvestor/Rwanda/Eucalyptus/Stoneaucd6_output/ITS/tree_177.las"

#bug
file = "/home/jr/Documents/r-lidar inc/arbor/Tree bank/bug/405_erreur.laz"
file = "/home/jr/Documents/r-lidar inc/arbor/Tree bank/bug/23_bad_alloc.laz"

# Buttress
file = "/home/jr/Documents/r-lidar inc/arbor/Tree bank/bad-dbh/DUC0001-02_44.las"

# Bad
file = "/home/jr/Documents/r-lidar inc/arbor/Tree bank/bad-dbh/MDD01_006_1.laz" # !!!!
file = "/home/jr/Documents/r-lidar inc/arbor/Tree bank/bad-dbh/MDD03_011_1.laz" # !!!!

# gap
file =  "/home/jr/Documents/r-lidar inc/arbor/Tree bank/gap/BD2_DO10.laz"

# mini
file =  "/home/jr/Documents/r-lidar inc/arbor/Tree bank/mini/1464.las"

# Conifer
file = "/home/jr/Documents/r-lidar inc/arbor/Validation/Havelange/BD1_DO11.las"
file = "/home/jr/Documents/r-lidar inc/arbor/Validation/Havelange/BD1_DO12.las"

# Broken
file = "/home/jr/Documents/r-lidar inc/arbor/Tree bank/broken/tree_119.las"
file = "/home/jr/Documents/r-lidar inc/arbor/Tree bank/broken/tree_217.las"
file = "/home/jr/Documents/r-lidar inc/arbor/Tree bank/broken/tree_269.las"
file = "/home/jr/Documents/r-lidar inc/arbor/Tree bank/broken/P05_TLS_tree_7.asc"

# Broken branch
file = "/home/jr/Documents/r-lidar inc/arbor/Tree bank/overestimate/P09_TLS_3.las"

# Zambia
file = "/home/jr/Documents/r-lidar/clients/fsinvestor/Zambia/JasonFarm/its//tree_236.las" # !!
file = "/home/jr/Documents/r-lidar inc/arbor/Validation/Special1/las/Zambia_164.laz"
file = "/home/jr/Documents/r-lidar inc/arbor/Validation/Special1/las/Zambia_126.laz"
file = "/home/jr/Documents/r-lidar inc/arbor/Validation/Special1/las/Zambia_030.laz"
file = "/home/jr/Documents/r-lidar inc/arbor/Validation/Special1/las/Zambia_060.laz"
file = "/home/jr/Documents/r-lidar inc/arbor/Validation/Special1/las/Zambia_025.laz"
file = "/home/jr/Documents/r-lidar inc/arbor/Validation/Special1/las/Zambia_028.laz"
file = "/home/jr/Documents/r-lidar inc/arbor/Validation/Special1/las/Zambia_086.laz"
file = "/home/jr/Documents/r-lidar inc/arbor/Validation/Special1/las/Zambia_153.laz"
file = "/home/jr/Documents/r-lidar inc/arbor/Validation/Special1/las/Zambia_065.laz"
file = "/home/jr/Documents/r-lidar inc/arbor/Validation/Special1/las/Zambia_117.laz"
file = "/home/jr/Documents/r-lidar inc/arbor/Validation/Special1/las/Zambia_084.laz"
file = "/home/jr/Documents/r-lidar inc/arbor/Validation/Special1/las/Zambia_096.laz"

# Cameroun
file = "/home/jr/Documents/r-lidar inc/arbor/Validation/Momo/tree_laz/3.laz"
file = "/home/jr/Documents/r-lidar inc/arbor/Validation/Momo/tree_laz/12.laz"
file = "/home/jr/Documents/r-lidar inc/arbor/Validation/Momo/tree_laz/64.laz"
file = "/home/jr/Documents/r-lidar inc/arbor/Validation/Momo/tree_laz/100.laz"

# Validation
file = "/home/jr/Documents/r-lidar inc/arbor/Validation/Special1/las/PRF193_663.laz"
file = "/home/jr/Documents/r-lidar inc/arbor/Validation/Special1/las/PRF193_331.laz"
file = "/home/jr/Documents/r-lidar inc/arbor/Validation/Special1/las/PRF193_418.laz"
file = "/home/jr/Documents/r-lidar inc/arbor/Validation/Special1/las/PRF193_47.laz"
file = "/home/jr/Documents/r-lidar inc/arbor/Validation/Special1/las/PRF193_425.laz"

file = "/home/jr/Documents/r-lidar inc/arbor/Validation/Demol/pointclouds_clean/LXDC1.laz" # branches surestimée bon test.

file = "/home/jr/Documents/r-lidar inc/arbor/Validation/Special1/las/StA_09_tree_5574.laz"
file = "/home/jr/Documents/r-lidar inc/arbor/Validation/Special1/las/StA_09_tree_5664.laz" # bad bad bad

file = "/home/jr/Documents/r-lidar inc/arbor/Validation/Demol/pointclouds_clean/FEXC1.laz"

file = "/home/jr/Documents/r-lidar inc/arbor/Validation/TLSvMLS//MLS15m/P03_MLS15m_28.las"

# validation
file = "/home/jr/Documents/r-lidar inc/arbor/Validation/Gonzales/4_LidarTreePoinCloudData/MDD01_006_1.laz"

# Big tree non circular
file = "/home/jr/Documents/r-lidar inc/arbor/Tree bank/non-circular/tree_1114.las"
file = "/home/jr/Documents/r-lidar inc/arbor/Tree bank/non-circular/tree_54.las"  # prolongation!!
file = "/home/jr/Documents/r-lidar inc/arbor/Tree bank/non-circular/tree_681.las" # prolongation!!
file = "/home/jr/Documents/r-lidar inc/arbor/Tree bank/non-circular/tree_705.las" # prolongation!!
file = "/home/jr/Documents/r-lidar inc/arbor/Tree bank/non-circular/tree_848.las"
file = "/home/jr/Documents/r-lidar inc/arbor/Tree bank/non-circular/tree_939.las"

file = "/home/jr/Documents/r-lidar inc/arbor/Tree bank/non-circular/wytham_15240.las"
file = "/home/jr/Documents/r-lidar inc/arbor/Tree bank/non-circular/wytham_15408.las"


file = "/home/jr/Documents/r-lidar inc/arbor/Tree bank/urban/urban-karl.las" ##!!!


# ==== Single QSM ======

# QSM for a random tree
i = sample(seq_along(files), 1)
file <- files[i]
files = files[-i]
file
#tree <- lidR::readLAS(file)
#plot_semantic(tree)

{
  if (tools::file_ext(file) == "asc")
  {
    tree = data.table::fread(file)
    tree = LAS(tree) |> suppressWarnings()
    tree@data$foliage = 0L
    tree@data$hag = tree$Z - min(tree$Z)
    tree = hybrid_homogeneization(tree)
  }
  else if (tools::file_ext(file) == "xyz")
  {
    tree = data.table::fread(file)
    tree = LAS(tree) |> suppressWarnings()
    tree@data$foliage = 0L
    tree@data$hag = tree$Z - min(tree$Z)
    tree = hybrid_homogeneization(tree)
  }
  else
  {
    tree <- lidR::readLAS(file)

    if ("WL" %in% names(tree))
    {
      data.table::setnames(tree@data, "WL", "foliage")
      tree@data[, foliage := foliage - 1]
    }

    if ("TreeID" %in% names(tree))
    {
      data.table::setnames(tree@data, "TreeID", "treeID")
    }


    if (all(is.na(tree@data$foliage)))
        tree@data$foliage = 0L

    if (is.null(tree$hag))
      tree@data$hag = tree$Z - min(tree$Z)
  }


  tree@data$treeID  = as.integer(tree@data$treeID)
  tree@data$foliage = as.integer(tree@data$foliage)
  tree@data$passage = as.integer(tree@data$passage)

  t0 = Sys.time()
  params= arbor_parameters_default
  #params$qsm$skeleton_node_distance = 0.1
  #params$qsm$min_measurable_radius = 0.025
  qsm = qsm(tree, params)
  tf = Sys.time()
  print(difftime(tf, t0))

  qsm$radius[is.na(qsm$radius)] = 0

  V = qsm_volume(qsm)

  x = plot_semantic(tree)
  plot_qsm(qsm, add = x, cylinder = T, skel = F)

  x = plot_semantic(tree)
  plot_qsm(qsm, add = x, cylinder = T, skel = F, color = "quality")
}

as.data.frame(qsm)

merch = qsm_merchantable(qsm)
merch = merch[merch$dist_to_root < 2.5,]
plot(merch)

qsm_write(qsm, "~/Téléchargements/test.csv")
qsm_write(qsm, "~/Téléchargements/test.qsm")
v = qsm_read( "/home/jr/Documents/r-lidar inc/arbor/ArborStudio/ArborStudio/ressources/data/qsm1.qsm")
plot(v)

as.data.frame(merch)
as.data.frame(v)

qsm = qsm(tree, params)
plot(qsm)
as.data.frame(qsm)

qsf_write()

x = plot_semantic(tree)
plot(qsm, add = x, color = "quality")

x = plot_semantic(tree)
plot_qsm(qsm, add = x, cylinder = T)
plot_qsm(qsm, cylinder = T)
p

qsm_stem(qsm) |> qsm_nostump()|> qsm_merchantable(0.19/2, 1) |> qsm_volume()

d = 10
qsm$startX = qsm$startX - d
qsm$endX = qsm$endX - d
qsm$startY = qsm$startY - d
qsm$endY = qsm$endY - d
x <- plot_semantic(tree, size = 1)
plot(qsm, add = x, skel = F)


x = lidR::plot(tree)
rgl::close3d()

rgl::open3d()
rgl::bg3d("black")
rgl::mfrow3d(1, 2, sharedMouse = TRUE)

# --- Left subplot ---
# Capture current subscene id before plotting
left <- rgl::currentSubscene3d()
plot_qsm(qsm2, add = x, color = "branch_order")
rgl::useSubscene3d(left)

# --- Right subplot ---
rgl::next3d()
right <- rgl::currentSubscene3d()
x <- lidR::plot(tree, pal = "chocolate4", add = x)
plot_qsm(qsm2, add = x, color = "branch_order")

stem = qsm_stem(qsm)
merch = qsm_merchantable(qsm, 0.045)
merch_stem  = qsm_stem(qsm) |> qsm_merchantable(0.045)
plot_qsm(qsm)
plot_qsm(stem)
plot_qsm(merch)
plot_qsm(merch_stem)

find_selected_treeID <- function(las) {

  # Check that treeID exists
  if (!"treeID" %in% names(las@data)) {
    stop("LAS does not contain a 'treeID' attribute.")
  }

  # Open 3D plot
  x = arbor::plot_instance(las)

  message("Click a point in the rgl window...")

  # User clicks a point in 3D
  click <- rgl::select3d(button = "left")

  # Extract xyz coordinates
  xyz <- as.matrix(las@data[, c("X", "Y", "Z")])
  xyz[,1] =   xyz[,1] - x[1]
  xyz[,2] =   xyz[,2] - x[2]

  # Find selected points
  selected <- click(xyz[,1], xyz[,2], xyz[,3])

  if (!any(selected)) {
    message("No point selected.")
    return(NULL)
  }

  # If multiple points selected, use first one
  idx <- which(selected)[1]

  tree_id <- las@data$treeID[idx]

  message(sprintf("Selected treeID: %s", tree_id))

  return(tree_id)
}


### TEST

wood = arbor:::filter_tree(tree)

# ground
mhag = min(wood$hag)
gnd = wood[wood$hag < mhag + 0.1]

p = arbor_parameters_default

wood@data$dist2root = arbor:::dist2root(wood@data, gnd@data, p)
wood@data$dgroup = as.integer(cut(wood$dist2root, seq(floor(min(wood$dist2root)), ceiling(max(wood$dist2root)), by = 0.1)))
wood@data$iter = wood@data$dgroup


eps_value <- 0.3      # Adjust based on your distance units
minPts_value <- 1     # Minimum points to form a dense region
clust = function(x,y,z)
{
  dbscan::dbscan(cbind(x,y,z), eps = eps_value, minPts = minPts_value)$cluster
}


wood@data[, dbcl := clust(X,Y,Z), by = dgroup]
wood@data[, cl := as.integer(as.factor(paste0("g", dgroup, "_c", dbcl)))]
wood@data$cluster = wood@data$cl

neg = lidR::filter_poi(wood, dist2root < 0)
x = lidR::plot(wood, color = "dist2root")
lidR::plot(neg, add = x, pal = "pink", size = 6)
lidR::plot(gnd, add = x, size = 6)
lidR::plot(wood, color = "dgroup", pal = pastel.colors(2500))
lidR::plot(wood, color = "cl", pal = lidR::pastel.colors(2500))

wood = lidR::filter_poi(wood, dist2root >= 0, !is.na(iter))

plot(wood, color ="iter")
plot(wood, color ="cluster", pal = pastel.colors(2500))

x = plot(wood)
plot(filter_poi(wood, cluster == 50), add = x, pal = "red", size = 4)

w = wood@data[, .(X = mean(X), Y = mean(Y), Z = mean(Z)), by = .(iter, cluster)]
plot(LAS(w))

u = cpp_build_skeleton(wood@data, default_arbor_parameters$qsm$max_d)
plot_qsm(u)

library(rgl)

x = plot(wood)

# Open ONE device manually first
open3d()
bg3d("black")
mfrow3d(1, 2, sharedMouse = TRUE)

# --- Left subplot ---
# Capture current subscene id before plotting
left <- currentSubscene3d()
x <- plot(wood, pal = "chocolate4", add = x)
plot_qsm(u, add = x, color = "cyl_ID", pal = "red")
useSubscene3d(left)

# --- Right subplot ---
next3d()
right <- currentSubscene3d()
x <- plot(wood, pal = "chocolate4", add = x)
plot_qsm(qsm, add = x, pal = "blue", cylinder = F)

#lidR::plot(wood, color = "dist2root")
#lidR::plot(wood, color = "dgroup", pal = lidR::pastel.colors(250))
#lidR::plot(wood, color = "cl", pal = lidR::pastel.colors(7300))

nodes = wood@data[, .(X = mean(X), Y = mean(Y), Z = mean(Z), dgroup = dgroup[1], dist2root = dist2root[1]), by = cl]
nodes = lidR::LAS(nodes)
nodes = lidR::filter_poi(nodes, dist2root >= 0)
lidR::plot(nodes, color = "cl", pal = lidR::pastel.colors(7300))
lidR::plot(nodes, color = "dist2root")
lidR::plot(nodes, color = "dgroup", pal = lidR::pastel.colors(7300))

params$path_finder$k_neighborhood_connectivity = 10
params$path_finder$max_gap = 0.5
params$path_finder$distance_power = 3
spt = spt_to_dataframe(nodes@data, params)

library(rgl)

x = plot(wood)

# Open ONE device manually first
open3d()
bg3d("black")
mfrow3d(1, 2, sharedMouse = TRUE)

# --- Left subplot ---
# Capture current subscene id before plotting
left <- currentSubscene3d()
x <- plot(wood, pal = "chocolate4", add = x)
plot_qsm(spt, add = x, color = "pid1", pal = "red")
useSubscene3d(left)

# --- Right subplot ---
next3d()
right <- currentSubscene3d()
x <- plot(wood, pal = "chocolate4", add = x)
plot_qsm(qsm, add = x, pal = "blue", cylinder = F)

mst = qsm_topology(mst)
mst = qsm_architecture(mst)
plot_qsm(mst, color = "subtree_length")
lidR:::.pan3d(2)

# 1. Compute KNN
k = 10
knn = lidR::knn(nodes, k = k)
n_nodes = nrow(nodes)

# 2. Create long-format edge list
from_nodes = rep(1:n_nodes, each = k)
to_nodes = as.vector(t(knn$nn.index))
edge_costs = as.vector(t(knn$nn.dist^2)) # Using squared distance as weight

edges_df = data.frame(from = from_nodes, to = to_nodes, cost = edge_costs)

# 3. Create the igraph object (Undirected for shortest path calculations)
g <- graph_from_data_frame(edges_df, directed = FALSE)

# 4. Find the lowest node (the node with the minimum Z coordinate)
root_node <- which.min(nodes$Z)

# 5. Compute shortest paths from the root to all other nodes
# This returns the sequence of vertices/edges for the shortest paths
paths <- shortest_paths(g, from = root_node, to = V(g), weights = E(g)$cost, output = "epath")

# 6. Extract all edge IDs that belong to the Shortest Path Tree
spt_edge_ids <- unique(unlist(paths$epath))

# 7. Create a new graph containing only these shortest path edges
spt_graph <- subgraph_from_edges(g, eids = spt_edge_ids, delete.vertices = FALSE)

# 8. Convert the SPT graph back to a data frame of edges
spt_df = as_data_frame(spt_graph, what = "edges")
spt_df$from = as.integer(spt_df$from)
spt_df$to = as.integer(spt_df$to)
qsm_message
# 9. Map coordinates back to the SPT edges
startX = nodes$X[spt_df$from]
startY = nodes$Y[spt_df$from]
startZ = nodes$Z[spt_df$from]
endX = nodes$X[spt_df$to]
endY = nodes$Y[spt_df$to]
endZ = nodes$Z[spt_df$to]
pid1 = spt_df$from
pid2 = spt_df$to

# 10. Create final data.table and plot
spt_final = data.table::data.table(startX, startY, startZ, endX, endY, endZ, pid1, pid2)

library(rgl)

x = plot(wood)

# Open ONE device manually first
open3d()
bg3d("black")
mfrow3d(1, 2, sharedMouse = TRUE)

# --- Left subplot ---
# Capture current subscene id before plotting
left <- currentSubscene3d()
x <- plot(wood, pal = "chocolate4", add = x)
plot_qsm(spt_final, add = x, color = "pid1", pal = "red")
useSubscene3d(left)

# --- Right subplot ---
next3d()
right <- currentSubscene3d()
x <- plot(wood, pal = "chocolate4", add = x)
plot_qsm(qsm, add = x, pal = "blue", cylinder = F)

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
