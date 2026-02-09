#include "Adaptor.h"
#include "Grid3D.h"

// =======================
// PRE-PROCESSING
// =======================

std::vector<bool> homogeneization(const PointCloud& pc, double res, bool hybrid = true);

//[[Rcpp::export(rng = false)]]
Rcpp::LogicalVector C_homogeneization(Rcpp::DataFrame df, double res, bool hybrid = true)
{
  PointCloud pc(df);
  auto ans = homogeneization(pc, res);
  return(Rcpp::wrap(ans));
}

std::vector<float> anisotropy(PointCloud& adaptor, int k, int ncpu = 1);

//[[Rcpp::export(rng = false)]]
Rcpp::NumericVector C_anisotropy(Rcpp::DataFrame df,  int k, int ncpu = 1)
{
  PointCloud pc(df);
  auto ans = anisotropy(pc, k, ncpu);
  return(Rcpp::wrap(ans));
}

//[[Rcpp::export(rng = false)]]
Rcpp::IntegerVector C_connected_component(Rcpp::DataFrame df, double res, int connectivity)
{
  Rcpp::NumericVector X = df["X"];
  Rcpp::NumericVector Y = df["Y"];
  Rcpp::NumericVector Z = df["Z"];
  const double* x = X.begin();
  const double* y = Y.begin();
  const double* z = Z.begin();
  int n = X.size();
  Grid3D grid(x, y, z, n, res);
  return Rcpp::wrap(grid.connected_components(connectivity));
}

std::vector<bool> sor(const PointCloud& pc, unsigned int k, double m, int ncpu);

//[[Rcpp::export(rng = false)]]
Rcpp::LogicalVector C_sor(Rcpp::DataFrame df, unsigned int k, double m, int ncpu = 1)
{
  PointCloud pc(df);
  auto ans = sor(pc, k, m, ncpu);
  return(Rcpp::wrap(ans));
}

// ========================
// SEGMENTATION
// ========================

//[[Rcpp::export(rng = false)]]
SEXP build_semantic_graph(Rcpp::DataFrame dec, Rcpp::DataFrame target, Rcpp::DataFrame gnd, Rcpp::DataFrame master_seed, Rcpp::List params);

//[[Rcpp::export(rng = false)]]
SEXP build_instance_graph(Rcpp::DataFrame dec, Rcpp::DataFrame seed, Rcpp::DataFrame master_seed, Rcpp::List params);

//[[Rcpp::export(rng = false)]]
Rcpp::IntegerVector accumulate_passages(SEXP graph_ptr, int start_node, Rcpp::IntegerVector goal_nodes, int num_points);

//[[Rcpp::export(rng = false)]]
Rcpp::IntegerVector find_closest_node(SEXP graph_ptr, Rcpp::IntegerVector ids);

// ========================
// QSM
// ========================

//[[Rcpp::export(rng = false)]]
Rcpp::DataFrame qsm_topology_cpp(Rcpp::DataFrame qsm);

//[[Rcpp::export(rng = false)]]
Rcpp::DataFrame qsm_architecture_cpp(Rcpp::DataFrame qsm, int root_id = 1, bool use_volume = false);

//[[Rcpp::export(rng = false)]]
Rcpp::DataFrame cpp_build_skeleton(Rcpp::DataFrame data, double max_d);

//[[Rcpp::export(rng = false)]]
Rcpp::DataFrame qsm_layers_cpp(Rcpp::DataFrame df, double D);

//[[Rcpp::export(rng = false)]]
Rcpp::DataFrame qsm_simplify_cpp(Rcpp::DataFrame qsm, double max_length = 0.3);

//[[Rcpp::export(rng = false)]]
void qsm_write_cpp(Rcpp::DataFrame df, std::string filename, bool binary);

//[[Rcpp::export(rng = false)]]
Rcpp::DataFrame qsm_smooth_cpp(Rcpp::DataFrame df, int niter = 1, double th = 0);

//[[Rcpp::export(rng = false)]]
Rcpp::DataFrame qsm_prolongation_cpp(Rcpp::DataFrame df, double d, double L = 0.1);

//[[Rcpp::export(rng = false)]]
Rcpp::DataFrame qsm_measure_cpp(Rcpp::DataFrame pc, Rcpp::DataFrame df, float sarc = 180, float sins = 0.2, float sinl = 0.3, float srmeas = 0.05);

//[[Rcpp::export(rng = false)]]
Rcpp::DataFrame qsm_polynomial_fitting_cpp(Rcpp::DataFrame df, double tip_radius);

//[[Rcpp::export(rng = false)]]
Rcpp::DataFrame qsm_reconstruction_cpp(Rcpp::DataFrame df, double tip_radius);

//[[Rcpp::export(rng = false)]]
Rcpp::DataFrame read_adtree_skeleton(std::string filename);

// ========================
// QSF
// ========================

//[[Rcpp::export(rng = false)]]
void qsf_write_cpp(Rcpp::List x, std::string dir, std::string format, bool binary);

// ========================
// EXPERIMENTAL
// ========================

//[[Rcpp::export(rng = false)]]
Rcpp::DataFrame qsm_distances_cpp(Rcpp::DataFrame qsm_df, Rcpp::DataFrame pts_df);

//[[Rcpp::export(rng = false)]]
Rcpp::IntegerVector extract_tree_context_cpp(Rcpp::DataFrame las, int tree_id,  bool exclude_tree = false, int k = 10);

