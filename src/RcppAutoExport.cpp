// ! This file exist because Rcpp does not export things that are not in root folder.
// Everything is defined in R/RcppApi_*.cpp
#ifdef USING_R

#include <Rcpp.h>

// =======================
// PRE-PROCESSING
// =======================

//[[Rcpp::export(rng = false)]]
Rcpp::LogicalVector C_homogeneization(Rcpp::DataFrame df, double res, bool hybrid = true);

//[[Rcpp::export(rng = false)]]
Rcpp::NumericVector C_anisotropy(Rcpp::DataFrame df,  int k, int ncpu = 1);

//[[Rcpp::export(rng = false)]]
Rcpp::IntegerVector C_connected_component(Rcpp::DataFrame df, double res, int connectivity);

//[[Rcpp::export(rng = false)]]
Rcpp::LogicalVector C_sor(Rcpp::DataFrame df, unsigned int k, double m, int ncpu = 1);

// ========================
// SEGMENTATION
// ========================

//[[Rcpp::export(rng = false)]]
void segment_semantic_cpp(Rcpp::DataFrame core, Rcpp::DataFrame ground, Rcpp::List params);

//[[Rcpp::export(rng = false)]]
void segment_instance_cpp(Rcpp::DataFrame core, Rcpp::DataFrame seeds, Rcpp::List params);

//[[Rcpp::export(rng = false)]]
Rcpp::DataFrame find_seeds_cpp(Rcpp::DataFrame core, Rcpp::List params);

//[[Rcpp::export(rng = false)]]
Rcpp::IntegerVector accumulate_passages_cpp(Rcpp::DataFrame core, Rcpp::DataFrame gnd, Rcpp::List params);

//[[Rcpp::export(rng = false)]]
Rcpp::LogicalVector assign_wood_from_passage_cpp(Rcpp::DataFrame core, Rcpp::List params);

//[[Rcpp::export(rng = false)]]
Rcpp::LogicalVector assign_wood_from_high_likelihood_cpp(Rcpp::DataFrame core, Rcpp::List params);

//[[Rcpp::export(rng = false)]]
Rcpp::LogicalVector assign_wood_from_medium_likelihood_cpp(Rcpp::DataFrame core, Rcpp::List params);

//[[Rcpp::export(rng = false)]]
Rcpp::LogicalVector assign_wood_from_wood_dilatation_cpp(Rcpp::DataFrame core, Rcpp::List params);

//[[Rcpp::export(rng = false)]]
SEXP build_semantic_graph(Rcpp::DataFrame dec, Rcpp::DataFrame target, Rcpp::DataFrame gnd, Rcpp::List params);

//[[Rcpp::export(rng = false)]]
SEXP build_instance_graph(Rcpp::DataFrame dec, Rcpp::DataFrame seed, Rcpp::List params);

//[[Rcpp::export(rng = false)]]
Rcpp::IntegerVector accumulate_passages_old(SEXP graph_ptr, int start_node, Rcpp::IntegerVector goal_nodes, int num_points);

//[[Rcpp::export(rng = false)]]
Rcpp::IntegerVector find_closest_node(SEXP graph_ptr, Rcpp::IntegerVector ids);

// ========================
// SEEDS
// ========================

//[[Rcpp::export(rng = false)]]
Rcpp::DataFrame generate_cage_cpp(Rcpp::DataFrame circles, double decimation);


  //[[Rcpp::export(rng = false)]]
Rcpp::DataFrame detect_tree_circles_cpp(Rcpp::DataFrame wood_df, double resolution = 0.05, int connectivity = 26, int num_ransac_iterations = 400, double inlier_threshold = 0.02, int min_cluster_size = 20);


// ========================
// QSM
// ========================

//[[Rcpp::export(rng = false)]]
Rcpp::DataFrame qsm_cpp(Rcpp::DataFrame tree, Rcpp::List params);

//[[Rcpp::export(rng = false)]]
Rcpp::List qsf_cpp(Rcpp::DataFrame scene, Rcpp::List params);

//[[Rcpp::export(rng = false)]]
Rcpp::DataFrame qsm_layers_cpp(Rcpp::DataFrame df, double D);

//[[Rcpp::export(rng = false)]]
Rcpp::DataFrame qsm_cluster_cpp(Rcpp::DataFrame df, double cl_dist);

//[[Rcpp::export(rng = false)]]
Rcpp::DataFrame qsm_clean_tree_butt_cpp(Rcpp::DataFrame tree);

//[[Rcpp::export(rng = false)]]
Rcpp::DataFrame qsm_topology_cpp(Rcpp::DataFrame qsm);

//[[Rcpp::export(rng = false)]]
Rcpp::DataFrame qsm_architecture_cpp(Rcpp::DataFrame qsm, int root_id = 1, bool use_volume = false);

//[[Rcpp::export(rng = false)]]
Rcpp::DataFrame cpp_build_skeleton(Rcpp::DataFrame data, double max_d);

//[[Rcpp::export(rng = false)]]
Rcpp::DataFrame qsm_simplify_cpp(Rcpp::DataFrame qsm, double max_length = 0.3);

//[[Rcpp::export(rng = false)]]
void qsm_write_cpp(Rcpp::DataFrame df, std::string filename, bool binary);

//[[Rcpp::export(rng = false)]]
Rcpp::DataFrame qsm_smooth_cpp(Rcpp::DataFrame df, int niter = 1, double th = 0);

//[[Rcpp::export(rng = false)]]
Rcpp::DataFrame qsm_conic_allometry_cpp(Rcpp::DataFrame df, double R0, double tip_radius = 0.0025);

//[[Rcpp::export(rng = false)]]
double qsm_estimate_prolongation_cpp(Rcpp::DataFrame tree, Rcpp::DataFrame df);

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

//[[Rcpp::export(rng = false)]]
Rcpp::List qsm_tmesh_cpp(Rcpp::DataFrame df, int resolution);

//[[Rcpp::export(rng = false)]]
Rcpp::List qsm_qmesh_cpp(Rcpp::DataFrame df, int resolution);

// ========================
// QSF
// ========================

//[[Rcpp::export(rng = false)]]
void qsf_write_cpp(Rcpp::List x, std::string dir, std::string format, bool binary);

// ========================
// FITTING
// ========================

//[[Rcpp::export(rng = false)]]
Rcpp::List ransac_circle_cpp(Rcpp::NumericMatrix x, int num_iterations = 100, double inlier_threshold = 0.01, double early_exit = 1.0);

// ========================
// EXPERIMENTAL
// ========================

//[[Rcpp::export(rng = false)]]
Rcpp::DataFrame qsm_distances_cpp(Rcpp::DataFrame qsm_df, Rcpp::DataFrame pts_df);

//[[Rcpp::export(rng = false)]]
Rcpp::IntegerVector extract_tree_context_cpp(Rcpp::DataFrame las, int tree_id,  bool exclude_tree = false, int k = 10);

#endif
