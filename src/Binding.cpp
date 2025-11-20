#include <Rcpp.h>

// =======================
// PRE-PROCESSING
// =======================

// [[Rcpp::export]]
Rcpp::LogicalVector C_voxel_barycenter_decimate(Rcpp::NumericVector X, Rcpp::NumericVector Y, Rcpp::NumericVector Z, double res);

// ========================
// SEGMENTATION
// ========================

// [[Rcpp::export]]
SEXP build_semantic_graph(Rcpp::DataFrame dec, Rcpp::DataFrame target, Rcpp::DataFrame gnd, Rcpp::DataFrame master_seed, Rcpp::List params);

// [[Rcpp::export]]
SEXP build_instance_graph(Rcpp::DataFrame dec, Rcpp::DataFrame seed, Rcpp::DataFrame master_seed, Rcpp::List params);


// [[Rcpp::export]]
Rcpp::IntegerVector accumulate_passages(SEXP graph_ptr, int start_node, Rcpp::IntegerVector goal_nodes, int num_points);

// [[Rcpp::export]]
Rcpp::IntegerVector find_closest_node(SEXP graph_ptr, Rcpp::IntegerVector ids);

// ========================
// QSM
// ========================

// [[Rcpp::export]]
Rcpp::DataFrame qsm_topology_cpp(Rcpp::DataFrame qsm);

// [[Rcpp::export]]
Rcpp::DataFrame qsm_architecture_cpp(Rcpp::DataFrame qsm, int root_id = 1);

// [[Rcpp::export]]
Rcpp::DataFrame cpp_build_skeleton(Rcpp::DataFrame data, double max_d);

// [[Rcpp::export]]
Rcpp::DataFrame cpp_compute_layers(Rcpp::NumericMatrix coords, double D);

// [[Rcpp::export]]
Rcpp::DataFrame qsm_simplify_cpp(Rcpp::DataFrame qsm, double max_length = 0.3);

// [[Rcpp::export]]
void qsm_write_cpp(Rcpp::DataFrame df, std::string filename, bool binary);

// [[Rcpp::export]]
Rcpp::DataFrame qsm_smooth_cpp(Rcpp::DataFrame df, int niter = 1, double th = 0);

// [[Rcpp::export]]
Rcpp::DataFrame read_adtree_skeleton(std::string filename);
