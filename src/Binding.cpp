#include <Rcpp.h>

// =======================
// PRE-PROCESSING
// =======================

// [[Rcpp::export]]
Rcpp::LogicalVector C_voxel_barycenter_decimate(Rcpp::NumericVector X, Rcpp::NumericVector Y, Rcpp::NumericVector Z, Rcpp::NumericVector id);

// ========================
// SEGMENTATION
// ========================

// [[Rcpp::export]]
SEXP init_graph();

// [[Rcpp::export]]
SEXP build_graph(Rcpp::DataFrame graph_df);

// [[Rcpp::export]]
Rcpp::IntegerVector accumulate_passages(SEXP graph_ptr, Rcpp::IntegerVector start_nodes, Rcpp::IntegerVector goal_nodes, int num_points);

// [[Rcpp::export]]
Rcpp::List find_closest_ground(SEXP graph_ptr, Rcpp::IntegerVector ground_node_ids);



// [[Rcpp::export]]
Rcpp::DataFrame compute_point_network_cpp(
    Rcpp::DataFrame dec,
    int k,
    double max_gap = 1.0,
    Rcpp::Nullable<Rcpp::LogicalVector> wood_mask = R_NilValue,
    Rcpp::Nullable<Rcpp::List> cost_factors = R_NilValue,
    double power = 3.0,
    bool downward = false);

// ========================
// QSM
// ========================

// [[Rcpp::export]]
Rcpp::DataFrame cpp_compute_architecture(Rcpp::DataFrame qsm, int root_id = 1);

// [[Rcpp::export]]
Rcpp::DataFrame cpp_build_skeleton(Rcpp::DataFrame data, double max_d);

// [[Rcpp::export]]
Rcpp::DataFrame cpp_build_skeleton_old(Rcpp::DataFrame data, double max_d);

// [[Rcpp::export]]
Rcpp::DataFrame cpp_compute_layers(Rcpp::NumericMatrix coords, double D);

// [[Rcpp::export]]
Rcpp::DataFrame cpp_compute_topology(Rcpp::DataFrame skel);

// [[Rcpp::export]]
Rcpp::DataFrame qsm_simplify_cpp(Rcpp::DataFrame qsm, double max_length = 0.3);

// [[Rcpp::export]]
Rcpp::List cpp_smooth_skeleton(Rcpp::DataFrame qsm, int niter = 1, double th = 0);

