#include "Adaptor.h"
#include "Grid3D.h"
#include "ransac.h"

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
  PointCloud pc(df);
  Grid3D grid(pc, res);
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
// FITTING
// ========================

//[[Rcpp::export(rng = false)]]
Rcpp::List ransac_circle_cpp(Rcpp::NumericMatrix x, int num_iterations = 100, double inlier_threshold = 0.01, double early_exit = 1.0)
{
  MatrixAdaptor pc(x);
  RansacCircle rc(num_iterations, inlier_threshold, early_exit);
  for (int i = 0 ; i < pc.point_count() ; i++)
    rc.add_point(pc.get_x(i), pc.get_y(i), pc.get_z(i));
  rc.find_circle();

  std::array<double, 3> center = rc.get_center();
  double radius = rc.get_radius();
  double inlier_pct = rc.get_inlier_percentage();
  double inside_pct = rc.get_inside_percentage();
  double arc_deg = rc.get_arc_coverage();
  const std::vector<int>& inliers = rc.get_inliers();

  // Calculate CFQI (matching R implementation)
  double arc_score = arc_deg / 360.0;
  double score_inside = (inside_pct == 0.0) ? 1.0 : (1.0 - inside_pct);

  double w_arc = 0.4;
  double w_inside = 0.1;
  double w_inlier = 0.4;
  double cfqi = (w_arc * arc_score + w_inlier * inlier_pct + w_inside / score_inside);

  Rcpp::IntegerVector r_inliers = Rcpp::wrap(inliers);

  return Rcpp::List::create(
    Rcpp::Named("center_x") = center[0],
    Rcpp::Named("center_y") = center[1],
    Rcpp::Named("radius") = radius,
    Rcpp::Named("z") = center[2],
    Rcpp::Named("covered_arc_degree") = arc_deg,
    Rcpp::Named("percentage_inlier") = inlier_pct,
    Rcpp::Named("percentage_inside") = inside_pct,
    Rcpp::Named("inliers") = r_inliers+1,
    Rcpp::Named("CFQI") = cfqi
  );
}

// ========================
// EXPERIMENTAL
// ========================

//[[Rcpp::export(rng = false)]]
Rcpp::DataFrame qsm_distances_cpp(Rcpp::DataFrame qsm_df, Rcpp::DataFrame pts_df);

//[[Rcpp::export(rng = false)]]
Rcpp::IntegerVector extract_tree_context_cpp(Rcpp::DataFrame las, int tree_id,  bool exclude_tree = false, int k = 10);

