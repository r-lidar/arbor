#include <Rcpp.h>

#include "myomp.h"
#include "arbor.h"
#include "Rwrappers.h"


Rcpp::DataFrame qsm_layers_cpp(Rcpp::DataFrame df, double D)
{
  PointCloud pc(df);
  std::vector<std::pair<int, double>> res_pairs = QSM::layers(pc, D);

  // Unzip the vector of pairs into two separate vectors for Rcpp compatibility
  size_t n = res_pairs.size();
  Rcpp::IntegerVector iter_out(n);
  Rcpp::NumericVector dist_out(n);

  for(size_t i = 0; i < n; ++i) {
    iter_out[i] = res_pairs[i].first;
    dist_out[i] = res_pairs[i].second;
  }

  return Rcpp::DataFrame::create(
    Rcpp::_["iter"] = iter_out,
    Rcpp::_["dist"] = dist_out
  );
}

Rcpp::DataFrame qsm_topology_cpp(Rcpp::DataFrame df)
{
  QSM qsm = as_qsm(df);
  qsm.compute_topology();

  Rcpp::IntegerVector parent_ID(qsm.size());
  for (const auto &[i, c] : qsm)
    parent_ID[i-1] = c.parent_ID;

  df["parent_ID"] = parent_ID;

  return df;
}

Rcpp::DataFrame qsm_architecture_cpp(Rcpp::DataFrame df, int root_id = 1, bool use_volume = false)
{
  QSM qsm = as_qsm(df);
  qsm.compute_architecture(root_id, use_volume); // Rcpp catches exceptions

  // Prepare new columns
  int n = (int)qsm.size();
  Rcpp::NumericVector subtree_length(n);
  Rcpp::IntegerVector axis_ID(n);
  Rcpp::IntegerVector branching_order(n);
  for (const auto &[i, c] : qsm)
  {
    //length[i-1]            = c.length();
    subtree_length[i-1]    = c.subtree_length;
    axis_ID[i-1]           = c.axis_ID;
    branching_order[i-1]   = c.branch_order;
  }

  df["axis_ID"]          = axis_ID;
  df["branch_order"]     = branching_order;
  df["subtree_length"]   = subtree_length;

  return df;
}

void qsm_write_cpp(Rcpp::DataFrame df, std::string filename, bool binary)
{
  QSM qsm = as_qsm(df);
  qsm.write(filename, binary);
}


Rcpp::DataFrame qsm_smooth_cpp(Rcpp::DataFrame df, int niter = 1, double th = 0)
{
  QSM qsm = as_qsm(df);
  qsm.smooth_skeleton(niter, th);
  Rcpp::DataFrame ans = as_dataframe(qsm);

  // Keep only columns present in original df
  Rcpp::CharacterVector keep = df.names();
  Rcpp::CharacterVector ans_names = ans.names();

  std::vector<std::string> to_keep;
  to_keep.reserve(keep.size());

  for (auto &nm : keep)
  {
    if (std::find(ans_names.begin(), ans_names.end(), std::string(nm)) != ans_names.end())
    {
      to_keep.push_back(std::string(nm));
    }
  }

  // Subset ans
  Rcpp::List out(to_keep.size());
  for (size_t i = 0; i < to_keep.size(); i++)
  {
    out[i] = ans[to_keep[i]];
  }
  out.attr("names") = to_keep;
  out.attr("class") = "data.frame";
  out.attr("row.names") = ans.attr("row.names");

  return out;
}

Rcpp::DataFrame qsm_prolongation_cpp(Rcpp::DataFrame df, double d, double L = 0.1)
{
  QSM qsm = as_qsm(df);
  qsm.prolongate(d, L);
  Rcpp::DataFrame ans = as_dataframe(qsm);

  // Keep only columns present in original df
  Rcpp::CharacterVector keep = df.names();
  Rcpp::CharacterVector ans_names = ans.names();

  std::vector<std::string> to_keep;
  to_keep.reserve(keep.size());

  for (auto &nm : keep)
  {
    if (std::find(ans_names.begin(), ans_names.end(), std::string(nm)) != ans_names.end())
    {
      to_keep.push_back(std::string(nm));
    }
  }

  // Subset ans
  Rcpp::List out(to_keep.size());
  for (size_t i = 0; i < to_keep.size(); i++)
  {
    out[i] = ans[to_keep[i]];
  }
  out.attr("names") = to_keep;
  out.attr("class") = "data.frame";
  out.attr("row.names") = ans.attr("row.names");

  return out;
}

Rcpp::DataFrame qsm_measure_cpp(Rcpp::DataFrame pc, Rcpp::DataFrame df, float sarc = 180, float sins = 0.2, float sinl = 0.3, float srmeas = 0.05)
{
  PointCloud tree(pc);
  QSM qsm = as_qsm(df);
  qsm.measure_radii(tree, sarc, sins, sinl, srmeas);
  Rcpp::DataFrame ans = as_dataframe(qsm);
  return ans;
}

Rcpp::DataFrame qsm_polynomial_fitting_cpp(Rcpp::DataFrame df, double tip_radius)
{
  QSM qsm = as_qsm(df);
  qsm.polynomial_fitting(tip_radius);
  Rcpp::DataFrame ans = as_dataframe(qsm);
  return ans;
}

Rcpp::DataFrame qsm_reconstruction_cpp(Rcpp::DataFrame df, double tip_radius)
{
  QSM qsm = as_qsm(df);
  qsm.reconstruct_missing_radii(tip_radius);
  Rcpp::DataFrame ans = as_dataframe(qsm);
  return ans;
}

Rcpp::List qsm_tmesh_cpp(Rcpp::DataFrame df, int resolution)
{
  QSM qsm = as_qsm(df);

  std::vector<std::array<double, 3>> verts_vec;
  std::vector<std::array<int, 3>> faces_vec;
  qsm.tmesh(verts_vec, faces_vec, resolution);

  Rcpp::NumericMatrix vertices(3, verts_vec.size());
  for (size_t i = 0; i < verts_vec.size(); ++i) {
    vertices(0, i) = verts_vec[i][0];
    vertices(1, i) = verts_vec[i][1];
    vertices(2, i) = verts_vec[i][2];
  }

  Rcpp::IntegerMatrix indices(3, faces_vec.size());
  for (size_t i = 0; i < faces_vec.size(); ++i) {
    indices(0, i) = faces_vec[i][0] + 1;
    indices(1, i) = faces_vec[i][1] + 1;
    indices(2, i) = faces_vec[i][2] + 1;
  }

  return Rcpp::List::create(Rcpp::Named("vertices") = vertices,  Rcpp::Named("indices")  = indices);
}

Rcpp::List qsm_qmesh_cpp(Rcpp::DataFrame df, int resolution)
{
  QSM qsm = as_qsm(df);

  std::vector<std::array<double, 3>> verts_vec;
  std::vector<std::array<int, 4>> faces_vec;
  qsm.qmesh(verts_vec, faces_vec, resolution);

  Rcpp::NumericMatrix vertices(3, verts_vec.size());
  for (size_t i = 0; i < verts_vec.size(); ++i) {
    vertices(0, i) = verts_vec[i][0];
    vertices(1, i) = verts_vec[i][1];
    vertices(2, i) = verts_vec[i][2];
  }

  Rcpp::IntegerMatrix indices(4, faces_vec.size());
  for (size_t i = 0; i < faces_vec.size(); ++i) {
    indices(0, i) = faces_vec[i][0] + 1;
    indices(1, i) = faces_vec[i][1] + 1;
    indices(2, i) = faces_vec[i][2] + 1;
    indices(3, i) = faces_vec[i][3] + 1;
  }

  return Rcpp::List::create(Rcpp::Named("vertices") = vertices,  Rcpp::Named("indices")  = indices);
}


void qsf_write_cpp(Rcpp::List x, std::string dir, std::string format, bool binary)
{
  QSF qsf = as_qsf(x);
  qsf.write(dir, format, binary);
}
