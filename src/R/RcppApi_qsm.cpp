#ifdef USING_R

#include <Rcpp.h>

#include "myomp.h"
#include "arbor.h"
#include "QSMbuilder.h"
#include "RcppApi_wrappers.h"
#include "RcppApi_params.h"

using QSM = arbor::qsm::QSM;
using QSF = arbor::qsm::QSF;
using QSMbuilder = arbor::qsm::QSMbuilder;

Rcpp::DataFrame qsm_cpp(Rcpp::DataFrame tree, Rcpp::List params)
{
  PointCloud pc(tree);
  arbor::settings::ArborParameters p = extract_arbor_params(params);
  QSM qsm = arbor::qsm::qsm(pc, p);
  return as_dataframe(qsm);
}

Rcpp::List qsf_cpp(Rcpp::DataFrame scene, double min_height, Rcpp::List params)
{

  PointCloud pc(scene);
  arbor::settings::ArborParameters p = extract_arbor_params(params);
  QSF qsf = arbor::qsm::qsf(pc, min_height, p);

  Rcpp::List output;
  for (const auto& item : qsf.get_qsm_map())
  {
    output[item.first] = as_dataframe(item.second);
  }

  return output;
}


void qsm_write_cpp(Rcpp::DataFrame df, std::string filename, bool binary)
{
  QSM qsm = as_qsm(df);
  qsm.write(filename, binary);
}


Rcpp::DataFrame qsm_dbh_cpp(Rcpp::DataFrame df, double d = 1.3)
{
  double xyz[3];
  double n[3];
  double dbh;
  QSM qsm = as_qsm(df);
  dbh = qsm.dbh(d, xyz, n);

  return Rcpp::DataFrame::create(
    Rcpp::Named("dbh") = dbh,
    Rcpp::Named("x")   = xyz[0],
    Rcpp::Named("y")   = xyz[1],
    Rcpp::Named("z")   = xyz[2],
    Rcpp::Named("nx")  = n[0],
    Rcpp::Named("ny")  = n[1],
    Rcpp::Named("nz")  = n[2]
  );
}

Rcpp::DataFrame qsm_stem_cpp(Rcpp::DataFrame df)
{
  QSM qsm = as_qsm(df);
  QSM stem = qsm.stem();
  return(as_dataframe(stem));
}

Rcpp::DataFrame qsm_merchantable_cpp(Rcpp::DataFrame df, double merchandable_radius)
{
  QSM qsm = as_qsm(df);
  QSM stem = qsm.merchantable(merchandable_radius);
  return(as_dataframe(stem));
}

Rcpp::DataFrame qsm_layers_cpp(Rcpp::DataFrame df, double D)
{
  PointCloud pc(df);
  std::vector<std::pair<int, double>> res_pairs = QSMbuilder::layers(pc, D);

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

Rcpp::DataFrame qsm_cluster_cpp(Rcpp::DataFrame df, double cl_dist)
{
  // Extract the iter and dist columns from the dataframe
  Rcpp::IntegerVector iter = df["iter"];
  Rcpp::NumericVector dist = df["dist"];

  // Convert to std::vector<std::pair<int, double>>
  size_t n = iter.size();
  std::vector<std::pair<int, double>> iter_dist(n);
  for(size_t i = 0; i < n; ++i) {
    iter_dist[i] = {iter[i], dist[i]};
  }

  // Create PointCloud from dataframe
  PointCloud pc(df);

  // Call the clustering function
  std::vector<std::pair<int, double>> res_pairs = QSMbuilder::clusters(pc, iter_dist, cl_dist);

  // Unzip the vector of pairs into two separate vectors for Rcpp compatibility
  Rcpp::IntegerVector cluster_out(n);
  Rcpp::NumericVector radius_out(n);

  for(size_t i = 0; i < n; ++i) {
    cluster_out[i] = res_pairs[i].first;
    radius_out[i] = res_pairs[i].second;
  }

  return Rcpp::DataFrame::create(
    Rcpp::_["cluster"] = cluster_out,
    Rcpp::_["radius"] = radius_out
  );
}

Rcpp::DataFrame cpp_build_skeleton(Rcpp::DataFrame data, double max_d)
{
  Rcpp::IntegerVector iter = data["iter"];
  Rcpp::IntegerVector cluster = data["cluster"];

  std::vector<std::pair<int, int>> iter_cluster;
  iter_cluster.reserve(iter.size());
  for (size_t i = 0; i < iter.size(); i++)
  {
    iter_cluster.push_back({iter[i], cluster[i]});
  }

  PointCloud pc(data);

  QSM qsm;
  QSMbuilder b(qsm);
  b.build_skeleton(pc, iter_cluster, max_d);
  return as_dataframe(qsm);
}


 /*Rcpp::DataFrame qsm_clean_tree_butt_cpp(Rcpp::DataFrame tree)
 {
   PointCloud pc(tree);
   PointCloud res = QSMbuilder::clean_tree_butt(pc);
   return(as_dataframe(res));
 }*/

 Rcpp::DataFrame qsm_topology_cpp(Rcpp::DataFrame df)
 {
   QSM qsm = as_qsm(df);
   QSMbuilder b(qsm);
   b.compute_topology();

   // Extract parent_ID for each cylinder (ordered by cyl_ID)
   int n = (int)qsm.edges().size();
   Rcpp::IntegerVector parent_ID(n);
   for (const auto& [eid, einfo] : qsm.edges())
   {
     int cid = einfo.data.cyl_ID;
     if (cid >= 1 && cid <= n)
       parent_ID[cid - 1] = einfo.data.parent_ID;
   }

   df["parent_ID"] = parent_ID;
   return df;
 }


Rcpp::DataFrame qsm_architecture_cpp(Rcpp::DataFrame df, bool use_volume = false)
{
  QSM qsm = as_qsm(df);
  QSMbuilder b(qsm);
  b.compute_architecture(use_volume);

  int n = (int)qsm.edges().size();
  Rcpp::NumericVector subtree_length(n);
  Rcpp::IntegerVector axis_ID(n);
  Rcpp::IntegerVector branching_order(n);

  for (const auto& [eid, einfo] : qsm.edges())
  {
    int cid = einfo.data.cyl_ID;
    if (cid >= 1 && cid <= n)
    {
      int idx = cid - 1;
      subtree_length[idx]   = einfo.data.subtree_length;
      axis_ID[idx]          = einfo.data.axis_ID;
      branching_order[idx]  = einfo.data.branch_order;
    }
  }

  df["axis_ID"]        = axis_ID;
  df["branch_order"]   = branching_order;
  df["subtree_length"] = subtree_length;

  return df;
}

/*Rcpp::DataFrame qsm_smooth_cpp(Rcpp::DataFrame df, int niter = 1, double th = 0)
{
  QSM qsm = as_qsm(df);
  QSMbuilder b(qsm);
  b.smooth_skeleton(niter, th);

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

double qsm_estimate_prolongation_cpp(Rcpp::DataFrame tree, Rcpp::DataFrame df)
{
  QSM qsm = as_qsm(df);
  PointCloud pc(tree);
  QSMbuilder b(qsm);
  b.estimate_prolongation(pc);
  return b.prolongation_distance;
}

Rcpp::DataFrame qsm_prolongation_cpp(Rcpp::DataFrame df, double d, double L = 0.1)
{
  QSM qsm = as_qsm(df);
  QSMbuilder b(qsm);
  b.prolongate(d, L);

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
  QSMbuilder b(qsm);
  b.measure_radii(tree, sarc, sins, sinl, srmeas);
  return as_dataframe(qsm);
}*/

Rcpp::DataFrame qsm_polynomial_fitting_cpp(Rcpp::DataFrame df, double tip_radius)
{
  QSM qsm = as_qsm(df);
  QSMbuilder b(qsm);
  b.polynomial_fitting(tip_radius);
  return as_dataframe(qsm);
}

/*Rcpp::DataFrame qsm_conic_allometry_cpp(Rcpp::DataFrame df, double R0, double tip_radius)
{
  QSM qsm = as_qsm(df);
  QSMbuilder b(qsm);
  b.conic_allometry(R0, tip_radius);
  return as_dataframe(qsm);
}

Rcpp::DataFrame qsm_reconstruction_cpp(Rcpp::DataFrame df, double tip_radius)
{
  QSM qsm = as_qsm(df);
  QSMbuilder b(qsm);
  b.reconstruct_missing_radii(tip_radius);
  return as_dataframe(qsm);
}*/

Rcpp::List qsm_mesh_cpp(Rcpp::DataFrame df, int resolution)
{
  QSM qsm = as_qsm(df);

  std::vector<std::array<double, 3>> verts_vec;
  std::vector<std::array<int, 4>> faces_vec;
  std::vector<int> node_ids;
  qsm.qmesh(verts_vec, faces_vec, node_ids, resolution);

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

  Rcpp::IntegerVector nodes = Rcpp::wrap(node_ids);

  return Rcpp::List::create(Rcpp::Named("vertices") = vertices,  Rcpp::Named("indices")  = indices, Rcpp::Named("NodeID") = nodes);
}


void qsf_write_cpp(Rcpp::List x, std::string dir, std::string format, bool binary)
{
  QSF qsf = as_qsf(x);
  qsf.write(dir, format, binary);
}

Rcpp::DataFrame qsm_simplify_cpp(Rcpp::DataFrame qsm, double max_length = 0.3)
{
  Rcpp::IntegerVector cyl_ID = qsm["cyl_ID"];
  Rcpp::IntegerVector parent_ID = qsm["parent_ID"];
  Rcpp::NumericVector startX = qsm["startX"];
  Rcpp::NumericVector startY = qsm["startY"];
  Rcpp::NumericVector startZ = qsm["startZ"];
  Rcpp::NumericVector endX = qsm["endX"];
  Rcpp::NumericVector endY = qsm["endY"];
  Rcpp::NumericVector endZ = qsm["endZ"];
  Rcpp::NumericVector radius = qsm["radius"];
  Rcpp::NumericVector cyl_length = qsm["cyl_length"];

  int n = cyl_ID.size();

  // Build parent-to-children map and cyl_ID -> row index
  std::unordered_map<int, std::vector<int>> children;
  std::unordered_map<int, int> id_to_row;
  std::unordered_map<int, int> child_count;

  for (int i = 0; i < n; ++i)
  {
    id_to_row[cyl_ID[i]] = i;
    children[parent_ID[i]].push_back(cyl_ID[i]);
    if (parent_ID[i] != 0)
      child_count[parent_ID[i]]++;
  }

  // Identify branching points and tips
  std::unordered_set<int> branching, tips;
  for (const auto& kv : child_count)
  {
    if (kv.second > 1) branching.insert(kv.first);
  }

  for (int i = 0; i < n; ++i)
  {
    if (children.find(cyl_ID[i]) == children.end())
      tips.insert(cyl_ID[i]);
  }

  std::unordered_set<int> important;
  for (int i = 0; i < n; ++i)
  {
    if (branching.count(cyl_ID[i]) || branching.count(parent_ID[i]) || tips.count(cyl_ID[i]))
      important.insert(cyl_ID[i]);
  }

  std::vector<bool> visited(n, false);
  std::vector<int> new_cyl_ID, new_parent_ID, new_original_row;
  std::vector<double> new_startX, new_startY, new_startZ, new_endX, new_endY, new_endZ, new_radius, new_length;

  std::unordered_map<int, int> old_to_new_id;
  int next_id = 1;

  for (int i = 0; i < n; ++i)
  {
    if (visited[i]) continue;

    int row = i;
    std::vector<int> chain = {row};
    visited[row] = true;

    // Grow forward
    int current = cyl_ID[row];
    while (children[current].size() == 1)
    {
      int child = children[current][0];
      int child_row = id_to_row[child];
      if (visited[child_row] || important.count(child)) break;
      chain.push_back(child_row);
      visited[child_row] = true;
      current = child;
    }

    int start = 0;
    while (start < (int)chain.size())
    {
      int idx = chain[start];
      int last = start;

      int s_idx = chain[start];
      double sx = startX[s_idx], sy = startY[s_idx], sz = startZ[s_idx];
      double ex = endX[s_idx], ey = endY[s_idx], ez = endZ[s_idx];
      double r = radius[s_idx];

      for (int j = start + 1; j < (int)chain.size(); ++j)
      {
        int k = chain[j];
        double dx = endX[k] - sx;
        double dy = endY[k] - sy;
        double dz = endZ[k] - sz;
        double d = std::sqrt(dx * dx + dy * dy + dz * dz);
        if (d > max_length) break;

        ex = endX[k];
        ey = endY[k];
        ez = endZ[k];
        r = (r + radius[k]) / 2.0;
        last = j;
      }

      double len = std::sqrt((ex - sx) * (ex - sx) + (ey - sy) * (ey - sy) + (ez - sz) * (ez - sz));

      new_cyl_ID.push_back(next_id);
      old_to_new_id[cyl_ID[chain[last]]] = next_id;

      new_parent_ID.push_back(0);  // temp, fix later
      new_startX.push_back(sx);
      new_startY.push_back(sy);
      new_startZ.push_back(sz);
      new_endX.push_back(ex);
      new_endY.push_back(ey);
      new_endZ.push_back(ez);
      new_radius.push_back(r);
      new_length.push_back(len);

      // Preserve original row index (1-based for R)
      new_original_row.push_back(chain[last] + 1);

      ++next_id;
      start = last + 1;
    }
  }

  // Fix parent_IDs
  for (size_t i = 0; i < new_cyl_ID.size(); ++i)
  {
    int old_idx = -1;
    for (int j = 0; j < n; ++j)
    {
      if (old_to_new_id[cyl_ID[j]] == new_cyl_ID[i])
        old_idx = j;
    }
    if (old_idx >= 0 && parent_ID[old_idx] != 0)
    {
      int old_pid = parent_ID[old_idx];
      auto it = old_to_new_id.find(old_pid);
      if (it != old_to_new_id.end())
        new_parent_ID[i] = it->second;
    }
  }

  return Rcpp::DataFrame::create(
    Rcpp::Named("startX") = new_startX,
    Rcpp::Named("startY") = new_startY,
    Rcpp::Named("startZ") = new_startZ,
    Rcpp::Named("endX") = new_endX,
    Rcpp::Named("endY") = new_endY,
    Rcpp::Named("endZ") = new_endZ,
    Rcpp::Named("radius") = new_radius,
    Rcpp::Named("cyl_ID") = new_cyl_ID,
    Rcpp::Named("parent_ID") = new_parent_ID,
    Rcpp::Named("original_row") = new_original_row  // <-- added column
  );
}

#endif
