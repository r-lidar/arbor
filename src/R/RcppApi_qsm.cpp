/**
 * @file RcppApi_qsm.cpp
 * Project: Arbor
 *
 * Copyright (C) 2026 Jean-Romain Roussel (r-lidar) <info @ r-lidar.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifdef USING_R

#include <Rcpp.h>

#include <array>

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

Rcpp::DataFrame qsm_read_cpp(std::string filename)
{
  QSM qsm;
  qsm.read(filename);
  return as_dataframe(qsm);
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

Rcpp::DataFrame qsm_merchantable_cpp(Rcpp::DataFrame df, double merchantable_radius, double merchantable_length)
{
  QSM qsm = as_qsm(df);
  QSM stem = qsm.merchantable(merchantable_radius, merchantable_length);
  return(as_dataframe(stem));
}

Rcpp::DataFrame qsm_cluster_cpp(Rcpp::DataFrame df, double cl_dist)
{
  // Extract the iter and dist columns from the dataframe
  Rcpp::IntegerVector iterdf = df["iter"];

  // Convert to std::vector<std::pair<int, double>>
  size_t n = iterdf.size();
  std::vector<int> iter = Rcpp::as<std::vector<int>>(iterdf);


  // Create PointCloud from dataframe
  PointCloud pc(df);

  // Call the clustering function
  std::vector<int> cl = QSMbuilder::cluster(pc, iter, cl_dist);

  // Unzip the vector of pairs into two separate vectors for Rcpp compatibility
  Rcpp::IntegerVector cluster_out = Rcpp::wrap(cl);

  return Rcpp::DataFrame::create(
    Rcpp::_["cluster"] = cluster_out
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
 /*  QSM qsm = as_qsm(df);
   QSMbuilder b(qsm);
   b.compute_topology();

   // Extract source for each cylinder (ordered by id)
   int n = (int)qsm.edges().size();
   Rcpp::IntegerVector source(n);
   for (const auto& [eid, einfo] : qsm.edges())
   {
     int cid = einfo.data.id;
     if (cid >= 1 && cid <= n)
       source[cid - 1] = einfo.data.source;
   }

   df["source"] = source;*/
   return df;
 }


Rcpp::DataFrame qsm_architecture_cpp(Rcpp::DataFrame df, bool use_volume = false)
{
/*  QSM qsm = as_qsm(df);
  QSMbuilder b(qsm);
  b.compute_architecture(use_volume);

  int n = (int)qsm.edges().size();
  Rcpp::NumericVector subtree_length(n);
  Rcpp::IntegerVector axis_id(n);
  Rcpp::IntegerVector branching_order(n);

  for (const auto& [eid, einfo] : qsm.edges())
  {
    int cid = einfo.data.id;
    if (cid >= 1 && cid <= n)
    {
      int idx = cid - 1;
      subtree_length[idx]   = einfo.data.subtree_length;
      axis_id[idx]          = einfo.data.axis_id;
      branching_order[idx]  = einfo.data.branch_order;
    }
  }

  df["axis_id"]        = axis_id;
  df["branch_order"]   = branching_order;
  df["subtree_length"] = subtree_length;*/

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
  Rcpp::IntegerVector id = qsm["id"];
  Rcpp::IntegerVector source = qsm["source"];
  Rcpp::NumericVector startX = qsm["startX"];
  Rcpp::NumericVector startY = qsm["startY"];
  Rcpp::NumericVector startZ = qsm["startZ"];
  Rcpp::NumericVector endX = qsm["endX"];
  Rcpp::NumericVector endY = qsm["endY"];
  Rcpp::NumericVector endZ = qsm["endZ"];
  Rcpp::NumericVector radius = qsm["radius"];
  Rcpp::NumericVector cyl_length = qsm["cyl_length"];

  int n = id.size();

  // Build parent-to-children map and id -> row index
  std::unordered_map<int, std::vector<int>> children;
  std::unordered_map<int, int> id_to_row;
  std::unordered_map<int, int> child_count;

  for (int i = 0; i < n; ++i)
  {
    id_to_row[id[i]] = i;
    children[source[i]].push_back(id[i]);
    if (source[i] != 0)
      child_count[source[i]]++;
  }

  // Identify branching points and tips
  std::unordered_set<int> branching, tips;
  for (const auto& kv : child_count)
  {
    if (kv.second > 1) branching.insert(kv.first);
  }

  for (int i = 0; i < n; ++i)
  {
    if (children.find(id[i]) == children.end())
      tips.insert(id[i]);
  }

  std::unordered_set<int> important;
  for (int i = 0; i < n; ++i)
  {
    if (branching.count(id[i]) || branching.count(source[i]) || tips.count(id[i]))
      important.insert(id[i]);
  }

  std::vector<bool> visited(n, false);
  std::vector<int> new_id, new_source, new_original_row;
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
    int current = id[row];
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

      new_id.push_back(next_id);
      old_to_new_id[id[chain[last]]] = next_id;

      new_source.push_back(0);  // temp, fix later
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

  // Fix sources
  for (size_t i = 0; i < new_id.size(); ++i)
  {
    int old_idx = -1;
    for (int j = 0; j < n; ++j)
    {
      if (old_to_new_id[id[j]] == new_id[i])
        old_idx = j;
    }
    if (old_idx >= 0 && source[old_idx] != 0)
    {
      int old_pid = source[old_idx];
      auto it = old_to_new_id.find(old_pid);
      if (it != old_to_new_id.end())
        new_source[i] = it->second;
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
    Rcpp::Named("id") = new_id,
    Rcpp::Named("source") = new_source,
    Rcpp::Named("original_row") = new_original_row  // <-- added column
  );
}

Rcpp::IntegerVector qsm_get_cylID(Rcpp::DataFrame qsmdf, Rcpp::DataFrame pcdf)
{
  QSM qsm = as_qsm(qsmdf);
  PointCloud pc(pcdf);
  auto map = QSMbuilder::group_points_by_edge(pc);

  Rcpp::IntegerVector res(pc.size());
  for (const auto& e : map)
  {
    for (auto i : e.second)
    {
      res[i] = e.first;
    }
  }

  return res;
}

#endif
