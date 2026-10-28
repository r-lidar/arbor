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
    output[std::to_string(item.first)] = as_dataframe(item.second);
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

// Test helper: build a QSM and write it directly from the C++ graph, as a C++
// host (e.g. ArborStudio) does. qsm_write_cpp() rebuilds a clean graph from the
// data.frame and would not exercise the writer on a builder-produced graph. A
// leaf is removed so that the node ids are non-contiguous.
Rcpp::DataFrame qsm_build_and_write_cpp(Rcpp::DataFrame tree, Rcpp::List params, std::string filename)
{
  PointCloud pc(tree);
  arbor::settings::ArborParameters p = extract_arbor_params(params);
  QSM qsm = arbor::qsm::qsm(pc, p);

  for (const auto& [nid, node] : qsm.nodes())
  {
    if (qsm.children(nid).empty() && qsm.parent_node(nid) != -1)
    {
      qsm.remove_node(nid);
      break;
    }
  }

  qsm.write(filename, true);
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

Rcpp::DataFrame qsm_autorepair_cpp(Rcpp::DataFrame df)
{
  QSM qsm = as_qsm(df, true);
  return(as_dataframe(qsm));
}

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

Rcpp::List qsf_read_cpp(std::string filename)
{
  QSF qsf = QSF::read(filename);

  Rcpp::List output;
  for (const auto& item : qsf.get_qsm_map())
  {
    output[std::to_string(item.first)] = as_dataframe(item.second);
  }

  return output;
}

#endif
