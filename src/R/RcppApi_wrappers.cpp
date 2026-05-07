/**
 * @file RcppApi_wrappers.cpp
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

#include <RcppApi_wrappers.h>

using QSM = arbor::qsm::QSM;
using QSF = arbor::qsm::QSF;
using QSMEdge = arbor::qsm::QSMEdge;
using QSMNode = arbor::qsm::QSMNode;

QSF as_qsf(Rcpp::List x)
{
  QSF qsf;

  Rcpp::CharacterVector names = x.names();

  for (int i = 0; i < x.size(); ++i)
  {
    std::string name = Rcpp::as<std::string>(names[i]);
    Rcpp::DataFrame df = Rcpp::as<Rcpp::DataFrame>(x[i]);
    QSM qsm = as_qsm(df);
    qsf.add_qsm(name, qsm);
  }

  return qsf;
}

Rcpp::DataFrame as_dataframe(const PointCloud& cloud)
{
  size_t n = cloud.size();

  Rcpp::NumericVector x(n);
  Rcpp::NumericVector y(n);
  Rcpp::NumericVector z(n);

  for (size_t i = 0; i < n; ++i)
  {
    x[i] = cloud.get_x(i);
    y[i] = cloud.get_y(i);
    z[i] = cloud.get_z(i);
  }

  Rcpp::List df_list = Rcpp::List::create(
    Rcpp::_["X"] = x,
    Rcpp::_["Y"] = y,
    Rcpp::_["Z"] = z
  );

  if (cloud.has_treeid())
  {
    Rcpp::IntegerVector treeid(n);
    for (size_t i = 0; i < n; ++i) treeid[i] = cloud.get_treeid(i);
    df_list["treeID"] = treeid;
  }

  if (cloud.has_foliage())
  {
    Rcpp::IntegerVector foliage(n);
    for (size_t i = 0; i < n; ++i) foliage[i] = cloud.get_foliage(i);
    df_list["foliage"] = foliage;
  }

  if (cloud.has_hag())
  {
    Rcpp::NumericVector hag(n);
    for (size_t i = 0; i < n; ++i) hag[i] = cloud.get_hag(i);
    df_list["hag"] = hag;
  }

  if (cloud.has_pwood())
  {
    Rcpp::NumericVector pwood(n);
    for (size_t i = 0; i < n; ++i) pwood[i] = cloud.get_pwood(i);
    df_list["pwood"] = pwood;
  }

  if (cloud.has_passage())
  {
    Rcpp::IntegerVector passage(n);
    for (size_t i = 0; i < n; ++i) passage[i] = cloud.get_passage(i);
    df_list["passage"] = passage;
  }

  Rcpp::DataFrame df(df_list);
  return df;
}



QSM as_qsm(Rcpp::DataFrame df)
{
  // --- mandatory columns ---
  const char* required[] = {
    "cyl_ID", "parent_ID", "startX", "startY", "startZ", "endX", "endY", "endZ"
  };
  for (const auto& col : required) {
    if (!df.containsElementNamed(col)) {
      Rcpp::stop("DataFrame must contain: cyl_ID, parent_ID, startX, startY, startZ, endX, endY, endZ");
    }
  }

  Rcpp::IntegerVector  cid = df["cyl_ID"];
  Rcpp::IntegerVector  pid = df["parent_ID"];
  Rcpp::NumericVector  sx  = df["startX"];
  Rcpp::NumericVector  sy  = df["startY"];
  Rcpp::NumericVector  sz  = df["startZ"];
  Rcpp::NumericVector  ex  = df["endX"];
  Rcpp::NumericVector  ey  = df["endY"];
  Rcpp::NumericVector  ez  = df["endZ"];

  // optional columns
  Rcpp::NumericVector radius         = df.containsElementNamed("radius")          ? df["radius"]          : Rcpp::NumericVector(cid.size(), arbor::qsm::RADIUS_UNSET);
  Rcpp::NumericVector theoric_radius = df.containsElementNamed("theoric_radius")  ? df["theoric_radius"]  : Rcpp::NumericVector(cid.size(), arbor::qsm::RADIUS_UNSET);
  Rcpp::IntegerVector axis_ID        = df.containsElementNamed("axis_ID")         ? df["axis_ID"]         : Rcpp::IntegerVector(cid.size(), 0);
  Rcpp::IntegerVector branch_order   = df.containsElementNamed("branch_order")    ? df["branch_order"]    : Rcpp::IntegerVector(cid.size(), 0);
  Rcpp::NumericVector dist_to_root   = df.containsElementNamed("dist_to_root")    ? df["dist_to_root"]    : Rcpp::NumericVector(cid.size(), arbor::qsm::DISTANCE_TO_ROOT_UNSET);
  Rcpp::NumericVector subtree_length = df.containsElementNamed("subtree_length")  ? df["subtree_length"]  : Rcpp::NumericVector(cid.size(), arbor::qsm::SUBTREE_LENGTH_UNSET);

  int n = cid.size();

  // Prepare the QSM
  QSM graph;

  if (df.hasAttribute("message"))
  {
    graph.messages = Rcpp::as<std::vector<std::string>>(df.attr("message"));
  }

  // Map coordinates to node IDs
  constexpr int digits = 6;
  const double  factor = std::pow(10.0, digits);

  struct CoordKey {
    int x, y, z;
    bool operator==(const CoordKey& o) const noexcept { return x == o.x && y == o.y && z == o.z; }
  };
  struct CoordKeyHash {
    std::size_t operator()(const CoordKey& k) const noexcept {
      std::size_t h1 = std::hash<int>{}(k.x);
      std::size_t h2 = std::hash<int>{}(k.y);
      std::size_t h3 = std::hash<int>{}(k.z);
      return h1 ^ (h2 << 1) ^ (h3 << 2);
    }
  };
  std::unordered_map<CoordKey, int, CoordKeyHash> coord_to_node;
  auto get_key = [&](double x, double y, double z) -> CoordKey {
    return CoordKey{
      static_cast<int>(std::llround(x * factor)),
      static_cast<int>(std::llround(y * factor)),
      static_cast<int>(std::llround(z * factor))
    };
  };
  auto get_or_create_node = [&](double x, double y, double z) -> int {
    auto key = get_key(x, y, z);
    auto it  = coord_to_node.find(key);
    if (it != coord_to_node.end()) return it->second;
    int nid = graph.add_node({x, y, z});
    coord_to_node[key] = nid;
    return nid;
  };

  for (int i = 0; i < n; ++i)
  {
    int src = get_or_create_node(sx[i], sy[i], sz[i]);
    int tgt = get_or_create_node(ex[i], ey[i], ez[i]);

    QSMEdge ed;
    ed.cyl_ID            = cid[i];
    ed.parent_ID         = pid[i];
    ed.radius            = Rcpp::NumericVector::is_na(radius[i])         ? arbor::qsm::RADIUS_UNSET          : radius[i];
    ed.conic_allometry   = Rcpp::NumericVector::is_na(theoric_radius[i]) ? arbor::qsm::RADIUS_UNSET          : theoric_radius[i];
    ed.axis_ID           = Rcpp::IntegerVector::is_na(axis_ID[i])        ? 0                                 : axis_ID[i];
    ed.branch_order      = Rcpp::IntegerVector::is_na(branch_order[i])   ? 0                                 : branch_order[i];
    ed.distance_to_root  = Rcpp::NumericVector::is_na(dist_to_root[i])   ? arbor::qsm::DISTANCE_TO_ROOT_UNSET: dist_to_root[i];
    ed.subtree_length    = Rcpp::NumericVector::is_na(subtree_length[i]) ? arbor::qsm::SUBTREE_LENGTH_UNSET  : subtree_length[i];
    ed.subtree_max_endZ  = arbor::qsm::SUBTREE_MAXZ_UNSET;
    ed.subtree_volume    = arbor::qsm::SUBTREE_VOLUME_UNSET;

    graph.add_edge(src, tgt, ed);
  }

  return graph;
}

Rcpp::DataFrame as_dataframe(const QSM& graph)
{
  const auto& edge_map = graph.edges();
  int n = edge_map.size();

  // Copy into a vector of pointers, then sort by cyl_ID for deterministic output
  std::vector<const std::pair<const int, QSM::EdgeInfo>*> vec;
  vec.reserve(n);
  for (const auto& kv : edge_map) vec.push_back(&kv);
  std::sort(vec.begin(), vec.end(), [](const auto* a, const auto* b)
  {
    return a->second.data.cyl_ID < b->second.data.cyl_ID;
  });

  // Allocate R vectors
  Rcpp::IntegerVector cid(n), pid(n), axis_id(n), branch_order(n);
  Rcpp::NumericVector sx(n), sy(n), sz(n), ex(n), ey(n), ez(n);
  Rcpp::NumericVector radius(n), subtree_length(n), dist_to_root(n);

  // Fill the vectors row by row
  for (int i = 0; i < n; i++)
  {
    const auto& einfo = vec[i]->second;

    // Node indices
    const QSMNode& src = graph.node(einfo.source);
    const QSMNode& tgt = graph.node(einfo.target);
    const QSMEdge& ed  = einfo.data;

    cid[i]            = ed.cyl_ID;
    pid[i]            = ed.parent_ID;
    sx[i]             = src.x;
    sy[i]             = src.y;
    sz[i]             = src.z;
    ex[i]             = tgt.x;
    ey[i]             = tgt.y;
    ez[i]             = tgt.z;
    radius[i]         = (ed.radius == arbor::qsm::RADIUS_UNSET)                 ? NA_REAL : ed.radius;
    dist_to_root[i]   = (ed.distance_to_root == arbor::qsm::DISTANCE_TO_ROOT_UNSET) ? NA_REAL : ed.distance_to_root;
    subtree_length[i] = (ed.subtree_length == arbor::qsm::SUBTREE_LENGTH_UNSET) ? NA_REAL : ed.subtree_length;
    axis_id[i]        = ed.axis_ID;
    branch_order[i]   = ed.branch_order;
  }

  Rcpp::DataFrame df = Rcpp::DataFrame::create(
    Rcpp::Named("startX") = sx,
    Rcpp::Named("startY") = sy,
    Rcpp::Named("startZ") = sz,
    Rcpp::Named("endX")   = ex,
    Rcpp::Named("endY")   = ey,
    Rcpp::Named("endZ")   = ez,
    Rcpp::Named("cyl_ID") = cid,
    Rcpp::Named("parent_ID") = pid,
    Rcpp::Named("axis_ID") = axis_id,
    Rcpp::Named("branch_order") = branch_order,
    Rcpp::Named("radius") = radius,
    Rcpp::Named("dist_to_root") = dist_to_root,
    Rcpp::Named("subtree_length") = subtree_length,
    Rcpp::Named("stringsAsFactors") = false
  );

  df.attr("message") = graph.messages;

  return df;
}

#endif
