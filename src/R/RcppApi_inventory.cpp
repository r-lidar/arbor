/**
 * @file RcppApi_inventory.cpp
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

// Bindings used by the arborInventoryMerger package. The point cloud is never modified:
// functions that change tree IDs return the new treeID vector.

#ifdef USING_R

#include <Rcpp.h>

#include <cmath>
#include <limits>

#include "arbor.h"
#include "InventoryMerger.h"

// Inventory from a data.frame with columns ID_Arbre, X, Y, Z, DHP, POM
static arbor::inventory::Inventory as_inventory(Rcpp::DataFrame df)
{
  for (const char* name : {"ID_Arbre", "X", "Y", "Z", "DHP", "POM"})
  {
    if (!df.containsElementNamed(name))
      Rcpp::stop(std::string("Inventory: missing required column '") + name + "'");
  }

  Rcpp::NumericVector id  = Rcpp::as<Rcpp::NumericVector>(df["ID_Arbre"]);
  Rcpp::NumericVector x   = Rcpp::as<Rcpp::NumericVector>(df["X"]);
  Rcpp::NumericVector y   = Rcpp::as<Rcpp::NumericVector>(df["Y"]);
  Rcpp::NumericVector z   = Rcpp::as<Rcpp::NumericVector>(df["Z"]);
  Rcpp::NumericVector dhp = Rcpp::as<Rcpp::NumericVector>(df["DHP"]);
  Rcpp::NumericVector pom = Rcpp::as<Rcpp::NumericVector>(df["POM"]);

  arbor::inventory::Inventory inventory;
  inventory.reserve(id.size());
  for (R_xlen_t i = 0; i < id.size(); ++i)
  {
    if (!R_finite(id[i]) || id[i] != std::floor(id[i]) || id[i] > std::numeric_limits<int>::max())
      Rcpp::stop("Inventory: ID_Arbre must contain integral numbers (row " + std::to_string(i+1) + ")");
    inventory.add_tree(static_cast<int>(id[i]), x[i], y[i], z[i], dhp[i]/2, pom[i]);
  }

  return inventory;
}

// Shallow copy of the data.frame in which only treeID is a new vector
static Rcpp::DataFrame with_treeid_copy(Rcpp::DataFrame las)
{
  if (!las.containsElementNamed("treeID"))
    Rcpp::stop("The point cloud does not contain a 'treeID' attribute.");

  Rcpp::List copy(Rf_shallow_duplicate(las));
  int i = las.findName("treeID");
  copy[i] = Rcpp::clone(Rcpp::as<Rcpp::IntegerVector>(las[i]));
  return Rcpp::DataFrame(copy);
}

static Rcpp::IntegerVector treeid_of(Rcpp::DataFrame las)
{
  return las["treeID"];
}

static arbor::inventory::MatchingTable as_matching_table(Rcpp::DataFrame df)
{
  if (!df.containsElementNamed("treeID") || !df.containsElementNamed("targetID"))
    Rcpp::stop("matching_table must contain the columns 'treeID' and 'targetID'.");

  Rcpp::IntegerVector tree   = Rcpp::as<Rcpp::IntegerVector>(df["treeID"]);
  Rcpp::IntegerVector target = Rcpp::as<Rcpp::IntegerVector>(df["targetID"]);

  arbor::inventory::MatchingTable table;
  table.reserve(tree.size());
  for (R_xlen_t i = 0; i < tree.size(); ++i)
  {
    if (target[i] == NA_INTEGER) Rcpp::stop("matching_table$targetID must not contain NA values.");
    int id = (tree[i] == NA_INTEGER) ? 0 : tree[i];
    table.push_back({id, target[i], 0});
  }
  return table;
}

static Rcpp::DataFrame as_dataframe(const arbor::inventory::MatchingTable& table)
{
  size_t n = table.size();
  Rcpp::IntegerVector tree(n), target(n), nfit(n);
  for (size_t i = 0; i < n; ++i)
  {
    tree[i]   = table[i].tree_id > 0 ? table[i].tree_id : NA_INTEGER;
    target[i] = table[i].target_id;
    nfit[i]   = table[i].nfit;
  }

  return Rcpp::DataFrame::create(Rcpp::Named("treeID") = tree, Rcpp::Named("targetID") = target, Rcpp::Named("nfit") = nfit);
}

static Rcpp::DataFrame as_dataframe(const std::vector<arbor::utils::fitting::Vec3>& pts)
{
  size_t n = pts.size();
  Rcpp::NumericVector x(n), y(n), z(n);
  for (size_t i = 0; i < n; ++i) { x[i] = pts[i].x; y[i] = pts[i].y; z[i] = pts[i].z; }
  return Rcpp::DataFrame::create(Rcpp::Named("X") = x, Rcpp::Named("Y") = y, Rcpp::Named("Z") = z);
}

static Rcpp::IntegerVector as_ids(const std::vector<int>& ids)
{
  if (ids.empty()) return Rcpp::IntegerVector::create(NA_INTEGER);
  return Rcpp::IntegerVector(ids.begin(), ids.end());
}

Rcpp::IntegerVector aim_disentangle_ids_cpp(Rcpp::DataFrame las, Rcpp::DataFrame inventory)
{
  arbor::inventory::Inventory inventory_ = as_inventory(inventory);
  Rcpp::DataFrame copy = with_treeid_copy(las);
  PointCloud pc(copy);
  arbor::inventory::InventoryMerger(inventory_).disentangle(pc);
  return treeid_of(copy);
}

Rcpp::DataFrame aim_matching_table_cpp(Rcpp::DataFrame las, Rcpp::DataFrame inventory)
{
  arbor::inventory::Inventory inventory_ = as_inventory(inventory);
  PointCloud pc(las);
  return as_dataframe(arbor::inventory::matching_table(arbor::inventory::InventoryMerger(inventory_).match(pc)));
}

Rcpp::IntegerVector aim_reassign_ids_cpp(Rcpp::DataFrame las, Rcpp::DataFrame matching_table)
{
  arbor::inventory::MatchingTable table = as_matching_table(matching_table);
  Rcpp::DataFrame copy = with_treeid_copy(las);
  PointCloud pc(copy);
  arbor::inventory::InventoryMerger::reassign(pc, table);
  return treeid_of(copy);
}

// The matches of the matching step are not available: the inventory trees involved are fitted again
Rcpp::IntegerVector aim_resolve_multi_matching_cpp(Rcpp::DataFrame las, Rcpp::DataFrame matching_table, Rcpp::DataFrame inventory)
{
  arbor::inventory::Inventory inventory_ = as_inventory(inventory);
  arbor::inventory::MatchingTable table = as_matching_table(matching_table);
  Rcpp::DataFrame copy = with_treeid_copy(las);
  PointCloud pc(copy);
  arbor::inventory::InventoryMerger(inventory_).resolve(pc, table);
  return treeid_of(copy);
}

Rcpp::IntegerVector aim_merge_inventory_cpp(Rcpp::DataFrame las, Rcpp::DataFrame inventory)
{
  arbor::inventory::Inventory inventory_ = as_inventory(inventory);
  Rcpp::DataFrame copy = with_treeid_copy(las);
  PointCloud pc(copy);
  arbor::inventory::InventoryMerger(inventory_).run(pc);
  return treeid_of(copy);
}

Rcpp::List aim_fit_tree_cpp(Rcpp::DataFrame las, Rcpp::DataFrame tree)
{
  arbor::inventory::Inventory inventory = as_inventory(tree);
  if (inventory.size() != 1) Rcpp::stop("Cannot process multiple trees");

  const arbor::inventory::InventoryTree& t = inventory[0];

  PointCloud pc(las);
  arbor::inventory::StemMatch res = arbor::inventory::InventoryMerger(inventory).match(pc, {t.id}).front();

  Rcpp::List out;

  if (res.status != arbor::inventory::StemMatch::Status::Fitted)
  {
    out.push_back(Rcpp::List::create(
      Rcpp::Named("center_x") = t.x,
      Rcpp::Named("center_y") = t.y,
      Rcpp::Named("center_z") = t.z,
      Rcpp::Named("to_be_merge_id") = as_ids(res.tree_ids),
      Rcpp::Named("seeds") = as_dataframe(res.seeds)));
    return out;
  }

  for (const arbor::inventory::StemSlice& slice : res.slices)
  {
    const arbor::utils::fitting::FittingResult& fit = slice.fit;

    Rcpp::NumericMatrix nodes(fit.contour.size(), 3);
    for (size_t i = 0; i < fit.contour.size(); ++i)
    {
      nodes(i, 0) = fit.contour[i].x;
      nodes(i, 1) = fit.contour[i].y;
      nodes(i, 2) = fit.contour[i].z;
    }

    std::vector<arbor::utils::fitting::Vec3> seeds;
    Rcpp::IntegerVector inliers(fit.inlier_indices.size());
    for (size_t i = 0; i < fit.inlier_indices.size(); ++i)
    {
      inliers[i] = fit.inlier_indices[i] + 1;
      seeds.push_back(slice.points[fit.inlier_indices[i]]);
    }

    Rcpp::List l = Rcpp::List::create(
      Rcpp::Named("success") = fit.success,
      Rcpp::Named("radius") = fit.radius,
      Rcpp::Named("shape_type") = fit.shape_type,
      Rcpp::Named("center_x") = fit.center.x,
      Rcpp::Named("center_y") = fit.center.y,
      Rcpp::Named("center_z") = fit.center.z,
      Rcpp::Named("covered_arc_degree") = fit.arc_coverage_deg,
      Rcpp::Named("percentage_inlier") = fit.inlier_percentage,
      Rcpp::Named("nodes") = nodes,
      Rcpp::Named("inliers") = inliers,
      Rcpp::Named("offset") = slice.offset,
      Rcpp::Named("to_be_merge_id") = as_ids(slice.tree_ids),
      Rcpp::Named("seeds") = as_dataframe(seeds));

    out.push_back(l);
  }

  return out;
}

#endif
