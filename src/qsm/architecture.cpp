#include <vector>
#include <string>
#include <tuple>
#include <stdexcept>

#include "QSM.h"

struct Result
{
  // in same order as input rows
  std::vector<int> cyl_ID;
  std::vector<double> subtree_length;
  std::vector<double> subtree_max_endZ;
  std::vector<int> axis_ID;
  std::vector<int> branching_order;
};

Result compute_architecture(const std::vector<int>& cyl_ID,
                            const std::vector<int>& parent_ID,
                            const std::vector<double>& length,
                            const std::vector<double>& startZ,
                            const std::vector<double>& endZ,
                            int root_id)
{
  if (cyl_ID.size() != parent_ID.size() ||
      cyl_ID.size() != length.size() ||
      cyl_ID.size() != startZ.size() ||
      cyl_ID.size() != endZ.size())
  {
    throw std::invalid_argument("Input vectors must all have the same length");
  }

  QSM qsm;
  qsm.build_from_vectors(cyl_ID, parent_ID, length, startZ, endZ);
  qsm.compute_architecture(root_id);

  // prepare result in input order
  Result res;
  size_t n = cyl_ID.size();
  res.cyl_ID.resize(n);
  res.subtree_length.resize(n);
  res.subtree_max_endZ.resize(n);
  res.axis_ID.resize(n);
  res.branching_order.resize(n);

  const auto& subtree_lengths = qsm.get_subtree_lengths();
  const auto& subtree_max_z = qsm.get_subtree_max_z();
  const auto& subtree_ids = qsm.get_subtree_ids();
  const auto& branching_orders = qsm.get_branching_orders();

  for (size_t i = 0; i < n; ++i)
  {
    int cid = cyl_ID[i];
    res.cyl_ID[i] = cid;

    if (subtree_lengths.count(cid))
      res.subtree_length[i] = subtree_lengths.at(cid);
    else
      res.subtree_length[i] = 0.0;

    if (subtree_max_z.count(cid))
      res.subtree_max_endZ[i] = subtree_max_z.at(cid);
    else
      res.subtree_max_endZ[i] = endZ[i]; // fallback to own endZ

    if (subtree_ids.count(cid))
      res.axis_ID[i] = subtree_ids.at(cid);
    else
      res.axis_ID[i] = -1;

    if (branching_orders.count(cid))
      res.branching_order[i] = branching_orders.at(cid);
    else
      res.branching_order[i] = -1;
  }

  return res;
}


// Rcpp wrapper

#include <Rcpp.h>

using namespace Rcpp;

DataFrame cpp_compute_architecture(DataFrame qsm, int root_id = 1)
{
  if (!qsm.containsElementNamed("cyl_ID") || !qsm.containsElementNamed("parent_ID") ||
      !qsm.containsElementNamed("length") || !qsm.containsElementNamed("startZ") ||
      !qsm.containsElementNamed("endZ"))
  {
    stop("qsm must contain columns: cyl_ID, parent_ID, length, startZ, endZ");
  }

  IntegerVector cyl_ID_r     = qsm["cyl_ID"];
  IntegerVector parent_ID_r  = qsm["parent_ID"];
  NumericVector length_r     = qsm["length"];
  NumericVector startZ_r     = qsm["startZ"];
  NumericVector endZ_r       = qsm["endZ"];

  std::vector<int> cyl_ID    = as<std::vector<int>>(cyl_ID_r);
  std::vector<int> parent_ID = as<std::vector<int>>(parent_ID_r);
  std::vector<double> length = as<std::vector<double>>(length_r);
  std::vector<double> startZ = as<std::vector<double>>(startZ_r);
  std::vector<double> endZ   = as<std::vector<double>>(endZ_r);

  Result res;
  try
  {
    res = compute_architecture(cyl_ID, parent_ID, length, startZ, endZ, root_id);
  }
  catch (std::exception &ex)
  {
    stop(ex.what());
  }

  return DataFrame::create(
    Named("cyl_ID") = wrap(res.cyl_ID),
    Named("subtree_length") = wrap(res.subtree_length),
    Named("axis_ID") = wrap(res.axis_ID),
    Named("branching_order") = wrap(res.branching_order)
  );
}

