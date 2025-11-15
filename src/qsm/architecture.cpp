#include <vector>
#include <stdexcept>

#include "Rcpp_cast.h"
#include "QSM.h"


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

Rcpp::DataFrame qsm_architecture_cpp(Rcpp::DataFrame df, int root_id = 1)
{
  QSM qsm = as_qsm(df);
  qsm.compute_architecture(root_id); // Rcpp catches exceptions

  // Prepare new columns
  int n = (int)qsm.size();
  Rcpp::NumericVector length(n);
  Rcpp::NumericVector subtree_length(n);
  Rcpp::IntegerVector axis_ID(n);
  Rcpp::IntegerVector branching_order(n);
  for (const auto &[i, c] : qsm)
  {
    length[i-1]            = c.length();
    subtree_length[i-1]    = c.subtree_length;
    axis_ID[i-1]           = c.axis_ID;
    branching_order[i-1]   = c.branch_order;
  }
  df["cyl_length"]       = length;
  df["subtree_length"]   = subtree_length;
  df["axis_ID"]          = axis_ID;
  df["branch_order"]     = branching_order;

  return df;
}

