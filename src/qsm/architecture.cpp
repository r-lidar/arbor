#include <vector>
#include <stdexcept>

#include "Rcpp_cast.h"
#include "QSM.h"

Rcpp::DataFrame qsm_architecture_cpp(Rcpp::DataFrame df, int root_id = 1)
{
  QSM qsm = as_qsm(df);
  qsm.compute_architecture(root_id); // Rcpp catches exceptions

  const auto& cylinders = qsm.cylinders();

  // Prepare new columns
  int n = (int)cylinders.size();
  Rcpp::NumericVector length(n);
  Rcpp::NumericVector subtree_length(n);
  Rcpp::IntegerVector axis_ID(n);
  Rcpp::IntegerVector branching_order(n);

  for (const auto &[i, c] : cylinders)
  {
    length[i]            = c.length();
    subtree_length[i]    = c.subtree_length;
    axis_ID[i]           = c.axis_ID;
    branching_order[i]   = c.branch_order;
  }

  df["cyl_length"]       = length;
  df["subtree_length"]   = subtree_length;
  df["axis_ID"]          = axis_ID;
  df["branch_order"]     = branching_order;

  return df;
}

