#ifndef RCPP_CAST_H
#define RCPP_CAST_H

#include <Rcpp.h>
#include "QSM.h"

inline QSM as_qsm(Rcpp::DataFrame df)
{
  if (!df.containsElementNamed("cyl_ID") ||
      !df.containsElementNamed("parent_ID") ||
      !df.containsElementNamed("startX") ||
      !df.containsElementNamed("startY") ||
      !df.containsElementNamed("startZ") ||
      !df.containsElementNamed("endX") ||
      !df.containsElementNamed("endY") ||
      !df.containsElementNamed("endZ"))
  {
    Rcpp::stop("DataFrame must contain: cyl_ID, parent_ID, startX, startY, startZ, endX, endY, endZ");
  }

  Rcpp::IntegerVector  cid = df["cyl_ID"];
  Rcpp::IntegerVector  pid = df["parent_ID"];
  Rcpp::NumericVector  sx  = df["startX"];
  Rcpp::NumericVector  sy  = df["startY"];
  Rcpp::NumericVector  sz  = df["startZ"];
  Rcpp::NumericVector  ex  = df["endX"];
  Rcpp::NumericVector  ey  = df["endY"];
  Rcpp::NumericVector  ez  = df["endZ"];

  int n = cid.size();
  std::vector<QSMcylinder> vec(n);

  for (int i = 0; i < n; ++i)
  {
    auto &c = vec[i];
    c.cyl_ID    = cid[i];
    c.parent_ID = pid[i];
    c.startX    = sx[i];
    c.startY    = sy[i];
    c.startZ    = sz[i];
    c.endX      = ex[i];
    c.endY      = ey[i];
    c.endZ      = ez[i];

    c.subtree_length   = SUBTREE_LENGTH_UNSET;
    c.subtree_max_endZ = SUBTREE_MAXZ_UNSET;
    c.axis_ID          = 0;
    c.branch_order     = 0;
  }

  QSM qsm;
  qsm.build_from_cylinders(vec);

  return qsm;
}

#endif
