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

// Convert QSM to Rcpp DataFrame
inline Rcpp::DataFrame as_dataframe(const QSM& qsm)
{
  const auto& cyl_map = qsm.cylinders();
  int n = cyl_map.size();

  Rcpp::IntegerVector cid(n);
  Rcpp::IntegerVector pid(n);
  Rcpp::NumericVector sx(n), sy(n), sz(n);
  Rcpp::NumericVector ex(n), ey(n), ez(n);
  Rcpp::NumericVector radius(n);
  Rcpp::NumericVector subtree_length(n);
  Rcpp::IntegerVector axis_id(n), branch_order(n);

  for (const auto& [i, c] : cyl_map)
  {
    cid[i]            = c.cyl_ID;
    pid[i]            = c.parent_ID;
    sx[i]             = c.startX;
    sy[i]             = c.startY;
    sz[i]             = c.startZ;
    ex[i]             = c.endX;
    ey[i]             = c.endY;
    ez[i]             = c.endZ;
    radius[i]         = c.radius;
    subtree_length[i] = c.subtree_length;
    axis_id[i]        = c.axis_ID;
    branch_order[i]   = c.branch_order;
  }

  return Rcpp::DataFrame::create(
    Rcpp::Named("cyl_ID") = cid,
    Rcpp::Named("parent_ID") = pid,
    Rcpp::Named("startX") = sx,
    Rcpp::Named("startY") = sy,
    Rcpp::Named("startZ") = sz,
    Rcpp::Named("endX") = ex,
    Rcpp::Named("endY") = ey,
    Rcpp::Named("endZ") = ez,
    Rcpp::Named("radius") = radius,
    Rcpp::Named("subtree_length") = subtree_length,
    Rcpp::Named("axis_ID") = axis_id,
    Rcpp::Named("branch_order") = branch_order,
    Rcpp::Named("stringsAsFactors") = false
  );
}

#endif
