#include <vector>
#include <stdexcept>
#include <Rcpp.h>
#include "QSM.h"

static inline QSM as_qsm(Rcpp::DataFrame df)
{
  // --- mandatory columns ---
  const char* required[] = {"cyl_ID", "parent_ID", "startX", "startY", "startZ", "endX", "endY", "endZ"};
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
  Rcpp::NumericVector radius = df.containsElementNamed("radius") ? df["radius"] : Rcpp::NumericVector(cid.size(), 0.0);
  Rcpp::IntegerVector axis_ID = df.containsElementNamed("axis_ID") ? df["axis_ID"] : Rcpp::IntegerVector(cid.size(), 0);
  Rcpp::IntegerVector branch_order = df.containsElementNamed("branch_order") ? df["branch_order"] : Rcpp::IntegerVector(cid.size(), 0);
  Rcpp::NumericVector subtree_length = df.containsElementNamed("subtree_length") ? df["subtree_length"] : Rcpp::NumericVector(cid.size(), SUBTREE_LENGTH_UNSET);

  int n = cid.size();

  QSM qsm;
  QSMcylinder c;

  for (int i = 0; i < n; ++i)
  {
    c.cyl_ID    = cid[i];
    c.parent_ID = pid[i];
    c.startX    = sx[i];
    c.startY    = sy[i];
    c.startZ    = sz[i];
    c.endX      = ex[i];
    c.endY      = ey[i];
    c.endZ      = ez[i];

    c.radius           = radius[i];
    c.axis_ID          = axis_ID[i];
    c.branch_order     = branch_order[i];
    c.subtree_length   = subtree_length[i];
    c.subtree_max_endZ = SUBTREE_MAXZ_UNSET;
    qsm.add_cylinder(c);
  }

  return qsm;
}

static inline Rcpp::DataFrame as_dataframe(const QSM& qsm)
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

void qsm_write_cpp(Rcpp::DataFrame df, std::string filename)
{
  QSM qsm = as_qsm(df);
  qsm.write(filename, 16);
}


