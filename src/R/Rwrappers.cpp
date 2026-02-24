#include <Rwrappers.h>

QSM as_qsm(Rcpp::DataFrame df)
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
  Rcpp::NumericVector radius = df.containsElementNamed("radius") ? df["radius"] : Rcpp::NumericVector(cid.size(), RADIUS_UNSET);
  Rcpp::NumericVector theoric_radius = df.containsElementNamed("theoric_radius") ? df["theoric_radius"] : Rcpp::NumericVector(cid.size(), RADIUS_UNSET);
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

    // Optional
    c.radius         = Rcpp::NumericVector::is_na(radius[i])        ? RADIUS_UNSET         : radius[i];
    c.conic_allometry= Rcpp::NumericVector::is_na(theoric_radius[i])? RADIUS_UNSET         : theoric_radius[i];
    c.axis_ID        = Rcpp::IntegerVector::is_na(axis_ID[i])       ? 0                    : axis_ID[i];
    c.branch_order   = Rcpp::IntegerVector::is_na(branch_order[i])  ? 0                    : branch_order[i];
    c.subtree_length = Rcpp::NumericVector::is_na(subtree_length[i])? SUBTREE_LENGTH_UNSET : subtree_length[i];
    c.subtree_max_endZ = SUBTREE_MAXZ_UNSET;
    qsm.add_cylinder(c);
  }

  return qsm;
}

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

Rcpp::DataFrame as_dataframe(const QSM& qsm)
{
  const auto& cyl_map = qsm.cylinders();
  int n = cyl_map.size();

  // Copy into a vector of pointers (or references), then sort by cyl_ID
  std::vector<const QSMcylinder*> vec;
  vec.reserve(n);
  for (const auto& kv : cyl_map) vec.push_back(&kv.second);
  std::sort(vec.begin(), vec.end(), [](const QSMcylinder* a, const QSMcylinder* b)
  {
    return a->cyl_ID < b->cyl_ID;
  });

  // Allocate R vectors
  Rcpp::IntegerVector cid(n);
  Rcpp::IntegerVector pid(n);
  Rcpp::IntegerVector axis_id(n), branch_order(n);
  Rcpp::NumericVector sx(n), sy(n), sz(n);
  Rcpp::NumericVector ex(n), ey(n), ez(n);
  Rcpp::NumericVector radius(n);
  Rcpp::NumericVector subtree_length(n);

  // Fill the vectors row by row
  for (int i = 0; i < n; i++)
  {
    const QSMcylinder* c = vec[i];

    cid[i]            = c->cyl_ID;
    pid[i]            = c->parent_ID;
    sx[i]             = c->startX;
    sy[i]             = c->startY;
    sz[i]             = c->startZ;
    ex[i]             = c->endX;
    ey[i]             = c->endY;
    ez[i]             = c->endZ;
    radius[i]         = (c->radius == RADIUS_UNSET)               ? NA_REAL : c->radius;
    subtree_length[i] = (c->subtree_length == SUBTREE_LENGTH_UNSET) ? NA_REAL : c->subtree_length;
    axis_id[i]        = c->axis_ID;
    branch_order[i]   = c->branch_order;
  }

  //Build DataFrame
  return Rcpp::DataFrame::create(
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
    Rcpp::Named("subtree_length") = subtree_length,
    Rcpp::Named("stringsAsFactors") = false
  );
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
