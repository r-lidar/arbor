#include <SpatialIndex.h>
using namespace Rcpp;
using namespace lidR;

#include <Rcpp.h>
#include <progress.hpp>     // Include Progress library
#include <progress_bar.hpp> // Include Progress bar
#include <omp.h>


#define ROUNDANY(x,m) round((x) / m) * m

// [[Rcpp::export]]
DataFrame cpp_smooth3d(S4 las, NumericVector radius, NumericVector weight, int ncpu, bool pgbar, bool verbose)
{
  SparsePartition3D tree(las);

  DataFrame data = as<DataFrame>(las.slot("data"));
  NumericVector X = data["X"];
  NumericVector Y = data["Y"];
  NumericVector Z = data["Z"];
  NumericVector Xsmooth = clone(X);
  NumericVector Ysmooth = clone(Y);
  NumericVector Zsmooth = clone(Z);
  int npoints = X.size();

  if (weight.size() != X.size()) stop("Invadid size for weights");
  if (radius.size() != X.size()) stop("Invadid size for radius");

  // Initialize the progress bar
  Progress pb(npoints, true);

  bool abort = false;

  #pragma omp parallel for num_threads(ncpu)
  for (unsigned int i = 0; i < npoints; i++)
  {
    if (Progress::check_abort())
    {
      abort = true;
      continue;
    }

    std::vector<PointXYZ> pts;

    Sphere s(X[i], Y[i], Z[i], radius[i]);
    tree.lookup(s, pts);

    double xtot = 0;
    double ytot = 0;
    double ztot = 0;
    double wtot = 0;

    for (unsigned int j = 0; j < pts.size(); j++)
    {
      double dx = X[i] - pts[j].x;
      double dy = Y[i] - pts[j].y;
      double dz = Z[i] - pts[j].z;
      double w = weight[i];

      xtot += w * pts[j].x;
      ytot += w * pts[j].y;
      ztot += w * pts[j].z;
      wtot += w;
    }

    Xsmooth[i] = xtot / wtot;
    Ysmooth[i] = ytot / wtot;
    Zsmooth[i] = ztot / wtot;

    if (pgbar) pb.increment();
  }

  if (abort) throw Rcpp::internal::InterruptedException();

  return DataFrame::create(Named("X") = Xsmooth, Named("Y") = Ysmooth, Named("Z") = Zsmooth);
}
