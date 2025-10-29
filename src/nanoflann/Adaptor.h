#include <Rcpp.h>
#include <vector>

class PointCloudAdaptor
{
public:
  std::vector<Rcpp::NumericVector> coords;
  size_t dim;

  PointCloudAdaptor(const Rcpp::DataFrame& df, std::vector<std::string> col_names = {"X", "Y", "Z"})
  {
    dim = col_names.size();
    coords.reserve(dim);
    for (const auto& name : col_names)
      coords.push_back(df[name]);
  }

  inline size_t kdtree_get_point_count() const { return coords[0].size(); }
  inline size_t point_count() const { return coords[0].size(); }
  inline double kdtree_get_pt(const size_t idx, const size_t d) const { return coords[d][idx]; }
  template <class BBOX> bool kdtree_get_bbox(BBOX&) const { return false; }
  inline void get_point(const size_t idx, double* q) const
  {
    q[0] = kdtree_get_pt(idx, 0);
    q[1] = kdtree_get_pt(idx, 1);
    q[2] = kdtree_get_pt(idx, 2);
  }
};
