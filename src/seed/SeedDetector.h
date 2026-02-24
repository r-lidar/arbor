#ifndef SEEDDETECTOR_H
#define SEEDDETECTOR_H

#include "api.h"

namespace SeedDetectorGeometries
{
  struct Circle
  {
    double X;
    double Y;
    double Z;
    double R;
    int id;

    Circle(double x, double y, double z, double r, int i) : X(x), Y(y), Z(z), R(r), id(i) {}
  };

  struct Point3D
  {
    double X;
    double Y;
    double Z;

    Point3D(double x, double y, double z) : X(x), Y(y), Z(z) {}
  };
}

class SeedDetector
{
public:
  SeedDetector(const ArborParameters& par) : params(par) {};
  void set_logger(Logger new_logger) { logger = std::move(new_logger); }
  void run(const PointCloud& scene);
  static std::vector<SeedDetectorGeometries::Circle> detect_tree_circles(const PointCloud& wood, double resolution = 0.05, int connectivity = 26, int num_ransac_iterations = 400, double inlier_threshold = 0.02, size_t min_cluster_size = 10);
  static std::vector<SeedDetectorGeometries::Point3D> generate_cage(const std::vector<SeedDetectorGeometries::Circle>& circles, double decimation);

  const PointCloud& get_long_passages() const { return long_passages; }
  const PointCloud& get_short_passages() const { return short_passages; }
  const PointCloud& get_wood() const { return wood; }
  const PointCloud& get_cages() const { return cages; }
  const PointCloud& get_primary_seeds() const { return primary_seeds; }
  const PointCloud& get_secondary_seeds() const { return secondary_seeds; }
  const PointCloud& get_seeds() const { return seeds; }
  PointCloud&& move_seeds() { return std::move(seeds); }   // Bonus: allow to move the data out (transfer ownership)

private:
  void compute_min_hag(const PointCloud& scene);
  void slice_wood(const PointCloud& scene);
  void extract_passages(const PointCloud& scene);
  void make_cages();
  void safe_zone();
  void find_primary_seeds();
  void find_secondary_seeds();
  void merge_short_passages();
  void filter_seeds();
  PointCloud densify_passages(const PointCloud& x);


private:
  std::vector<SeedDetectorGeometries::Circle> circles;

  PointCloud long_passages;
  PointCloud short_passages;
  PointCloud wood;
  PointCloud cages;

  PointCloud primary_seeds;
  PointCloud secondary_seeds;
  PointCloud seeds;

  double min_hag;

  ArborParameters params;
  Logger logger = [](const std::string&) {};
};

#endif
