#ifndef API_H
#define API_H

#include "Adaptor.h"

#include <vector>

struct WoodlikelihoodParameters
{
  int k = 80;
};

struct GraphParameters
{
  bool downward = false;
  int k = 10;
  int k_seed = 100;
  double decimation = 0.05;
  double space_res = 0.2;
  double max_gap = 0.2;
  double power = 3.0;
  double wood2wood = 0.1;
  double leaf2leaf = 20.0;
  double wood2leaf = 100.0;
  std::vector<float> angle_penalty;
  GraphParameters(): angle_penalty(181)
  {
    for (int x = 0; x <= 180; ++x)
    {
      float y = std::exp(0.046051f * x);
      angle_penalty[x] = (x > 100) ? 100.0f : y;
    }
  }
};

struct SemanticParameters
{
  int   min_passage = 3;
  double high_pwood_threshold   = 0.9;
  double medium_pwood_threshold = 0.75;
  double connected_components_res = 0.05;
  int    connected_components_min = 2000;
  int    wood_assignation_k = 50;
  double wood_assignation_dist = 0.05f;
  int    wood_extra_reasignation_k = 10;
  double wood_extra_reasignation_dist = 0.03;
  int    medium_pwood_sor_k = 50;
  double medium_pwood_sor_m = 0.05;
  double ground_res = 0.2;
};

struct SeedParameters
{
  std::vector<double> slice_at = {0.7, 0.9};
  double slice_thickness = 0.05;
  int min_passage = 15;
  double safe_zone = 0.2;
};

struct ArborParameters
{
  WoodlikelihoodParameters woodlikelihood;
  GraphParameters pathfinder;
  SemanticParameters semantic;
  SeedParameters seeds;
};

using Logger = std::function<void(const std::string&)>;

namespace arbor
{
  namespace segment
  {
    void segment_instance(PointCloud& core, const PointCloud& seeds, const ArborParameters& params, const Logger& logger = [](const std::string&) {});
    void segment_semantic(PointCloud& core, const PointCloud& ground, const ArborParameters& params, const Logger& logger = [](const std::string&) {});
  }
  namespace seeds
  {
    PointCloud find_seeds(const PointCloud&,  const ArborParameters& params, const Logger& logger = [](const std::string&) {});
  }
  namespace utils
  {
    std::vector<bool> homogeneization(const PointCloud& pc, double res, bool hybrid = true);
    std::vector<bool> sor(const PointCloud& pc, unsigned int k, double m, int ncpu);
    std::vector<float> anisotropy(const PointCloud& pc, int k, int ncpu);
  }
}

#endif
