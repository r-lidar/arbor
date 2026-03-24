#include "SeedDetector.h"

namespace arbor::seeds {

PointCloud find_seeds(const PointCloud& scene, const settings::ArborParameters& params)
{
  arbor::seeds::SeedDetector sd(params);
  sd.run(scene);
  return sd.move_seeds();
}

}
