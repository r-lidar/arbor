#include "SeedDetector.h"

namespace arbor::seeds {

PointCloud find_seeds(const PointCloud& scene, const ArborParameters& params, const Logger& logger)
{
  arbor::seeds::SeedDetector sd(params);
  sd.set_logger(logger);
  sd.run(scene);
  return sd.move_seeds();
}

}
