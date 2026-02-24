#include "SeedDetector.h"

PointCloud find_seeds(const PointCloud& scene, const ArborParameters& params, const Logger& logger)
{
  SeedDetector sd(params);
  sd.set_logger(logger);
  sd.run(scene);
  return sd.move_seeds();
}
