#include "arbor.h"

namespace arbor::qsm {

QSM qsm(const PointCloud& pc, const settings::ArborParameters& params, const Logger& logger)
{
  auto layers = QSM::layers(pc, 0.2);
  return QSM();

}

}
