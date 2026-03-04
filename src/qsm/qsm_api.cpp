#include "arbor.h"
#include "QSMbuilder.h"

namespace arbor::qsm {

QSM qsm(const PointCloud& tree, const settings::ArborParameters& params, const Logger& logger)
{
  QSM qsm;
  QSMbuilder builder(qsm, params);
  builder.build(tree);
  return qsm;
}

}
