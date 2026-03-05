#include "arbor.h"
#include "QSMbuilder.h"
#include "QSMConversion.h"

namespace arbor::qsm {

QSM qsm(const PointCloud& tree, const settings::ArborParameters& params, const Logger& logger)
{
  QSMGraph graph;
  QSMbuilder builder(graph, params);
  builder.set_logger(logger);
  builder.build(tree);
  return graph_to_qsm(graph);
}

}
