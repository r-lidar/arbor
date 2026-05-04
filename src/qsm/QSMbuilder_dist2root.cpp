#include "QSMbuilder.h"

#include <queue>

namespace arbor::qsm {

void QSMbuilder::distance_to_root()
{
  ServiceLocator::logger()("Compute distance to root");

  for (auto& [eid, einfo] : graph.edges())
    einfo.data.distance_to_root = DISTANCE_TO_ROOT_UNSET;

  NodeID root = graph.find_root_node();
  if (root < 0) return; // empty or degenerate graph

  // pair<NodeID, accumulated distance to that node's source>
  std::queue<std::pair<NodeID, double>> q;
  q.push({root, 0.0});

  while (!q.empty())
  {
    auto [nid, dist] = q.front();
    q.pop();

    const QSMNode& src_node = graph.node(nid);

    for (EdgeID eid : graph.outgoing_edges(nid))
    {
      QSMEdge& e        = graph.edge_data(eid);
      const QSMNode& tgt_node = graph.node(graph.edges().at(eid).target);

      e.distance_to_root = dist;

      double next_dist = dist + e.length(src_node, tgt_node);
      q.push({graph.edges().at(eid).target, next_dist});
    }
  }
}

}
