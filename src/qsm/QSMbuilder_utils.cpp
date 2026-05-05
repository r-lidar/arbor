#include "QSMbuilder.h"

namespace arbor::qsm {

void QSMbuilder::shift(double tx, double ty, double tz)
{
  ServiceLocator::logger()("Shift back to geographic coordinates");

  // Shift all node positions (a single update per node covers all incident edges)
  for (auto& [nid, ndata] : graph.nodes())
  {
    ndata.x += tx;
    ndata.y += ty;
    ndata.z += tz;
  }
}

int QSMbuilder::count_nodes_connected_to_root() const
{
  // Iterate through the public nodes map
  for (const auto& [id, data] : graph.nodes())
  {
    // The root is defined as a node with zero incoming edges
    if (graph.incoming_edges(id).empty())
    {
      // Return the count of its outgoing edges
      return graph.outgoing_edges(id).size();
    }
  }

  return 0;
}

int QSMbuilder::find_root_edge() const
{
  for (const auto& [eid, einfo] : graph.edges())
  {
    if (graph.incoming_edges(einfo.source).empty())
      return eid;
  }
  return -1;
}

}
