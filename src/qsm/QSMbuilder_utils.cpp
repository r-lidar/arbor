/**
 * @file QSMbuilder_utils.cpp
 * Project: Arbor
 * 
 * Copyright (C) 2026 Jean-Romain Roussel (r-lidar) <info @ r-lidar.com>
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

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
