#ifndef QSMCONVERSION_H
#define QSMCONVERSION_H

#include <cmath>
#include <algorithm>
#include <vector>
#include <unordered_map>

#include "QSM.h"
#include "QSMGraph.h"

namespace arbor::qsm {

// Convert a legacy QSM (edge-centric, with cyl_ID / parent_ID) to a
// QSMGraph (explicit nodes + edges).
//
// Unique 3-D junction nodes are extracted from cylinder endpoints by
// coordinate matching.  Cylinders are processed in ascending cyl_ID
// order so that the graph edge IDs are monotonically increasing and
// predictable.
//
// The original cyl_ID and parent_ID are preserved in QSMEdge for
// lossless round-trip conversion.
inline QSMGraph qsm_to_graph(const QSM& qsm)
{
  QSMGraph graph;

  // ---- Build a coordinate → NodeID lookup ----
  constexpr int    digits = 6;
  const double     factor = std::pow(10.0, digits);

  struct CoordKey {
    int x, y, z;
    bool operator==(const CoordKey& o) const noexcept
    { return x == o.x && y == o.y && z == o.z; }
  };
  struct CoordKeyHash {
    std::size_t operator()(const CoordKey& k) const noexcept
    {
      std::size_t h1 = std::hash<int>{}(k.x);
      std::size_t h2 = std::hash<int>{}(k.y);
      std::size_t h3 = std::hash<int>{}(k.z);
      return h1 ^ (h2 << 1) ^ (h3 << 2);
    }
  };

  std::unordered_map<CoordKey, int, CoordKeyHash> coord_to_node;
  coord_to_node.reserve(qsm.size() * 2);

  auto get_key = [&](double x, double y, double z) -> CoordKey {
    return CoordKey{
      static_cast<int>(std::llround(x * factor)),
      static_cast<int>(std::llround(y * factor)),
      static_cast<int>(std::llround(z * factor))
    };
  };

  auto get_or_create_node = [&](double x, double y, double z) -> int {
    auto key = get_key(x, y, z);
    auto it  = coord_to_node.find(key);
    if (it != coord_to_node.end()) return it->second;
    int nid = graph.add_node({x, y, z});
    coord_to_node[key] = nid;
    return nid;
  };

  // ---- Sort cylinders by cyl_ID for deterministic edge ordering ----
  std::vector<const QSMcylinder*> sorted;
  sorted.reserve(qsm.size());
  for (const auto& kv : qsm.cylinders()) sorted.push_back(&kv.second);
  std::sort(sorted.begin(), sorted.end(),
    [](const QSMcylinder* a, const QSMcylinder* b) { return a->cyl_ID < b->cyl_ID; });

  // ---- Create one graph edge per cylinder ----
  for (const QSMcylinder* cyl : sorted)
  {
    int src = get_or_create_node(cyl->startX, cyl->startY, cyl->startZ);
    int tgt = get_or_create_node(cyl->endX,   cyl->endY,   cyl->endZ);

    QSMEdge ed;
    ed.cyl_ID           = cyl->cyl_ID;
    ed.parent_ID        = cyl->parent_ID;
    ed.radius           = cyl->radius;
    ed.conic_allometry  = cyl->conic_allometry;
    ed.subtree_length   = cyl->subtree_length;
    ed.subtree_max_endZ = cyl->subtree_max_endZ;
    ed.subtree_volume   = cyl->subtree_volume;
    ed.axis_ID          = cyl->axis_ID;
    ed.branch_order     = cyl->branch_order;

    graph.add_edge(src, tgt, ed);
  }

  return graph;
}

// Convert a QSMGraph back to a legacy QSM.
//
// cyl_ID and parent_ID are read directly from QSMEdge (populated either
// by qsm_to_graph or by the QSMbuilder during construction).
inline QSM graph_to_qsm(const QSMGraph& graph)
{
  QSM qsm;

  for (const auto& [eid, einfo] : graph.edges())
  {
    const QSMNode& src = graph.node(einfo.source);
    const QSMNode& tgt = graph.node(einfo.target);
    const QSMEdge& ed  = einfo.data;

    QSMcylinder cyl;
    cyl.cyl_ID          = ed.cyl_ID  != 0 ? ed.cyl_ID  : eid;
    cyl.parent_ID       = ed.parent_ID;
    cyl.startX          = src.x;
    cyl.startY          = src.y;
    cyl.startZ          = src.z;
    cyl.endX            = tgt.x;
    cyl.endY            = tgt.y;
    cyl.endZ            = tgt.z;
    cyl.radius          = ed.radius;
    cyl.conic_allometry = ed.conic_allometry;
    cyl.subtree_length  = ed.subtree_length;
    cyl.subtree_max_endZ = ed.subtree_max_endZ;
    cyl.subtree_volume  = ed.subtree_volume;
    cyl.axis_ID         = ed.axis_ID;
    cyl.branch_order    = ed.branch_order;

    qsm.add_cylinder(cyl);
  }

  return qsm;
}

} // namespace arbor::qsm

#endif // QSMCONVERSION_H
