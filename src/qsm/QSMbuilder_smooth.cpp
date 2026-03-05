#include "QSMbuilder.h"

#include <algorithm>

namespace arbor::qsm {

// Compute distance from point b to line ac
static double dist2line(const std::array<double, 3>& b, const std::array<double, 3>& a, const std::array<double, 3>& c)
{
  std::array<double, 3> ab, ac;
  for (int i = 0; i < 3; ++i)
  {
    ab[i] = b[i] - a[i];
    ac[i] = c[i] - a[i];
  }

  double t_num = 0.0, t_den = 0.0;
  for (int i = 0; i < 3; ++i)
  {
    t_num += ab[i] * ac[i];
    t_den += ac[i] * ac[i];
  }

  double t = t_num / t_den;
  std::array<double, 3> proj;
  for (int i = 0; i < 3; ++i)
  {
    proj[i] = a[i] + t * ac[i];
  }

  double d2 = 0.0;
  for (int i = 0; i < 3; ++i)
  {
    double diff = b[i] - proj[i];
    d2 += diff * diff;
  }

  return std::sqrt(d2);
}

void QSMbuilder::smooth_skeleton(int niter, double th)
{
  if (graph.edge_count() == 0) return;

  logger("Smoothing skeleton");

  // Build axis map: axis_ID -> ordered list of edge IDs (sorted by cyl_ID for root→tip order)
  std::unordered_map<int, std::vector<int>> axis_map;
  for (const auto& [eid, einfo] : graph.edges())
    axis_map[einfo.data.axis_ID].push_back(eid);

  for (auto& [axis, vec] : axis_map)
    std::sort(vec.begin(), vec.end(), [this](int a, int b) {
      return graph.edge_data(a).cyl_ID < graph.edge_data(b).cyl_ID;
    });

  // Iterative smoothing
  for (int iter = 0; iter < niter; ++iter)
  {
    for (const auto& [axis_id, edge_ids] : axis_map)
    {
      if (edge_ids.size() < 2) continue;

      for (size_t j = 1; j < edge_ids.size(); ++j)
      {
        int prev_eid = edge_ids[j - 1];
        int curr_eid = edge_ids[j];

        QSMGraph::NodeID prev_src = graph.edge(prev_eid).source;
        QSMGraph::NodeID prev_tgt = graph.edge(prev_eid).target; // == curr_src
        QSMGraph::NodeID curr_tgt = graph.edge(curr_eid).target;

        const QSMNode& prev_src_n = graph.node(prev_src);
        const QSMNode& prev_tgt_n = graph.node(prev_tgt);
        const QSMNode& curr_tgt_n = graph.node(curr_tgt);

        std::array<double, 3> a = {curr_tgt_n.x,  curr_tgt_n.y,  curr_tgt_n.z};
        std::array<double, 3> b = {prev_src_n.x,  prev_src_n.y,  prev_src_n.z};
        std::array<double, 3> c = {prev_tgt_n.x,  prev_tgt_n.y,  prev_tgt_n.z};

        double d = dist2line(b, a, c);
        if (d <= th) continue;

        // Move the shared junction node to the midpoint
        std::array<double, 3> mid = {
          0.5 * (prev_src_n.x + curr_tgt_n.x),
          0.5 * (prev_src_n.y + curr_tgt_n.y),
          0.5 * (prev_src_n.z + curr_tgt_n.z)
        };

        // Updating the node automatically propagates to ALL edges that reference it
        // (both prev_eid's target and curr_eid's source, and any sibling branches)
        QSMNode& junction = graph.node(prev_tgt);
        junction.x = mid[0];
        junction.y = mid[1];
        junction.z = mid[2];
      }
    }
  }
}

}
