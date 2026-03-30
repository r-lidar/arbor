#include "QSMbuilder.h"

#include <algorithm>
#include <map>

namespace arbor::qsm {

void QSMbuilder::smooth_skeleton(int steps, double lambda, double mu)
{
  if (graph.edge_count() == 0) return;

  ServiceLocator::logger()("Smoothing skeleton");

  // Build axis map: axis_ID -> ordered list of edge IDs (sorted by subtree_length for root->tip order)
  std::map<int, std::vector<int>> axis_map;
  for (const auto& [eid, einfo] : graph.edges())
  {
    axis_map[einfo.data.axis_ID].push_back(eid);
  }

  // Sort each axis by subtree_length (descending = root to tip)
  for (auto& [axis, vec] : axis_map)
  {
    std::sort(vec.begin(), vec.end(), [this](int a, int b) {
      return graph.edge_data(a).subtree_length > graph.edge_data(b).subtree_length;
    });
  }

  // Smooth each axis
  for (const auto& [axis_id, edge_ids] : axis_map)
  {
    // 'steps' controls smoothing strength. We don't want to smooth too much small diameter branches
    // We don't have diameters yet. So 'steps' is reduced by 3/4 at each branch order as a workaround
    int curr_steps = steps;
    int branch_order = graph.edge(edge_ids[0]).data.branch_order;
    curr_steps = std::floor(static_cast<float>(curr_steps) * std::pow(3/4, static_cast<float>(branch_order-1)));
    curr_steps = std::max(curr_steps, 5);

    if (edge_ids.size() < 3) continue;  // Need at least 3 edges to smooth

    int n = edge_ids.size();

    // Extract coordinates from source nodes of each edge in the axis
    // (plus the target node of the last edge to complete the chain)
    std::vector<double> x(n + 1), y(n + 1), z(n + 1);

    for (int i = 0; i < n; ++i)
    {
      int eid = edge_ids[i];
      int source_node = graph.edge(eid).source;
      const QSMNode& node = graph.node(source_node);

      x[i] = node.x;
      y[i] = node.y;
      z[i] = node.z;
    }

    // Add the target node of the last edge
    int last_eid = edge_ids[n - 1];
    int target_node = graph.edge(last_eid).target;
    const QSMNode& last_node = graph.node(target_node);
    x[n] = last_node.x;
    y[n] = last_node.y;
    z[n] = last_node.z;

    // Apply Taubin smoothing
    for (int s = 0; s < curr_steps; ++s)
    {
      // Pass 1: Shrink (Lambda)
      // Keep first point (root) fixed: start at i=1
      // Keep last point (tip) fixed: end at i < n
      for (int i = 1; i < n; ++i)
      {
        double dx = 0.5 * (x[i-1] + x[i+1]) - x[i];
        double dy = 0.5 * (y[i-1] + y[i+1]) - y[i];
        double dz = 0.5 * (z[i-1] + z[i+1]) - z[i];

        x[i] += lambda * dx;
        y[i] += lambda * dy;
        z[i] += lambda * dz;
      }

      // Pass 2: Inflate (Mu)
      for (int i = 1; i < n; ++i)
      {
        double dx = 0.5 * (x[i-1] + x[i+1]) - x[i];
        double dy = 0.5 * (y[i-1] + y[i+1]) - y[i];
        double dz = 0.5 * (z[i-1] + z[i+1]) - z[i];

        x[i] += mu * dx;
        y[i] += mu * dy;
        z[i] += mu * dz;
      }
    }

    // Write smoothed coordinates back to the graph nodes
    for (int i = 0; i < n; ++i)
    {
      int eid = edge_ids[i];
      int source_node = graph.edge(eid).source;
      QSMNode& node = graph.node(source_node);

      node.x = x[i];
      node.y = y[i];
      node.z = z[i];
    }

    // Update the target node of the last edge
    last_eid = edge_ids[n - 1];
    target_node = graph.edge(last_eid).target;
    QSMNode& last_node_ = graph.node(target_node);
    last_node_.x = x[n];
    last_node_.y = y[n];
    last_node_.z = z[n];
  }
}

void QSMbuilder::smooth_radii(int steps, double lambda, double mu)
{
  if (graph.edge_count() == 0) return;
  ServiceLocator::logger()("Smoothing radii");

  // Build axis map: axis_ID -> edge IDs sorted root->tip (descending subtree_length)
  std::map<int, std::vector<int>> axis_map;
  for (const auto& [eid, einfo] : graph.edges())
    axis_map[einfo.data.axis_ID].push_back(eid);

  for (auto& [axis, vec] : axis_map)
  {
    std::sort(vec.begin(), vec.end(), [this](int a, int b) {
      return graph.edge_data(a).subtree_length > graph.edge_data(b).subtree_length;
    });
  }

  for (const auto& [axis_id, edge_ids] : axis_map)
  {
    int n = edge_ids.size();
    if (n < 3) continue;  // Need at least 3 cylinders to smooth

    // Extract radii into a working buffer
    std::vector<double> r(n);
    for (int i = 0; i < n; ++i)
      r[i] = graph.edge_data(edge_ids[i]).radius;

    // Minimum radius floor: half the smallest observed radius on this axis,
    // so smoothing never collapses a cylinder to zero
    double r_min = *std::min_element(r.begin(), r.end()) * 0.5;
    r_min = std::max(r_min, 1e-4);  // Absolute floor (metres)

    // Taubin smoothing on radii
    // Convention: lambda > 0 (shrink pass), mu < 0 (inflate pass)
    // To bias toward inflation: |mu| > lambda
    // The caller controls this via the mu parameter (pass a negative value
    // whose magnitude exceeds lambda, e.g. lambda=0.5, mu=-0.6)
    for (int s = 0; s < steps; ++s)
    {
      // Pass 1: Shrink (lambda > 0)
      // Fix index 0 (root of axis): its radius is the best-known value
      // Let tip (index n-1) move freely so thin tips can self-correct
      for (int i = 1; i < n - 1; ++i)
      {
        double delta = 0.5 * (r[i-1] + r[i+1]) - r[i];
        r[i] += lambda * delta;
        r[i] = std::max(r[i], r_min);
      }

      // Pass 2: Inflate (mu < 0 ->negative delta moves radius up)
      for (int i = 1; i < n - 1; ++i)
      {
        double delta = 0.5 * (r[i-1] + r[i+1]) - r[i];
        r[i] += mu * delta;
        r[i] = std::max(r[i], r_min);
      }
    }

    // Write smoothed radii back
    for (int i = 0; i < n; ++i)
      graph.edge(edge_ids[i]).data.radius = r[i];
  }
}
}
