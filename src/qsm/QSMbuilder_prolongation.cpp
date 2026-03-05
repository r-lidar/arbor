#include "QSMbuilder.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace arbor::qsm {

void QSMbuilder::prolongate(double d, double L)
{
  if (d <= 0.0) return;

  logger("Prolongation to the ground");

  // Collect main axis edges (axis_ID == 1), ordered root → tip (descending subtree_length)
  std::vector<int> axis_eids;
  for (const auto& [eid, einfo] : graph.edges())
    if (einfo.data.axis_ID == 1) axis_eids.push_back(eid);

    const size_t n = axis_eids.size();
    if (n < 2) return;

    std::sort(axis_eids.begin(), axis_eids.end(), [this](int a, int b) {
      return graph.edge_data(a).subtree_length > graph.edge_data(b).subtree_length;
    });

    // Find root edge: axis_ID==1 edge with the maximum subtree_length
    // (the root of the main axis always has the highest subtree_length)
    int root_eid = axis_eids[0];  // already sorted descending by subtree_length

    // Compute cylinder lengths and cumulative lengths
    std::vector<double> lens(n);
    std::vector<double> cum(n);
    double total = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
      const auto& einfo = graph.edge(axis_eids[i]);
      lens[i] = graph.edge_data(axis_eids[i]).length(
        graph.node(einfo.source), graph.node(einfo.target));
      total += lens[i];
      cum[i] = total;
    }
    if (total == 0.0) return;

    // Take only the first ~10% of the axis or 30 cm
    double cutoff = 0.1 * total;
    size_t k = 0;
    while (k < n && cum[k] <= cutoff) k++;
    if (k == 0) k = 1;
    if (k < n)
    {
      double s = 0.0;
      for (size_t i = 0; i < k; i++) s += lens[i];
      if (s < 0.3) k++;
    }
    if (k > n) k = n;

    // Root node position
    QSM::NodeID root_src_nid = graph.edge(root_eid).source;
    const QSMNode& root_node    = graph.node(root_src_nid);
    const QSMNode& last_node    = graph.node(graph.edge(axis_eids[k - 1]).target);

    // Estimated orientation of the trunk
    double dx = last_node.x - root_node.x;
    double dy = last_node.y - root_node.y;
    double dz = last_node.z - root_node.z;
    double N  = std::sqrt(dx*dx + dy*dy + dz*dz);
    if (N <= 0.0) return;
    double ox = dx / N;
    double oy = dy / N;
    double oz = dz / N;

    // Adjust prolongation distance by angle
    double d_adj;
    if (std::abs(oz) < 1e-9)
      d_adj = d;
    else
      d_adj = d / oz;

    // Create subdivided segments going from root.start downward
    int nseg = std::max(1, int(std::ceil(d_adj / L)));

    double root_subtree = graph.edge_data(root_eid).subtree_length;
    if (root_subtree == SUBTREE_LENGTH_UNSET)
      throw std::runtime_error("Invalid QSM, the root has no subtree length");

    // Start with the most negative ID for the deepest (new root) segment
    int next_id  = -nseg;
    int prev_cyl_id = 0;  // The deepest segment has no parent (it's the new root)

    // We'll build from the deepest point upward, then connect to the former root
    QSM::NodeID prev_node_id = -1;  // Will be set in the loop
    QSM::NodeID first_prolongation_node = -1;  // The deepest node (new root)

    // Create segments from deepest (new root) to shallowest (connects to former root)
    for (int i = nseg; i >= 1; i--)
    {
      double f1 = double(i) / nseg;
      double f2 = double(i - 1) / nseg;

      double x1 = root_node.x - ox * d_adj * f1;
      double y1 = root_node.y - oy * d_adj * f1;
      double z1 = root_node.z - oz * d_adj * f1;

      double x2 = root_node.x - ox * d_adj * f2;
      double y2 = root_node.y - oy * d_adj * f2;
      double z2 = root_node.z - oz * d_adj * f2;

      QSM::NodeID node_id_1, node_id_2;

      if (i == nseg)
      {
        // First iteration: create the deepest node (new root)
        node_id_1 = graph.add_node({x1, y1, z1});
        first_prolongation_node = node_id_1;
      }
      else
      {
        // Reuse the previous target node as source
        node_id_1 = prev_node_id;
      }

      if (i == 1)
      {
        // Last iteration: connect to the original root node
        node_id_2 = root_src_nid;
      }
      else
      {
        // Create intermediate node
        node_id_2 = graph.add_node({x2, y2, z2});
      }

      QSMEdge ed;
      ed.cyl_ID        = next_id;
      ed.parent_ID     = prev_cyl_id;
      ed.axis_ID       = 1;
      ed.branch_order  = 1;
      ed.subtree_length = root_subtree + d_adj - (d_adj / nseg) * (i - 1);

      graph.add_edge(node_id_1, node_id_2, ed);

      prev_cyl_id   = next_id;
      prev_node_id  = node_id_2;
      next_id++;
    }

    // Update the former root edge to connect it to the prolongation chain
    // The former root edge should now have parent_ID pointing to the last prolongation segment
    graph.edge_data(root_eid).parent_ID = prev_cyl_id;  // prev_cyl_id is now -1 (the last segment before root)
}

void QSMbuilder::estimate_prolongation(const PointCloud& tree)
{
  prolongation_distance = 0.0;

  if (!tree.has_hag() || graph.edge_count() == 0) return;

  // Compute the minimum startZ across all edges (source node Z) –
  // matches the original behaviour that used min(startZ) over all cylinders.
  double min_start_z = std::numeric_limits<double>::max();
  for (const auto& [eid, einfo] : graph.edges())
  {
    double sz = graph.node(einfo.source).z;
    if (sz < min_start_z) min_start_z = sz;
  }

  // Find max(hag) for points where Z <= min_start_z
  double max_hag = -std::numeric_limits<double>::max();
  bool found = false;

  for (size_t i = 0; i < tree.size(); ++i)
  {
    if (tree.get_z(i) <= min_start_z)
    {
      double hag = tree.get_hag(i);
      if (hag > max_hag) { max_hag = hag; found = true; }
    }
  }

  if (found) prolongation_distance = max_hag;
}

}
