#ifndef QSM_H
#define QSM_H

#include <unordered_map>
#include <vector>
#include <limits>

using NodeId = int;

class QSM
{
public:
  QSM() = default;

  // build internal maps from vectors
  void build_from_vectors(const std::vector<int>& cyl_ID,
                          const std::vector<int>& parent_ID,
                          const std::vector<double>& length,
                          const std::vector<double>& startZ,
                          const std::vector<double>& endZ);

  // high-level compute entry (fills subtree_lengths, subtree_max_z, subtree_ids, branching_orders)
  void compute_architecture(NodeId root_id = 1);

  // accessors for results (node id -> value)
  const std::unordered_map<NodeId,double>& get_subtree_lengths() const { return subtree_lengths_; }
  const std::unordered_map<NodeId,double>& get_subtree_max_z() const { return subtree_max_z_; }
  const std::unordered_map<NodeId,int>& get_subtree_ids() const { return subtree_ids_; }
  const std::unordered_map<NodeId,int>& get_branching_orders() const { return branching_orders_; }

private:
  // raw input storage (maps node id -> value)
  std::unordered_map<NodeId, std::vector<NodeId>> children_map_;
  std::unordered_map<NodeId, double> length_map_;
  std::unordered_map<NodeId, double> startz_map_;
  std::unordered_map<NodeId, double> endz_map_;

  // results
  std::unordered_map<NodeId, double> subtree_lengths_;
  std::unordered_map<NodeId, double> subtree_max_z_;
  std::unordered_map<NodeId, int> subtree_ids_;
  std::unordered_map<NodeId, int> branching_orders_;

  // recursive helpers (pure STL)
  double compute_subtree_length(NodeId node_id);
  double compute_subtree_max_z(NodeId node_id);
  void assign_subtree_ids(NodeId node_id, int current_subtree_id, int current_branching_order, int &next_subtree_id);
};

#endif // QSM_H
