#include <Rcpp.h>
#include <unordered_map>
#include <vector>
#include <algorithm>

using namespace Rcpp;
using namespace std;

// Recursive function to compute subtree lengths
double compute_subtree_length(
    int node_id,
    const unordered_map<int, vector<int>> &children_map,
    const unordered_map<int, double> &length_map,
    unordered_map<int, double> &subtree_lengths
)
{
  if (subtree_lengths.find(node_id) != subtree_lengths.end())
    return subtree_lengths[node_id];

  auto it = children_map.find(node_id);
  if (it == children_map.end())
  {
    subtree_lengths[node_id] = 0.0;
    return 0.0;
  }

  double max_length = 0.0;
  for (int child_id : it->second)
  {
    double child_length = compute_subtree_length(child_id, children_map, length_map, subtree_lengths);
    max_length = max(max_length, child_length + length_map.at(child_id));
  }

  subtree_lengths[node_id] = max_length;
  return max_length;
}

// Recursive function to assign subtree IDs and branching orders
void assign_subtree_ids(
    int node_id,
    int current_subtree_id,
    int current_branching_order,
    const unordered_map<int, vector<int>> &children_map,
    const unordered_map<int, double> &length_map,
    const unordered_map<int, double> &subtree_lengths,
    unordered_map<int, int> &subtree_ids,
    unordered_map<int, int> &branching_orders,
    int &next_subtree_id
)
{
  subtree_ids[node_id] = current_subtree_id;
  branching_orders[node_id] = current_branching_order;

  auto it = children_map.find(node_id);
  if (it == children_map.end())
    return;

  // Identify the main child (on the longest path)
  int main_child = -1;
  double max_path = -1;

  for (int child_id : it->second)
  {
    double path_length = subtree_lengths.at(child_id) + length_map.at(child_id);
    if (path_length > max_path)
    {
      max_path = path_length;
      main_child = child_id;
    }
  }

  for (int child_id : it->second)
  {
    if (child_id == main_child)
    {
      assign_subtree_ids(
        child_id, current_subtree_id, current_branching_order,
        children_map, length_map, subtree_lengths,
        subtree_ids, branching_orders, next_subtree_id
      );
    }
    else
    {
      int new_subtree_id = next_subtree_id++;
      assign_subtree_ids(
        child_id, new_subtree_id, current_branching_order + 1,
        children_map, length_map, subtree_lengths,
        subtree_ids, branching_orders, next_subtree_id
      );
    }
  }
}

DataFrame cpp_compute_architecture(DataFrame qsm, int root_id = 1)
{
  IntegerVector cyl_ID = qsm["cyl_ID"];
  IntegerVector parent_ID = qsm["parent_ID"];
  NumericVector length = qsm["length"];
  int n = cyl_ID.size();

  unordered_map<int, vector<int>> children_map;
  unordered_map<int, double> length_map;

  for (int i = 0; i < n; ++i)
  {
    int cid = cyl_ID[i];
    int pid = parent_ID[i];
    children_map[pid].push_back(cid);
    length_map[cid] = length[i];
  }

  unordered_map<int, double> subtree_lengths;
  compute_subtree_length(root_id, children_map, length_map, subtree_lengths);

  unordered_map<int, int> subtree_ids;
  unordered_map<int, int> branching_orders;
  int next_subtree_id = 2;
  assign_subtree_ids(
    root_id, 1, 1,
    children_map, length_map, subtree_lengths,
    subtree_ids, branching_orders, next_subtree_id
  );

  NumericVector lengths(n);
  IntegerVector ids(n);
  IntegerVector orders(n);
  for (int i = 0; i < n; ++i)
  {
    int cid = cyl_ID[i];
    lengths[i] = subtree_lengths.count(cid) ? subtree_lengths.at(cid) : 0.0;
    ids[i] = subtree_ids.count(cid) ? subtree_ids.at(cid) : NA_INTEGER;
    orders[i] = branching_orders.count(cid) ? branching_orders.at(cid) : NA_INTEGER;
  }

  return DataFrame::create(
    Named("cyl_ID") = cyl_ID,
    Named("subtree_length") = lengths,
    Named("axis_ID") = ids,
    Named("branching_order") = orders
  );
}
