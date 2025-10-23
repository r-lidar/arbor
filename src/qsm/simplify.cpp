#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <unordered_map>
#include <unordered_set>
using namespace Rcpp;

// [[Rcpp::export]]
DataFrame qsm_simplify_cpp(DataFrame qsm, double max_length = 0.3)
{
  IntegerVector cyl_ID = qsm["cyl_ID"];
  IntegerVector parent_ID = qsm["parent_ID"];
  NumericVector startX = qsm["startX"];
  NumericVector startY = qsm["startY"];
  NumericVector startZ = qsm["startZ"];
  NumericVector endX = qsm["endX"];
  NumericVector endY = qsm["endY"];
  NumericVector endZ = qsm["endZ"];
  NumericVector radius = qsm["radius"];
  NumericVector cyl_length = qsm["cyl_length"];

  int n = cyl_ID.size();

  // Build parent-to-children map and cyl_ID -> row index
  std::unordered_map<int, std::vector<int>> children;
  std::unordered_map<int, int> id_to_row;
  std::unordered_map<int, int> child_count;

  for (int i = 0; i < n; ++i)
  {
    id_to_row[cyl_ID[i]] = i;
    children[parent_ID[i]].push_back(cyl_ID[i]);
    if (parent_ID[i] != 0)
      child_count[parent_ID[i]]++;
  }

  // Identify branching points and tips
  std::unordered_set<int> branching, tips;
  for (const auto& kv : child_count)
  {
    if (kv.second > 1) branching.insert(kv.first);
  }

  for (int i = 0; i < n; ++i)
  {
    if (children.find(cyl_ID[i]) == children.end())
      tips.insert(cyl_ID[i]);
  }

  std::unordered_set<int> important;
  for (int i = 0; i < n; ++i)
  {
    if (branching.count(cyl_ID[i]) || branching.count(parent_ID[i]) || tips.count(cyl_ID[i]))
      important.insert(cyl_ID[i]);
  }

  std::vector<bool> visited(n, false);
  std::vector<int> new_cyl_ID, new_parent_ID, new_original_row;
  std::vector<double> new_startX, new_startY, new_startZ, new_endX, new_endY, new_endZ, new_radius, new_length;

  std::unordered_map<int, int> old_to_new_id;
  int next_id = 1;

  for (int i = 0; i < n; ++i)
  {
    if (visited[i]) continue;

    int row = i;
    std::vector<int> chain = {row};
    visited[row] = true;

    // Grow forward
    int current = cyl_ID[row];
    while (children[current].size() == 1)
    {
      int child = children[current][0];
      int child_row = id_to_row[child];
      if (visited[child_row] || important.count(child)) break;
      chain.push_back(child_row);
      visited[child_row] = true;
      current = child;
    }

    int start = 0;
    while (start < (int)chain.size())
    {
      int idx = chain[start];
      int last = start;

      int s_idx = chain[start];
      double sx = startX[s_idx], sy = startY[s_idx], sz = startZ[s_idx];
      double ex = endX[s_idx], ey = endY[s_idx], ez = endZ[s_idx];
      double r = radius[s_idx];

      for (int j = start + 1; j < (int)chain.size(); ++j)
      {
        int k = chain[j];
        double dx = endX[k] - sx;
        double dy = endY[k] - sy;
        double dz = endZ[k] - sz;
        double d = std::sqrt(dx * dx + dy * dy + dz * dz);
        if (d > max_length) break;

        ex = endX[k];
        ey = endY[k];
        ez = endZ[k];
        r = (r + radius[k]) / 2.0;
        last = j;
      }

      double len = std::sqrt((ex - sx) * (ex - sx) + (ey - sy) * (ey - sy) + (ez - sz) * (ez - sz));

      new_cyl_ID.push_back(next_id);
      old_to_new_id[cyl_ID[chain[last]]] = next_id;

      new_parent_ID.push_back(0);  // temp, fix later
      new_startX.push_back(sx);
      new_startY.push_back(sy);
      new_startZ.push_back(sz);
      new_endX.push_back(ex);
      new_endY.push_back(ey);
      new_endZ.push_back(ez);
      new_radius.push_back(r);
      new_length.push_back(len);

      // Preserve original row index (1-based for R)
      new_original_row.push_back(chain[last] + 1);

      ++next_id;
      start = last + 1;
    }
  }

  // Fix parent_IDs
  for (size_t i = 0; i < new_cyl_ID.size(); ++i)
  {
    int old_idx = -1;
    for (int j = 0; j < n; ++j)
    {
      if (old_to_new_id[cyl_ID[j]] == new_cyl_ID[i])
        old_idx = j;
    }
    if (old_idx >= 0 && parent_ID[old_idx] != 0)
    {
      int old_pid = parent_ID[old_idx];
      auto it = old_to_new_id.find(old_pid);
      if (it != old_to_new_id.end())
        new_parent_ID[i] = it->second;
    }
  }

  return DataFrame::create(
    Named("startX") = new_startX,
    Named("startY") = new_startY,
    Named("startZ") = new_startZ,
    Named("endX") = new_endX,
    Named("endY") = new_endY,
    Named("endZ") = new_endZ,
    Named("radius") = new_radius,
    Named("cyl_ID") = new_cyl_ID,
    Named("parent_ID") = new_parent_ID,
    Named("original_row") = new_original_row  // <-- added column
  );
}
