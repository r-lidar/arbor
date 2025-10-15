#include <Rcpp.h>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <vector>
#include <algorithm>

using namespace Rcpp;

DataFrame cpp_compute_topology(DataFrame skel)
{
  IntegerVector seg_ID = skel["seg_ID"];
  IntegerVector bearer_ID = skel["bearer_ID"];
  NumericVector length = skel["length"];

  int n = seg_ID.size();
  std::vector<double> bear_length(n, 0.0);
  std::vector<int> axis_ID(n, 0);
  std::vector<int> section(n, 0);

  // Step 1: Build reverse bearer map: parent -> children
  std::unordered_map<int, std::vector<int>> children_map;
  for (int i = 0; i < n; ++i)
  {
    if (bearer_ID[i] > 0)
      children_map[bearer_ID[i]].push_back(seg_ID[i]);
  }

  // Step 2: Build seg_ID to index map
  std::unordered_map<int, int> id_to_index;
  for (int i = 0; i < n; ++i) id_to_index[seg_ID[i]] = i;

  // Step 3: Compute bear_length in reverse order
  for (int i = n - 1; i >= 0; --i)
  {
    int s = seg_ID[i];
    double sum_child_length = 0.0;
    auto it = children_map.find(s);
    if (it != children_map.end())
    {
      for (int cid : it->second)
        sum_child_length += bear_length[id_to_index[cid]];
    }
    bear_length[i] = length[i] + sum_child_length;
  }

  // Step 4: Assign axis_ID and section
  std::queue<int> q;
  int cur_ID = 1;
  int cur_sec = 1;

  // Start with root segments (bearer_ID == 0)
  for (int i = 0; i < n; ++i)
  {
    if (bearer_ID[i] == 0) {
      axis_ID[i] = cur_ID;
      section[i] = cur_sec;
    }
  }

  // Initialize traversal from first root segment
  int cur_seg = -1;
  for (int i = 0; i < n; ++i)
  {
    if (bearer_ID[i] == 0) {
      cur_seg = seg_ID[i];
      break;
    }
  }

  while (std::any_of(axis_ID.begin(), axis_ID.end(), [](int id) { return id == 0; }))
  {
    int cur_idx = id_to_index[cur_seg];
    axis_ID[cur_idx] = cur_ID;
    section[cur_idx] = cur_sec;

    std::vector<int>& childs = children_map[cur_seg];

    if (childs.size() == 0)
    {
      if (!q.empty())
      {
        cur_seg = q.front();
        q.pop();
        ++cur_ID;
        ++cur_sec;
        continue;
      }
      else
      {
        break;
      }
    }
    else if (childs.size() == 1)
    {
      cur_seg = childs[0];
    }
    else
    {
      // Find max bear_length child
      int max_idx = 0;
      double max_len = bear_length[id_to_index[childs[0]]];
      for (size_t j = 1; j < childs.size(); ++j)
      {
        double len = bear_length[id_to_index[childs[j]]];
        if (len > max_len)
        {
          max_idx = j;
          max_len = len;
        }
      }

      cur_seg = childs[max_idx];
      for (size_t j = 0; j < childs.size(); ++j)
      {
        if ((int)j != max_idx)
          q.push(childs[j]);
      }

      ++cur_sec; // increase section only on branch
    }
  }

  return DataFrame::create(
    Named("seg_ID") = seg_ID,
    Named("bearer_ID") = bearer_ID,
    Named("length") = length,
    Named("bear_length") = bear_length,
    Named("axis_ID") = axis_ID,
    Named("section") = section
  );
}
