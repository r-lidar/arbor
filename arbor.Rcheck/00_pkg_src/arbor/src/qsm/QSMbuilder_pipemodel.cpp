#include "QSMbuilder.h"
#include <vector>
#include <cmath>
#include <algorithm>

namespace {

// Struct to hold the state that we previously passed via function arguments
struct WorkItem
{
  NodeID node;
  float r_anchor;
  float s_anchor;
};

} // anonymous namespace

namespace arbor::qsm {

// ─────────────────────────────────────────────────────────────────────────────
void QSMbuilder::pipe_model_reconstruction(double tip_radius)
{
  float t_rad = static_cast<float>(tip_radius);
  NodeID root = graph.find_root_node();

  // Pre-allocate an explicit stack to prevent reallocations.
  // 1024 is an arbitrary safe buffer; it will grow dynamically if needed.
  std::vector<WorkItem> stack;
  stack.reserve(1024);
  stack.push_back({root, RADIUS_UNSET, 0.0f});

  while (!stack.empty())
  {
    // Pop the current state off the stack
    WorkItem current = stack.back();
    stack.pop_back();

    const auto& out_eids = graph.outgoing_edges(current.node);
    if (out_eids.empty()) continue; // leaf node

    // ── Classify outgoing edges without heap allocations ─────────────────
    int unmeasured_count = 0;
    int measured_count = 0;
    float total_s2_all = 0.0f;

    for (EdgeID eid : out_eids)
    {
      const QSMEdge& e = graph.edge_data(eid);
      if (e.radius != RADIUS_UNSET) {
        measured_count++;
        // Measured edges carry their own anchor; push to stack unconditionally
        stack.push_back({graph.edge(eid).target, e.radius, e.subtree_length});
      } else {
        unmeasured_count++;
      }

      // Rule B requires sum of squared subtree lengths for ALL edges
      total_s2_all += (e.subtree_length * e.subtree_length);
    }

    if (unmeasured_count == 0) continue;

    // ── Guard: no valid anchor yet ────────────────────────────────────────
    if (current.r_anchor == RADIUS_UNSET || current.s_anchor <= 0.0f)
    {
      for (EdgeID eid : out_eids)
      {
        QSMEdge& e = graph.edge_data(eid);
        if (e.radius == RADIUS_UNSET)
        {
          e.radius  = t_rad;
          e.quality = CONICALLOM;
          stack.push_back({graph.edge(eid).target, t_rad, e.subtree_length});
        }
      }
      continue;
    }

    // ── Choose rule ───────────────────────────────────────────────────────
    const bool is_fork = (unmeasured_count > 1) || (measured_count > 0);

    if (!is_fork)
    {
      // ── Rule A: straight segment — linear taper ───────────────────────
      for (EdgeID eid : out_eids)
      {
        QSMEdge& e = graph.edge_data(eid);
        if (e.radius == RADIUS_UNSET) // There will only be 1
        {
          float r = t_rad + (current.r_anchor - t_rad) * (e.subtree_length / current.s_anchor);
          r = std::max(r, t_rad);

          e.radius  = r;
          e.quality = CONICALLOM;

          stack.push_back({graph.edge(eid).target, r, e.subtree_length});
        }
      }
    }
    else
    {
      // ── Rule B: fork — direct pipe model ──────────────────────────────
      for (EdgeID eid : out_eids)
      {
        QSMEdge& e = graph.edge_data(eid);
        if (e.radius == RADIUS_UNSET)
        {
          float r = (total_s2_all > 0.0f)
          ? current.r_anchor * (e.subtree_length / std::sqrt(total_s2_all))
            : t_rad;

          r = std::max(r, t_rad);

          e.radius  = r;
          e.quality = CONICALLOM;

          stack.push_back({graph.edge(eid).target, r, e.subtree_length});
        }
      }
    }
  }
}

} // namespace arbor::qsm
