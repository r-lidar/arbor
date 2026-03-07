#include "QSM.h"
#include <algorithm>
#include <map>
#include <set>
#include <limits>

namespace arbor::qsm {

static inline std::array<double,3> normalize(const std::array<double,3>& v)
{
  double n = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  return {v[0]/n, v[1]/n, v[2]/n};
}

static inline std::array<double,3> cross(const std::array<double,3>& a, const std::array<double,3>& b)
{
  return { a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0] };
}

static inline double dot(const std::array<double,3>& a, const std::array<double,3>& b)
{
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

/**
 * Utility: Make an orthonormal frame given a direction.
 */
static inline void make_frame(const std::array<double,3>& dir, std::array<double,3>& ortho1, std::array<double,3>& ortho2)
{
  std::array<double,3> z_axis{0,0,1};
  // special case: dir is aligned with z-axis
  if (std::abs(dir[0]) < 1e-6 && std::abs(dir[1]) < 1e-6 && std::abs(dir[2]-1) < 1e-6)
  {
    ortho1 = {1,0,0};
  }
  else
  {
    ortho1 = {
      z_axis[1]*dir[2] - z_axis[2]*dir[1],
      z_axis[2]*dir[0] - z_axis[0]*dir[2],
      z_axis[0]*dir[1] - z_axis[1]*dir[0]
    };
    double n = std::sqrt(ortho1[0]*ortho1[0] + ortho1[1]*ortho1[1] + ortho1[2]*ortho1[2]);
    for(int i = 0; i < 3; ++i) ortho1[i] /= n;
  }
  ortho2 = {
    dir[1]*ortho1[2] - dir[2]*ortho1[1],
    dir[2]*ortho1[0] - dir[0]*ortho1[2],
    dir[0]*ortho1[1] - dir[1]*ortho1[0]
  };
  double n = std::sqrt(ortho2[0]*ortho2[0] + ortho2[1]*ortho2[1] + ortho2[2]*ortho2[2]);
  for(int i = 0; i < 3; ++i) ortho2[i] /= n;
}

/**
 * Rotate a frame from old_dir to new_dir using parallel transport.
 * This minimizes twist by rotating the frame around the axis perpendicular to both directions.
 */
static inline void propagate_frame(const std::array<double,3>& old_dir, const std::array<double,3>& new_dir, std::array<double,3>& ortho1, std::array<double,3>& ortho2)
{
  const double threshold = 1e-6;
  double d = dot(old_dir, new_dir);

  // Nearly parallel: no rotation needed
  if (d > 1.0 - threshold)
  {
    return;
  }

  // Nearly opposite: rotate 180° around ortho1
  if (d < -1.0 + threshold)
  {
    for(int i = 0; i < 3; ++i)
    {
      ortho2[i] = -ortho2[i];
    }
    return;
  }

  // General case: Rotate frame using Rodrigues' rotation formula
  std::array<double,3> axis = cross(old_dir, new_dir);
  double axis_len = std::sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]);
  for(int i=0; i<3; ++i) axis[i] /= axis_len;

  double angle = std::acos(std::max(-1.0, std::min(1.0, d)));
  double c = std::cos(angle);
  double s = std::sin(angle);
  double t = 1.0 - c;

  // Rotation matrix
  auto rotate = [&](const std::array<double,3>& v) -> std::array<double,3>
  {
    double dot_av = dot(axis, v);
    std::array<double,3> cross_av = cross(axis, v);
    return {
      t * axis[0] * dot_av + c * v[0] + s * cross_av[0],
                                                    t * axis[1] * dot_av + c * v[1] + s * cross_av[1],
                                                                                                  t * axis[2] * dot_av + c * v[2] + s * cross_av[2]
    };
  };

  ortho1 = rotate(ortho1);
  ortho2 = rotate(ortho2);

  // Re-orthogonalize to fix numerical drift
  double proj = dot(ortho1, new_dir);
  for(int i=0; i<3; ++i) ortho1[i] -= proj * new_dir[i];
  double n1 = std::sqrt(ortho1[0]*ortho1[0] + ortho1[1]*ortho1[1] + ortho1[2]*ortho1[2]);
  for(int i=0; i<3; ++i) ortho1[i] /= n1;

  ortho2 = cross(new_dir, ortho1);
}

/**
 * Build a continuous mesh along each axis by connecting rings of vertices at each node.
 */
void QSM::mesh(std::vector<std::array<double,3>>& vertices, std::vector<std::array<int,4>>& faces, std::vector<int>& cyl_ids, int resolution) const
{
  vertices.clear();
  faces.clear();
  cyl_ids.clear();

  // 1. Group edges by axis_ID
  std::map<int, std::vector<EdgeID>> axis_map;
  for (const auto& kv : edges())
    axis_map[kv.second.data.axis_ID].push_back(kv.first);

  // Global: for a given (axis_id, node_index) a ring index
  using AxisNode = std::pair<int, NodeID>;
  std::map<AxisNode, int> ring_offset; // value = starting vertex index in `vertices`

  for (const auto& axis_pair : axis_map)
  {
    int axis_id = axis_pair.first;
    const auto& edge_ids = axis_pair.second;

    // 2. For ordering: Build a node graph for this axis, so we can walk node->node
    std::set<NodeID> axis_nodes;
    std::map<NodeID, std::vector<std::pair<NodeID, EdgeID>>> out_edges, in_edges;
    for (EdgeID eid : edge_ids)
    {
      const auto& einfo = edges().at(eid);
      NodeID src = einfo.source, tgt = einfo.target;
      axis_nodes.insert(src);
      axis_nodes.insert(tgt);
      out_edges[src].emplace_back(tgt, eid);
      in_edges[tgt].emplace_back(src, eid);
    }

    // 3. Find the root node for this axis
    // The root is the node with NO incoming edges WITHIN this axis
    // Check each node in axis_nodes to see if it has incoming edges

    std::vector<NodeID> root_candidates;
    for (NodeID nid : axis_nodes)
    {
      // Check if this node has any incoming edges within this axis
      auto it = in_edges.find(nid);
      if (it == in_edges.end() || it->second.empty())
      {
        // No incoming edges within this axis - this is a root candidate
        root_candidates.push_back(nid);
      }
    }

    NodeID start;
    if (root_candidates.empty())
    {
      // Shouldn't happen in a well-formed tree, but fallback to smallest ID
      start = *std::min_element(axis_nodes.begin(), axis_nodes.end());
    }
    else if (root_candidates.size() == 1)
    {
      // Perfect - exactly one root
      start = root_candidates[0];
    }
    else
    {
      // Multiple roots (shouldn't happen in a connected axis)
      // Choose the one with smallest ID (most negative for axis 1)
      start = *std::min_element(root_candidates.begin(), root_candidates.end());
    }

    // 4. Walk along axis: build list of nodes (in axis order), and associate edge radii
    std::vector<NodeID> node_path;
    std::vector<EdgeID> edge_path;
    NodeID current_node = start;
    std::set<NodeID> visited;

    while(true)
    {
      node_path.push_back(current_node);
      visited.insert(current_node);

      auto it = out_edges.find(current_node);
      if (it == out_edges.end() || it->second.empty())
        break; // reached end of axis

      // Pick the first outgoing edge (in case of branching within axis)
      NodeID next = it->second[0].first;
      EdgeID next_eid = it->second[0].second;

      if (visited.count(next))
        break; // cycle detected - shouldn't happen

      edge_path.push_back(next_eid);
      current_node = next;
    }

    // 5. For each node in path, build a ring of vertices with propagated frames
    int N = node_path.size();
    if (N == 0) continue; // skip empty axes

    std::vector<int> ring_starts(N, -1);

    std::array<double,3> prev_dir{0,0,1};
    std::array<double,3> ortho1, ortho2;

    for (int idx = 0; idx < N; ++idx)
    {
      NodeID nid = node_path[idx];
      QSMNode qn = this->node(nid);

      // Local direction for this node
      std::array<double,3> dir{0,0,1};
      double radius = RADIUS_UNSET;
      int cyl_id = 0;

      if (idx==0 && edge_path.size() >= 1)
      {
        // First node: use direction to next node
        EdgeID eid = edge_path[0];
        QSMNode qchild = this->node(node_path[1]);
        dir = {
          qchild.x - qn.x,
          qchild.y - qn.y,
          qchild.z - qn.z
        };
        double dlen = std::sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
        for(int i=0;i<3;++i) dir[i] /= (dlen>1e-12? dlen : 1.0);
        radius = edges().at(eid).data.radius;
        cyl_id = edges().at(eid).data.cyl_ID;

        // Initialize frame for first node
        make_frame(dir, ortho1, ortho2);
      }
      else if (idx == N-1 && edge_path.size() >= 1)
      {
        // Last node: use direction from previous node
        EdgeID eid = edge_path[idx-1];
        QSMNode qpar = this->node(node_path[idx-1]);
        dir = {
          qn.x - qpar.x,
          qn.y - qpar.y,
          qn.z - qpar.z
        };
        double dlen = std::sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
        for(int i = 0; i < 3;++i) dir[i] /= (dlen>1e-12? dlen : 1.0);
        radius = edges().at(eid).data.radius;
        cyl_id = edges().at(eid).data.cyl_ID;

        // Propagate frame from previous node
        propagate_frame(prev_dir, dir, ortho1, ortho2);
      }
      else if (idx>0 && idx<N-1)
      {
        // Middle node: average direction
        EdgeID eid_prev = edge_path[idx-1];
        EdgeID eid_next = edge_path[idx];
        QSMNode qpar = this->node(node_path[idx-1]);
        QSMNode qchild = this->node(node_path[idx+1]);
        std::array<double,3> dir_in {
          qn.x - qpar.x,
          qn.y - qpar.y,
          qn.z - qpar.z
        };
        std::array<double,3> dir_out {
          qchild.x - qn.x,
          qchild.y - qn.y,
          qchild.z - qn.z
        };
        double len_in = std::sqrt(dir_in[0]*dir_in[0]+dir_in[1]*dir_in[1]+dir_in[2]*dir_in[2]);
        double len_out= std::sqrt(dir_out[0]*dir_out[0]+dir_out[1]*dir_out[1]+dir_out[2]*dir_out[2]);
        for(int i = 0; i < 3; ++i)
        {
          dir_in[i] /= (len_in>1e-12? len_in:1.0);
          dir_out[i]/= (len_out>1e-12?len_out:1.0);
        }
        for(int i = 0; i < 3; ++i)
        {
          dir[i] = 0.5*(dir_in[i]+dir_out[i]);
        }
        double dlen = std::sqrt(dir[0]*dir[0]+dir[1]*dir[1]+dir[2]*dir[2]);
        for(int i = 0; i < 3; ++i) dir[i]/=(dlen>1e-12?dlen:1.0);

        double r0 = edges().at(eid_prev).data.radius;
        double r1 = edges().at(eid_next).data.radius;
        radius = 0.5*(r0+r1);

        // Use the cyl_ID from the incoming edge (eid_prev) for this node's ring
        cyl_id = edges().at(eid_prev).data.cyl_ID;

        // Propagate frame from previous node
        propagate_frame(prev_dir, dir, ortho1, ortho2);
      }
      else
      {
        // Single node in axis (edge case)
        radius = 0.1;
        cyl_id = 0;
        make_frame(dir, ortho1, ortho2);
      }

      prev_dir = dir; // Store for next iteration

      int vstart = vertices.size();
      ring_starts[idx] = vstart;
      ring_offset[{axis_id, nid}] = vstart;

      // Create a ring of vertices
      for (int j = 0; j < resolution; ++j)
      {
        double theta = 2.0 * M_PI * j / resolution;
        double ct = std::cos(theta);
        double st = std::sin(theta);
        std::array<double,3> p {
          qn.x + radius * (ct*ortho1[0] + st*ortho2[0]),
          qn.y + radius * (ct*ortho1[1] + st*ortho2[1]),
          qn.z + radius * (ct*ortho1[2] + st*ortho2[2])
        };
        vertices.push_back(p);
        cyl_ids.push_back(cyl_id);  // Each vertex in the ring gets the same cyl_ID
      }
    }

    // Connect consecutive rings as quads
    for(int idx = 0; idx+1 < N; ++idx)
    {
      int a_start = ring_starts[idx];
      int b_start = ring_starts[idx+1];
      for(int j = 0; j < resolution; ++j)
      {
        int a0 = a_start + j;
        int a1 = a_start + (j+1)%resolution;
        int b0 = b_start + j;
        int b1 = b_start + (j+1)%resolution;
        faces.push_back({a0, a1, b1, b0});
      }
    }
  }
}

void QSM::qmesh(std::vector<std::array<double,3>>& vertices, std::vector<std::array<int,4>>& faces, std::vector<int>& node_ids, int resolution) const
{
  mesh(vertices, faces, node_ids, resolution);
}

void QSM::tmesh(std::vector<std::array<double,3>>& vertices, std::vector<std::array<int,3>>& faces, std::vector<int>& node_ids, int resolution) const
{
  std::vector<std::array<int,4>> quad_faces;
  mesh(vertices, quad_faces, node_ids, resolution);

  faces.clear();
  for (const auto& q : quad_faces) {
    faces.push_back({q[0], q[1], q[2]});
    faces.push_back({q[0], q[2], q[3]});
  }
}

}
