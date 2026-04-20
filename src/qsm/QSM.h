#ifndef QSMGRAPH_H
#define QSMGRAPH_H

#include <cmath>

#include "DirectedGraph.h"

namespace arbor::qsm {

static constexpr double SUBTREE_LENGTH_UNSET    = -1.0;
static constexpr double SUBTREE_MAXZ_UNSET      = -1e300;
static constexpr double SUBTREE_VOLUME_UNSET    = -1;
static constexpr double RADIUS_UNSET            = -1.0;
static constexpr double Z_EPS                   = 1e-9;

// A node in the QSM graph: a 3-D junction point in the tree structure.
struct QSMNode
{
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
};

// An edge in the QSM graph: a cylinder connecting two junction points.
// Legacy fields (cyl_ID, parent_ID) are retained for backward-compatible
// conversion to and from the flat QSM / QSMcylinder representation.
struct QSMEdge
{
  // Backward-compatibility IDs populated during conversion or build.
  int cyl_ID    = 0;
  int parent_ID = 0;

  // Geometric / structural properties
  double radius           = RADIUS_UNSET;
  double conic_allometry  = RADIUS_UNSET;
  double subtree_length   = SUBTREE_LENGTH_UNSET;
  double subtree_max_endZ = SUBTREE_MAXZ_UNSET;
  double subtree_volume   = SUBTREE_VOLUME_UNSET;
  int    axis_ID          = 0;
  int    branch_order     = 0;

  // ---- Computed geometry (requires source/target node positions) ----

  inline double length(const QSMNode& src, const QSMNode& tgt) const noexcept
  {
    double dx = tgt.x - src.x;
    double dy = tgt.y - src.y;
    double dz = tgt.z - src.z;
    return std::sqrt(dx*dx + dy*dy + dz*dz);
  }

  inline double volume(const QSMNode& src, const QSMNode& tgt) const noexcept
  {
    return M_PI * radius * radius * length(src, tgt);
  }

  inline double angle(const QSMNode& src, const QSMNode& tgt) const noexcept
  {
    double dx = tgt.x - src.x;
    double dy = tgt.y - src.y;
    double dz = tgt.z - src.z;
    double len = std::sqrt(dx*dx + dy*dy + dz*dz);
    if (len < 1e-6) return 0.0;
    return std::acos(dz / len) * 180.0 / M_PI;
  }

  inline double distance(const QSMNode& src, const QSMNode& tgt,
                         double px, double py, double pz) const
  {
    double ab_x = tgt.x - src.x;
    double ab_y = tgt.y - src.y;
    double ab_z = tgt.z - src.z;
    double ab_sq = ab_x*ab_x + ab_y*ab_y + ab_z*ab_z;
    if (ab_sq < 1e-12) ab_sq = 1e-12;

    double ap_x = px - src.x;
    double ap_y = py - src.y;
    double ap_z = pz - src.z;
    double t = (ap_x*ab_x + ap_y*ab_y + ap_z*ab_z) / ab_sq;
    if (t < 0.0) t = 0.0;
    else if (t > 1.0) t = 1.0;

    double c_x = src.x + t * ab_x;
    double c_y = src.y + t * ab_y;
    double c_z = src.z + t * ab_z;

    double dx = px - c_x;
    double dy = py - c_y;
    double dz = pz - c_z;
    double dist_axis = std::sqrt(dx*dx + dy*dy + dz*dz);
    return (dist_axis > radius) ? (dist_axis - radius) : 0.0;
  }
};

class QSM : public DirectedGraph<QSMNode, QSMEdge>
{
public:
  enum class MeshMode { Cylinders, Continuous };
  void tmesh(std::vector<std::array<double,3>>& vertices, std::vector<std::array<int,3>>& faces, std::vector<int>& cyl_ids, int sides = 16) const;
  void qmesh(std::vector<std::array<double,3>>& vertices, std::vector<std::array<int,4>>& faces, std::vector<int>& cyl_ids, int sides = 16) const;
  void write(const std::string& filename, bool binary = true) const;
  std::vector<std::string> messages;
  //void dump(std::ostream& os = std::cout, bool detailed = false) const;
private:
  void mesh(std::vector<std::array<double,3>>& vertices, std::vector<std::array<int,4>>& faces, std::vector<int>& node_ids, int resolution) const;
  void write_ply(const std::string& filename, bool binary) const;
  void write_stl(const std::string& filename, bool binary) const;
  void write_obj(const std::string& filename) const;
  void write_csv(const std::string& filename) const;
};

} // namespace arbor::qsm

#endif // QSMGRAPH_H
