/**
 * @file QSM.h
 * Project: Arbor
 *
 * Copyright (C) 2026 Jean-Romain Roussel (r-lidar) <info @ r-lidar.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef QSMGRAPH_H
#define QSMGRAPH_H

#include <cmath>

#include "DirectedGraph.h"

namespace arbor::qsm {

static constexpr float SUBTREE_LENGTH_UNSET    = -1.0;
static constexpr float SUBTREE_MAXZ_UNSET      = -1e30;
static constexpr float SUBTREE_VOLUME_UNSET    = -1;
static constexpr float RADIUS_UNSET            = -1.0;
static constexpr float DISTANCE_TO_ROOT_UNSET  = -1.0;
static constexpr float Z_EPS                   = 1e-9;

enum class EdgeQuality
{
  UNKNOWN    = 0,
    PROLONG    = 1,
    CONICALLOM = 2,
    POLYNOMIAL = 3,
    MEASURED   = 4,
    REFINED    = 5,
};

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
  // Core geometric / structural properties
  float radius           = RADIUS_UNSET;
  float subtree_length   = SUBTREE_LENGTH_UNSET;
  float distance_to_root = DISTANCE_TO_ROOT_UNSET;
  int   axis_ID          = 0;
  int   cyl_ID           = 0;
  int   parent_ID        = 0;
  int   branch_order     = 0;
  EdgeQuality quality    = EdgeQuality::UNKNOWN;


  // Temporary properties needed only while building the QSM
  float conic_allometry  = RADIUS_UNSET;
  float subtree_max_endZ = SUBTREE_MAXZ_UNSET;
  float subtree_volume   = SUBTREE_VOLUME_UNSET;

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
    float R = radius;
    if (R == RADIUS_UNSET)
      R = conic_allometry;
    if (R == RADIUS_UNSET)
      R = 0;

    return M_PI * R * R * length(src, tgt);
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

  inline double distance(const QSMNode& src, const QSMNode& tgt, double px, double py, double pz) const
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
  void tmesh(std::vector<std::array<double,3>>& vertices, std::vector<std::array<int,3>>& faces, std::vector<int>& cyl_ids, int sides = 16) const;
  void qmesh(std::vector<std::array<double,3>>& vertices, std::vector<std::array<int,4>>& faces, std::vector<int>& cyl_ids, int sides = 16) const;
  void read(const std::string& filename);
  void write(const std::string& filename, bool binary = true) const;
  void validate() const;
  NodeID find_root_node() const;
  QSM stem() const;
  QSM merchantable(double merch_radius) const;
  double dbh(double d, double* xyz = nullptr, double* n = nullptr) const;
  //void dump(std::ostream& os = std::cout, bool detailed = false) const;

public:
  std::vector<std::string> messages;

private:
  void mesh(std::vector<std::array<double,3>>& vertices, std::vector<std::array<int,4>>& faces, std::vector<int>& node_ids, int resolution) const;
  void read_qsm (const std::string& filename);
  void write_ply(const std::string& filename, bool binary) const;
  void write_stl(const std::string& filename, bool binary) const;
  void write_obj(const std::string& filename) const;
  void write_csv(const std::string& filename) const;
  void write_qsm(const std::string& filename) const;
};

} // namespace arbor::qsm

#endif // QSMGRAPH_H
