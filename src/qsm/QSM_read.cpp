/**
 * @file QSM_read.cpp
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

#include "QSM.h"

#include <fstream>
#include <string>
#include <cstdint>
#include <stdexcept>
#include <algorithm>

namespace arbor::qsm {

void QSM::read(const std::string& filename)
{
  auto pos = filename.find_last_of('.');
  if (pos == std::string::npos)
    throw std::runtime_error("QSM::read: filename has no extension: " + filename);

  std::string ext = filename.substr(pos + 1);
  std::transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char c) { return std::tolower(c); });

  if (ext == "qsm")
    read_qsm(filename);
  else
    throw std::runtime_error("Unsupported extension '." + ext + "'. Supported read formats: .qsm");
}

void QSM::read_qsm(const std::string& filename)
{
  static constexpr uint32_t MIN_NODE_SIZE = 16;
  static constexpr uint32_t MIN_EDGE_SIZE = 30;

  std::ifstream in(filename, std::ios::binary);
  if (!in.is_open())
    throw std::runtime_error("Cannot open file: " + filename);

  // ---- Parse ASCII header --------------------------------------------------
  uint64_t n_nodes   = 0;
  uint64_t n_edges   = 0;
  uint32_t node_size = 0;
  uint32_t edge_size = 0;
  double   xoffset   = 0.0;
  double   yoffset   = 0.0;
  double   zoffset   = 0.0;
  bool     data_seen = false;

  messages.clear();

  std::string line;
  while (std::getline(in, line))
  {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.empty() || line[0] == '#') continue;

    const auto        space = line.find(' ');
    const std::string key   = (space == std::string::npos) ? line              : line.substr(0, space);
    const std::string val   = (space == std::string::npos) ? std::string{}     : line.substr(space + 1);

    if (key == "DATA")
    {
      if (val != "binary")
        throw std::runtime_error("Only 'DATA binary' is supported, got: DATA " + val);
      data_seen = true;
      break;
    }
    else if (key == "VERSION")
    {
      if (val.rfind("1.", 0) != 0)
        throw std::runtime_error("File VERSION is '" + val + "'; reader supports 1.x");
    }
    else if (key == "NODES")     { n_nodes   = std::stoull(val); }
    else if (key == "EDGES")     { n_edges   = std::stoull(val); }
    else if (key == "NODE_SIZE") { node_size = static_cast<uint32_t>(std::stoul(val)); }
    else if (key == "EDGE_SIZE") { edge_size = static_cast<uint32_t>(std::stoul(val)); }
    else if (key == "X_OFFSET")  { xoffset   = std::stod(val); }
    else if (key == "Y_OFFSET")  { yoffset   = std::stod(val); }
    else if (key == "Z_OFFSET")  { zoffset   = std::stod(val); }
    else if (key == "MESSAGE")
    {
      std::string msg;
      msg.reserve(val.size());
      for (std::size_t i = 0; i < val.size(); ++i)
      {
        if (val[i] == '\\' && i + 1 < val.size() && val[i + 1] == 'n')
        {
          msg += '\n';
          i++;
        }
        else
        {
          msg += val[i];
        }
      }
      messages.push_back(std::move(msg));
    }
    // Unknown keys silently ignored (forward-compatibility, §2.1 of spec)
  }

  if (!data_seen)
    throw std::runtime_error("'DATA binary' sentinel not found in header.");

  // ---- Validate record sizes -----------------------------------------------
  if (node_size == 0)
    throw std::runtime_error("NODE_SIZE not specified in header.");
  if (edge_size == 0)
    throw std::runtime_error("EDGE_SIZE not specified in header.");
  if (node_size < MIN_NODE_SIZE)
    throw std::runtime_error("NODE_SIZE " + std::to_string(node_size) + " is below the minimum of " + std::to_string(MIN_NODE_SIZE) + " bytes.");
  if (edge_size < MIN_EDGE_SIZE)
    throw std::runtime_error("EDGE_SIZE " + std::to_string(edge_size) + " is below the minimum of " + std::to_string(MIN_EDGE_SIZE) + " bytes.");

  const uint32_t node_extra = node_size - MIN_NODE_SIZE; // bytes to skip per node
  const uint32_t edge_extra = edge_size - MIN_EDGE_SIZE; // bytes to skip per edge

  // ---- Reset graph ---------------------------------------------------------
  clear();

  // ---- Read node records ---------------------------------------------------
  // Mandatory layout (16 bytes): int32 node_id | float x | float y | float z
  // Any extra bytes (node_extra) are implementation-defined and must be skipped.
  for (uint64_t i = 0; i < n_nodes; ++i)
  {
    int32_t nid = 0;
    float   fx  = 0.0f, fy = 0.0f, fz = 0.0f;

    in.read(reinterpret_cast<char*>(&nid), sizeof(int32_t));
    in.read(reinterpret_cast<char*>(&fx),  sizeof(float));
    in.read(reinterpret_cast<char*>(&fy),  sizeof(float));
    in.read(reinterpret_cast<char*>(&fz),  sizeof(float));
    if (node_extra > 0) in.ignore(node_extra);

    if (!in)
      throw std::runtime_error("Unexpected EOF reading node record " + std::to_string(i) + " of " + std::to_string(n_nodes) + ".");

    // Reconstruct full double-precision geographic coordinates
    QSMNode node = {
      static_cast<double>(fx) + xoffset,
      static_cast<double>(fy) + yoffset,
      static_cast<double>(fz) + zoffset
    };
    insert_node(static_cast<NodeID>(nid), node);
  }

  // ---- Read edge records ---------------------------------------------------
  // Mandatory layout (30 bytes): see write_qbf above.
  // Any extra bytes (edge_extra) are implementation-defined and must be skipped.
  for (uint64_t i = 0; i < n_edges; ++i)
  {
    int32_t eid    = 0, source = 0, target = 0;
    uint8_t bo_u8  = 0, q_u8  = 0;
    QSMEdge e{};

    in.read(reinterpret_cast<char*>(&eid),             sizeof(int32_t));
    in.read(reinterpret_cast<char*>(&source),          sizeof(int32_t));
    in.read(reinterpret_cast<char*>(&target),          sizeof(int32_t));
    in.read(reinterpret_cast<char*>(&e.radius),           sizeof(float));
    in.read(reinterpret_cast<char*>(&e.subtree_length),   sizeof(float));
    in.read(reinterpret_cast<char*>(&e.distance_to_root), sizeof(float));
    in.read(reinterpret_cast<char*>(&e.axis_ID),          sizeof(int32_t));
    in.read(reinterpret_cast<char*>(&bo_u8),              sizeof(uint8_t));
    in.read(reinterpret_cast<char*>(&q_u8),               sizeof(uint8_t));
    if (edge_extra > 0) in.ignore(edge_extra);

    if (!in)
      throw std::runtime_error("Unexpected EOF reading edge record " +  std::to_string(i) + " of " + std::to_string(n_edges) + ".");

    e.branch_order     = static_cast<int>(bo_u8);
    e.quality          = static_cast<EdgeQuality>(q_u8);

    // Runtime-only scratch fields reset to their sentinel values
    e.conic_allometry  = RADIUS_UNSET;
    e.subtree_max_endZ = SUBTREE_MAXZ_UNSET;
    e.subtree_volume   = SUBTREE_VOLUME_UNSET;

    insert_edge(static_cast<EdgeID>(eid), static_cast<NodeID>(source), static_cast<NodeID>(target), e);
  }

  // ---- Post-read integrity check -------------------------------------------
  if (node_count() != n_nodes)
    throw std::runtime_error("Node count mismatch: header declared " +  std::to_string(n_nodes) + ", graph contains " + std::to_string(node_count()) + ".");
  if (edge_count() != n_edges)
    throw std::runtime_error("Edge count mismatch: header declared " +  std::to_string(n_edges) + ", graph contains " + std::to_string(edge_count()) + ".");

  // ---- Reconstruct legacy flat-format fields from graph topology ----------
  // cyl_ID   = edge_id (the graph key, not stored separately)
  // parent_ID = edge_id of the single incoming edge of source, or 0 at root
  for (auto& [eid, einfo] : edges())
  {
    QSMEdge& e         = einfo.data;
    e.cyl_ID           = static_cast<int>(eid);
    const auto& inc    = incoming_edges(einfo.source);
    e.parent_ID        = inc.empty() ? 0 : static_cast<int>(inc[0]);
  }
}

} // namespace arbor::qsm
