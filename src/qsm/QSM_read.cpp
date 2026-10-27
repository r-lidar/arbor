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
#include "libqsm.h"

#include <fstream>
#include <string>
#include <cstdint>
#include <stdexcept>
#include <algorithm>
#include <vector>

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
  libqsm::QSMreader reader(filename);

  for (int i = 0 ; i < reader.get_message_count() ; i++)
    messages.push_back(reader.get_message(i));

  // Edges refer to nodes by their 1-based position in the file. Map them to
  // the NodeIDs assigned by add_node() rather than assuming they coincide.
  std::vector<NodeID> node_ids;
  node_ids.reserve(reader.nodes().size());
  for (const auto& node : reader.nodes())
    node_ids.push_back(add_node(node));

  for (const auto& edge : reader.edges())
  {
    if (edge.source < 1 || edge.source > node_ids.size() || edge.target < 1 || edge.target > node_ids.size())
      throw std::runtime_error("Invalid QSM file: edge refers to a non-existent node");

    QSMEdge data(edge);
    data.source = node_ids[edge.source - 1];
    data.target = node_ids[edge.target - 1];
    add_edge(data.source, data.target, data);
  }

  uint32_t tid = reader.get_treeid();
  if (tid < static_cast<uint32_t>(std::numeric_limits<int>::max()))
    id = static_cast<int>(tid);

  name = reader.get_treename();
  crs  = reader.get_crs();
}

} // namespace arbor::qsm
