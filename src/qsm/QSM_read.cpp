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

  for (const auto& node : reader.nodes())
    add_node(node);

  for (const auto& edge : reader.edges())
    add_edge(static_cast<NodeID>(edge.source), static_cast<NodeID>(edge.target), QSMEdge(edge));

  uint32_t tid = reader.get_treeid();
  if (tid < static_cast<uint32_t>(std::numeric_limits<int>::max()))
    id = static_cast<int>(tid);

  name = reader.get_treename();
  crs  = reader.get_crs();
}

} // namespace arbor::qsm
