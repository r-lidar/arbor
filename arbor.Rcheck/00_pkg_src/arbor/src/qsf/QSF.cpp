/**
 * @file QSF.cpp
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

#include "QSF.h"

#include <string>
#include <filesystem>
#include <stdexcept>

namespace fs = std::filesystem;

namespace arbor::qsm {

void QSF::add_qsm(const QSM& q)
{
  if (q.nodes().size() <= 1) return;
  qsm_[q.id] = q;
}

void QSF::write(const std::string& dir, const std::string& format, bool binary) const
{
  if (dir.empty())
    throw std::invalid_argument("QSF::write: output directory is empty");

  if (format.empty())
    throw std::invalid_argument("QSF::write: format is empty");

  fs::path base_dir(dir);

  // Create base directory if needed
  if (!fs::exists(base_dir))
  {
    fs::create_directories(base_dir);
  }
  else if (!fs::is_directory(base_dir))
  {
    throw std::runtime_error("QSF::write: path exists but is not a directory");
  }

  // Optional subfolder per format
  fs::path out_dir = base_dir / format;

  if (!fs::exists(out_dir))
    fs::create_directories(out_dir);

  // Write each QSM
  for (const auto& [key, qsm] : qsm_)
  {
    fs::path filename;
    if (qsm.name.empty())
      filename = out_dir / (std::to_string(qsm.id) + "." + format);
    else
      filename = out_dir / (qsm.name + "." + format);

    qsm.write(filename.string(), binary);
  }
}

}
