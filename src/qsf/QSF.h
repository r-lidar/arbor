/**
 * @file QSF.h
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

#ifndef QSF_H
#define QSF_H

#include "QSM.h"
#include <string>
#include <unordered_map>

namespace arbor::qsm {

class QSF
{
public:
  QSF() = default;
  void add_qsm(int id, const QSM& q);
  void write(const std::string& dir, const std::string& format, bool binary = true) const;
  const std::unordered_map<int, QSM>& get_qsm_map() const { return qsm_; }

private:
  std::unordered_map<int, QSM> qsm_;
};

}

#endif
