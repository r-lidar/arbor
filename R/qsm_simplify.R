# @file qsm_simplify.R
# Project: Arbor
# 
# Copyright (C) 2026 Jean-Romain Roussel (r-lidar) <info @ r-lidar.com>
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# Used to simply AdTree skeleton by MRNF only
qsm_simplify <- function(qsm, max_length = 0.3)
{
  stop("Function disabled. Why is it in use?")
  # qsm$cyl_length = sqrt((qsm$startX - qsm$endX)^2 + (qsm$startY- qsm$endY)^2 + (qsm$startZ - qsm$endZ)^2)
  # qsm = qsm_simplify_cpp(qsm, max_length)
  # qsm = qsm_topology(qsm)
  # qsm = qsm_architecture_cpp(qsm)
  # qsm
}
