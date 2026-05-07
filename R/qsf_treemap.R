# @file qsf_treemap.R
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

#' Convert a Quantitative Structure Forest (QSF) to an sf treemap
#'
#' Extracts global statistics from each Quantitative Structure Model (QSM)
#' contained in a Quantitative Structure Forest (QSF) object and returns
#' them as an \pkg{sf} point object. Each point corresponds to the root
#' position of an individual tree.
#'
#' @param qsf A Quantitative Structure Forest (QSF)
#'
#' @return An \pkg{sf} object with one feature per QSM. Attributes correspond
#'   to global statistics returned by \link{qsm_stats}, and geometry is
#'   defined by the root coordinates (\code{X_root}, \code{Y_root}, \code{Z_root}).
#'
#' @export
qsf_treemap = function(qsf)
{
  ans = Filter(function(x) inherits(x, "qsm"), qsf)
  ans = lapply(ans, function(x) qsm_stats(x)$stats_global)
  ans = data.table::rbindlist(ans)
  ans = ans[,-c("X_bh", "Y_bh", "Z_bh")]
  ans = sf::st_as_sf(ans, coords = c("X_root", "Y_root", "Z_root"))
  ans
}
