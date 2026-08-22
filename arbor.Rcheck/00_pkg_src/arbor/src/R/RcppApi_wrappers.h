/**
 * @file RcppApi_wrappers.h
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

#ifndef WRAPPERS_H
#define WRAPPERS_H

#ifdef USING_R

#include <Rcpp.h>

#include "arbor.h"

arbor::qsm::QSM as_qsm(Rcpp::DataFrame df);
arbor::qsm::QSF as_qsf(Rcpp::List x);
Rcpp::DataFrame as_dataframe(const arbor::qsm::QSM& qsm);
Rcpp::DataFrame as_dataframe(const PointCloud& cloud);

#endif

#endif
