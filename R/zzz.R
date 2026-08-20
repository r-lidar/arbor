# @file zzz.R
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

#' @useDynLib arbor, .registration = TRUE
.datatable.aware = TRUE

.onLoad <- function(libname, pkgname)
{
  arbor_parameters_default <<- default_arbor_params_cpp()

  #occ <- arbor_parameters_default
  #occ$path_finder$k_neighborhood_connectivity = 30
  #occ$path_finder$max_gap = 1
  #occ$qsm$step = 0.2
  #occ$qsm$cl_dist = 0.2
  #arbor_parameters_occlusion <<- occ
}

.onAttach <- function(libname, pkgname)
{
  check_update()
  threads <- data.table::getDTthreads()
  if (threads == 1) message("This version of Arbor has no multicore support.")
}

# Check if the package has more recent version
check_update = function()
{
  msg <- NULL

  last <- get_latest_version()
  curr <- utils::packageVersion("arbor")
  new_version = !is.null(last) && last[1] > curr[1]
  new_version = !is.null(last) && last[1] > curr
  dev_version = is_dev_version(curr)

  # nocov start
  if (new_version)
  {
    if (dev_version)
      msg = paste("arbor", last, "is now available. You are using", curr, "(unstable) \ninstall.packages('arbor', repos = 'https://r-lidar.r-universe.dev')")
    else
      msg = paste("arbor", last, "is now available. You are using", curr, "\ninstall.packages('arbor', repos = 'https://r-lidar.r-universe.dev')")
  }
  else if (dev_version)
  {
    msg = paste("arbor", curr, "is an unstable development version")
  }

  if (!is.null(msg) & interactive()) packageStartupMessage(msg)
  # nocov end

  return(NULL)
}

get_latest_version = function()
{
  nullcon = NULL

  ans <- tryCatch(
  {
    nullcon <- file(nullfile(), open = "wb")
    sink(nullcon, type = "message")
    res <- utils::old.packages(repos = "https://r-lidar.r-universe.dev")
    sink(type = "message")
    close(nullcon)
    res
  },
  error = function(e)
  {
    sink(NULL, type = "message") # nocov
    close(nullcon) # nocov
    return(NULL) # nocov
  })

  if (is.null(ans)) return(NULL) # nocov

  ind = which(ans[,1] == "arbor")

  if (length(ind) == 0) return(NULL)

  version <- ans[ind, 5] # nocov
  version <- package_version(version) # nocov
  return(version) # nocov
}

is_dev_version = function(version)
{
  class(version) <- "list"
  version = version[[1]]
  return(length(version) == 4)
}
