# @file cmd_usage.R
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

# Main Usage Help
cmd_usage <- function() {
  cat("
Usage:
  arbor <command> [arguments]

Commands:
  segment      Segment a point cloud (LAS/LAZ)
  qsm          Run QSM on a folder or file
  report       Produce a pdf report

Options:
  -h, --help   Show help for a specific command

Examples:
  arbor segment plot.laz --no-dtm
  arbor qsm ./my_trees/ -ncores 4
")
  quit(save = "no", status = 0)
}

# Parse a value associated with a flag (e.g., --buffer 5)
get_arg <- function(args, flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 1 && idx < length(args)) {
    return(args[idx + 1])
  }
  default
}

# Check if a flag exists (boolean)
has_flag <- function(args, flag) {
  flag %in% args
}

# Fail safely with a message
fail <- function(msg) {
  cat("Error:", msg, "\n", file = stderr())
  quit(status = 1)
}
