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
  qsf          Compute QSMs of segmented trees
  report       Produce an html report
  pipeline     Full pipeline

Options:
  -h, --help   Show help for a specific command

Examples:
  arbor segment plot.laz -no-dtm
  arbor qsf plot_output/plot_trees.laz
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

# Validate a command line: every token must be the leading positional input, a known
# boolean flag, or a known flag followed by a value.
check_args <- function(args, value_flags, bool_flags)
{
  i <- if (length(args) > 0 && !startsWith(args[1], "-")) 2L else 1L
  while (i <= length(args))
  {
    a <- args[i]
    if (a %in% value_flags)
    {
      if (i == length(args)) stop(paste0("Option '", a, "' requires a value."), call. = FALSE)
      i <- i + 2L
    }
    else if (a %in% bool_flags)
    {
      i <- i + 1L
    }
    else if (startsWith(a, "-"))
    {
      stop(paste0("Unknown option '", a, "'. Use -h for help."), call. = FALSE)
    }
    else
    {
      stop(paste0("Unexpected argument '", a, "'. Use -h for help."), call. = FALSE)
    }
  }
  invisible(TRUE)
}
