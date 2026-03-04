# Main Usage Help
usage_main <- function() {
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
