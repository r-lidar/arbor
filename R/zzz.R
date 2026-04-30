#' @useDynLib arbor, .registration = TRUE
.datatable.aware = TRUE

#' Use some arbor features from a terminal
#'
#' Install the required files to be able to run some simple arbor commands from a terminal.
#'
#' @export
#' @md
install_cmd_tools <- function() {

  # 1. Robust OS detection
  get_os <- function() {
    if (.Platform$OS.type == "windows") return("windows")
    sysname <- Sys.info()[["sysname"]]
    if (sysname == "Darwin") return("osx")
    if (sysname == "Linux") return("linux")
    return("unknown")
  }

  os <- get_os()

  # 2. Locate the source script correctly (Post-install location)
  # Note: 'inst/bash' moves to 'bash' in the installed package
  script_path <- system.file("bash/arbor.R", package = "arbor")

  if (script_path == "") {
    stop("Could not find 'bash/arbor.R' in the installed arbor package.")
  }

  script_path <- normalizePath(script_path, mustWork = TRUE)
  bin_dir     <- normalizePath(dirname(script_path))

  # 3. Platform Specific Logic
  if (os == "linux")
  {
    # Check if target directory exists
    target_dir <- path.expand("~/.local/bin")
    if (!dir.exists(target_dir)) {
      stop(sprintf("The target directory '%s' does not exist. Please create it or add it to your PATH.", target_dir))
    }

    target_link <- file.path(target_dir, "arbor")

    # Remove existing link if needed to avoid error
    if (file.exists(target_link)) file.remove(target_link)

    success <- file.symlink(from = script_path, to = target_link)

    if (success) {
      message(sprintf("Symlink created at %s", target_link))
    } else {
      stop("Failed to create symlink. Check permissions.")
    }
    return(invisible())
  }

  if (os == "windows")
  {
    rscript_exe <- Sys.which("Rscript")
    if (rscript_exe == "") stop("Cannot find Rscript.exe in the system PATH.")

    # Create .bat file next to the .R file
    # Warning: This still might fail if the library folder is read-only
    batch_file_path <- sub("\\.R$", ".bat", script_path)

    # Use standard Windows batch quoting
    batch_content <- c(
      "@echo off",
      paste0('"', rscript_exe, '" "', script_path, '" %*')
    )

    tryCatch({
      writeLines(batch_content, con = batch_file_path)
    }, error = function(e) {
      stop("Could not write .bat file. You may need to run R as Administrator.\nError: ", e$message)
    })

    cat("To run the 'qsm' command from anywhere in the Command Prompt, you need to add its location to your system's PATH:

  1. In the windows search bar type: 'env'
  2. Click on: Edit the system environment variables
  3. Go to the 'Advanced' tab and click on 'Environment variables'
  4. Click on 'Path' and edit the settings.
  5. Add this path:", bin_dir, "

For more details see https://www.eukhost.com/kb/how-to-add-to-the-path-on-windows-10-and-windows-11/")


    return(invisible())
  }

  if (os == "osx")
  {
    # OSX usually requires sudo for /usr/local/bin, or the user might use ~/bin
    target_link <- "/usr/local/bin/arbor"

    # file.symlink is cleaner, but triggers permission errors silently sometimes.
    # We will use system command to capture stderr if it fails.
    cmd <- sprintf("ln -sf %s %s", shQuote(script_path), shQuote(target_link))
    exit_code <- system(cmd, ignore.stderr = FALSE)

    if (exit_code != 0) {
      stop(sprintf("Failed to link to %s. You likely need 'sudo' privileges or the directory does not exist.", target_link))
    }

    message("Symlink successfully created at ", target_link)
    return(invisible())
  }

  stop("Operating system not supported or detected.")
}

.onLoad <- function(libname, pkgname)
{
  default_arbor_parameters <<- default_arbor_params_cpp()
}
