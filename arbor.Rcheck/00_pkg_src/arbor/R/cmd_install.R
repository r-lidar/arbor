#' Use some arbor features from a terminal
#'
#' Install the required files to be able to run some simple arbor commands from a terminal.
#'
#' @export
#' @md
install_cmd_tools <- function() 
{
  get_os <- function() {
    if (.Platform$OS.type == "windows") return("windows")
    sysname <- Sys.info()[["sysname"]]
    if (sysname == "Darwin") return("osx")
    if (sysname == "Linux") return("linux")
    return("unknown")
  }

  os <- get_os()

  script_path <- system.file("bash/arbor.R", package = "arbor")

  if (script_path == "") {
    stop("Could not find 'bash/arbor.R' in the installed arbor package.")
  }

  script_path <- normalizePath(script_path, mustWork = TRUE)
  bin_dir     <- normalizePath(dirname(script_path))

  if (os == "linux") 
  {
    target_dir <- path.expand("~/.local/bin")
    if (!dir.exists(target_dir)) {
      ans <- utils::askYesNo(sprintf("%s does not exist. Create it?", target_dir))
      if (!isTRUE(ans)) stop("Aborted: target directory missing.", call. = FALSE)
      dir.create(target_dir, recursive = TRUE)
    }

    # ensure the script itself is executable
    Sys.chmod(script_path, mode = "0755")

    target_link <- file.path(target_dir, "arbor")
    if (file.exists(target_link)) file.remove(target_link)
    if (!file.symlink(from = script_path, to = target_link)) {
      stop("Failed to create symlink. Check permissions.", call. = FALSE)
    }

    # warn if not on PATH
    path_dirs <- strsplit(Sys.getenv("PATH"), .Platform$path.sep)[[1]]
    if (!normalizePath(target_dir, mustWork = FALSE) %in% normalizePath(path_dirs, mustWork = FALSE)) {
      message(sprintf(
        "Note: %s is not on your PATH. Add this to your shell profile:\n  export PATH=\"%s:$PATH\"",
        target_dir, target_dir
      ))
    }

    message(sprintf("Symlink created at %s", target_link))
    return(invisible(target_link))
  }

  if (os == "windows") 
  {
    windows_manual_instructions <- function(bin_dir) 
    {
      paste0(
        "To run 'arbor' from anywhere in the Command Prompt, add its location to PATH manually:\n\n",
        "  1. In the Windows search bar, type: 'env'\n",
        "  2. Click: Edit the system environment variables\n",
        "  3. Go to the 'Advanced' tab, click 'Environment variables'\n",
        "  4. Select 'Path' under User variables, click Edit, then New\n",
        "  5. Add this path: ", bin_dir, "\n\n",
        "More details: https://www.eukhost.com/kb/how-to-add-to-the-path-on-windows-10-and-windows-11/\n"
      )
    }
    rscript_exe <- Sys.which("Rscript")
    if (rscript_exe == "") stop("Cannot find Rscript.exe in the system PATH.", call. = FALSE)

    batch_file_path <- sub("\\.R$", ".bat", script_path)

    batch_content <- c(
      "@echo off",
      paste0('"', rscript_exe, '" "', script_path, '" %*')
    )

    tryCatch({
      writeLines(batch_content, con = batch_file_path)
    }, error = function(e) {
      stop("Could not write .bat file. You may need to run R as Administrator.\nError: ", e$message, call. = FALSE)
    })

    # check whether bin_dir is already on PATH
    path_dirs <- strsplit(Sys.getenv("PATH"), ";")[[1]]
    on_path <- normalizePath(bin_dir, mustWork = FALSE) %in%
      normalizePath(path_dirs, mustWork = FALSE)

    if (on_path) 
    {
      message("'arbor' is ready to use. Its directory is already on your PATH.")
      return(invisible(batch_file_path))
    }

    ans <- utils::askYesNo(sprintf(
      "Add %s to your user PATH now? (This changes your permanent user environment.)",
      bin_dir
    ))

    if (isTRUE(ans)) 
    {
      # setx persists to the user's registry-based PATH; does not affect current session
      current_user_path <- system2("reg", c(
        "query", "HKCU\\Environment", "/v", "Path"
      ), stdout = TRUE, stderr = FALSE)

      # best-effort append; if reg query fails, just fall back to setx with new dir
      new_path <- paste0(bin_dir, ";", Sys.getenv("PATH"))
      exit_code <- system2("setx", c("PATH", shQuote(new_path)))

      if (exit_code == 0) 
      {
        message(sprintf(
          "PATH updated. Restart your terminal, then 'arbor' will be available.\n(Directory added: %s)",
          bin_dir
        ))
      } 
      else {
      
        cat(windows_manual_instructions(bin_dir))
      }
    } 
    else
    {
      cat(windows_manual_instructions(bin_dir))
    }

    return(invisible(batch_file_path))
  }

  if (os == "osx") 
  {
    # Try user-writable locations first, avoid requiring sudo when possible
    candidates <- c(
      path.expand("~/.local/bin"),
      "/opt/homebrew/bin",   # Apple Silicon Homebrew
      "/usr/local/bin"       # Intel Homebrew / traditional, often needs sudo
    )

    # ensure script is executable regardless of which path we use
    Sys.chmod(script_path, mode = "0755")

    target_dir <- NULL
    for (d in candidates) {
      if (dir.exists(d) && file.access(d, mode = 2) == 0) {
        target_dir <- d
        break
      }
    }

    if (is.null(target_dir)) {
      # nothing writable found — offer to create ~/.local/bin
      fallback <- path.expand("~/.local/bin")
      ans <- utils::askYesNo(sprintf(
        "No writable bin directory found. Create %s?", fallback
      ))
      if (!isTRUE(ans)) {
        stop("Aborted: no writable target directory.", call. = FALSE)
      }
      dir.create(fallback, recursive = TRUE)
      target_dir <- fallback
    }

    target_link <- file.path(target_dir, "arbor")
    if (file.exists(target_link)) file.remove(target_link)

    success <- file.symlink(from = script_path, to = target_link)

    if (!success) {
      stop(sprintf(
        "Failed to create symlink at %s. If this is /usr/local/bin, you may need:\n  sudo ln -sf %s %s",
        target_link, shQuote(script_path), shQuote(target_link)
      ), call. = FALSE)
    }

    path_dirs <- strsplit(Sys.getenv("PATH"), .Platform$path.sep)[[1]]
    on_path <- normalizePath(target_dir, mustWork = FALSE) %in%
      normalizePath(path_dirs, mustWork = FALSE)

    if (!on_path) {
      message(sprintf(
        "Note: %s is not on your PATH. Add this to ~/.zshrc or ~/.bash_profile:\n  export PATH=\"%s:$PATH\"",
        target_dir, target_dir
      ))
    }

    message("Symlink successfully created at ", target_link)
    return(invisible(target_link))
  }

  stop("Operating system not supported or detected.")
}