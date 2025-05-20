#' Check if a binary is available in the system path
#'
#' @param bin_name Name or path of the binary to check
#' @return TRUE if the binary is available, otherwise stops with an error
#' @keywords internal
check_bin <- function(bin_name) {
  # Use 'which' on Unix/Linux/Mac or 'where' on Windows
  cmd <- if (.Platform$OS.type == "windows") {
    paste0("where ", bin_name, " 2>NUL")
  } else {
    paste0("which ", bin_name, " 2>/dev/null")
  }
  
  # Run command and check result
  result <- system(cmd, intern = TRUE)
  
  if (length(result) == 0) {
    stop("Required binary '", bin_name, "' not found in system path. ",
         "Please install it or provide the full path.")
  }
  
  return(TRUE)
}
