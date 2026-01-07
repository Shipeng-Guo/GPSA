#' @keywords internal
"_PACKAGE"

## Package-level options and environment
.gpsa_env <- new.env(parent = emptyenv())

.onLoad <- function(libname, pkgname) {
  ## Initialize default options if not set
  op <- options()
  op.gpsa <- list(
    GPSA.h5 = NULL,
    GPSA.verbose = TRUE
  )
  toset <- !(names(op.gpsa) %in% names(op))
  if (any(toset)) options(op.gpsa[toset])
  invisible()
}

#' Get GPSA H5 file path
#'
#' Returns the HDF5 file path from options, environment variable, or parameter.
#' Priority: parameter > options("GPSA.h5") > Sys.getenv("GPSA_H5")
#'
#' @param h5file Optional explicit path to HDF5 file
#' @return Character path to HDF5 file
#' @export
gpsa_get_h5path <- function(h5file = NULL) {
  if (!is.null(h5file) && nzchar(h5file)) return(h5file)
  
  opt_path <- getOption("GPSA.h5")
  if (!is.null(opt_path) && nzchar(opt_path)) return(opt_path)
  
  env_path <- Sys.getenv("GPSA_H5", unset = "")
  if (nzchar(env_path)) return(env_path)
  
  stop("HDF5 path not specified. Set via:\n",
       "  - Function parameter: h5file = '/path/to/file.h5'\n",
       "  - R options: options(GPSA.h5 = '/path/to/file.h5')\n",
       "  - Environment: Sys.setenv(GPSA_H5 = '/path/to/file.h5')")
}

#' Set GPSA H5 file path
#'
#' Convenience function to set the default HDF5 file path.
#'
#' @param h5file Path to HDF5 file
#' @return Previous value (invisibly)
#' @export
gpsa_set_h5path <- function(h5file) {
  old <- getOption("GPSA.h5")
  options(GPSA.h5 = h5file)
  invisible(old)
}

