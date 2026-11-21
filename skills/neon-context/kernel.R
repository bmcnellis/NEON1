NEON_DATA_ROOT <- "/media/bem/data/NEON1"

neon_path <- function(..., root = NULL) {
  if (is.null(root)) root <- NEON_DATA_ROOT
  file.path(root, ...)
}

neon_load_model <- function(which = 0, diag_only = F, root = NULL) {
  stopifnot(which %in% c(0, 1))
  nm <- if (diag_only) "m_p_diag" else "m_p_mod"
  p <- neon_path("results", "model_results", paste0(nm, which, ".rda"), root = root)
  if (!file.exists(p)) stop(paste("not found:", p))
  e <- new.env(parent = emptyenv())
  load(p, envir = e)
  as.list(e)
}
