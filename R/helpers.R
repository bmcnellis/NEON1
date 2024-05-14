#' @title Helper functions for NEON1
#'
#' @description Self-explanatory
#'
#' @rdname helpers
#' @name helpers
NULL
#' @rdname helpers
#' @export
logit_Warton <- function(y, epsilon = NULL) {

  # Warton & Hui (2011). Ecology 92(1) pp.3-10

  # Take as 0 the minimum y that isn't zero, or if y is large (mean > 0.75) then min(1 - y)
  # BUT:
  # "Users are encouraged to experiment with different values of e (epsiolon)"
  # AND:
  # "Whatever transform is used, it is important to check diagnostic plots"
  # ESPECIALLY:
  # "... no evidence of nonlinearity or heteroscedasticity. In small sample sizes (n < 30),
  # a normal probability plot should also be considered"
  # For GLMM, check that random effect is normally distributed

  if (all(na.omit(y < 1), na.omit(y > 0))) {
    return(log(y / (1 - y)))
  }

  if (is.null(epsilon)) {
    if (mean(y, na.rm = T) > 0.75) {
      y[which(y == 0)] <- min(1 - y, na.rm = T)
      y[which(y == 1)] <- 1 - min(1 - y, na.rm = T)
    } else {
      z <- min(y[which(y > 0)], na.rm = T)

      y[which(y == 0)] <- z
      y[which(y == 1)] <- 1 - z
    } # end ifelse
  } else {
    if (inherits(epsilon, 'numeric')) {
      stopifnot(length(epsilon) == 2)

      y[which(y == 0)] <- epsilon[1]
      y[which(y == 1)] <- epsilon[2]
    } else if (inherits(epsilon, 'function')) {
      y <- epsilon(y)
    } # end ifelsif
  } # end ifelse

  return(log(y / (1 - y)))

}
#' @rdname helpers
#' @export
check_residuals <- function(fit, dir, grp = NULL, nm = NULL) {

  require(glmmTMB)
  require(DHARMa)

  r0 <- resid(fit)
  s0 <- DHARMa::simulateResiduals(fit)

  nm0 <- paste0(ifelse(is.null(nm), deparse(substitute(fit)), nm), '_')
  f0 <- paste0(dir, '/check_residuals_', nm0, format(Sys.Date(), "%Y%m%d"), '.pdf')

  pdf(file = paste0(dir, '/check_residuals_',  deparse(substitute(fit)), '_', format(Sys.Date(), "%Y%m%d"), '.pdf'))

  plot(s0, rank = T)
  plot(s0, asFactor = T)

  if (!is.null(grp)) {

    stopifnot(is.list(grp))
    # plotResiduals should (more or less) have a mean at 0.5, boxes from 0.25-0.75, and whiskers from 0-1
    lapply(seq_along(grp), \(xx) DHARMa::plotResiduals(s0, form = grp[[xx]]))

  }

  dev.off()
  invisible()

}
#' @rdname helpers
#' @export
standardize_df <- function(df_in) {

  require(tidyr)
  require(vegan)

  df0 <- df_in
  df0 <- tidyr::pivot_wider(df0, names_from = taxonID, values_from = cover)

  cind <- c((ncol(df_in) - 2 + 1):ncol(df0))

  mat <- df0[, cind]
  mat <- as.matrix(mat)
  mat <- apply(mat, c(1, 2), \(xx) ifelse(is.na(xx), 0, xx))
  mat <- vegan::decostand(mat, method = 'total', MARGIN = 1)

  df1[, cind] <- mat
  df1 <- tidyr::pivot_longer(df1, cind, names_to = taxonID, values_to = cover)

  return(df1)
}
