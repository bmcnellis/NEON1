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
