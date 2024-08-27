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

  on.exit(try(dev.off(), silent = T))

  require(glmmTMB)
  require(DHARMa)

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
#' @rdname helpers
#' @export
div_to_long <- function(div, type = 'plant', add_zeros = T) {

  # imports
  require(tidyr)
  require(dplyr)

  # magics
  df_ind <- 4
  c0 <- c(
    'siteID', 'decimalLatitude', 'decimalLongitude', 'elevation', 'plotType',
    'plotID', 'subplotID', 'endDate', 'divDataType', 'taxonID', 'scientificName',
    'nativeStatusCode', 'otherVariables', 'percentCover'
  )
  ntax <- 'taxonID'
  ncov <- 'percentCover'
  t0 <- ifelse(type == 'plant', 'plantSpecies', 'otherVariables')
  cov_NA <- 0.5
  df_ind <- 4

  # process
  df_out <- div[[df_ind]] |>
    dplyr::select(c0) |>
    dplyr::filter(divDataType == t0) |>
    # put in a trace value for NA
    dplyr::mutate(percentCover = ifelse(is.na(percentCover), cov_NA, percentCover)) |>
    # remove non-uniquely-identified subplots
    dplyr::group_by(dplyr::across(c(-percentCover))) |>
    dplyr::summarize(percentCover = sum(percentCover)) |>
    dplyr::ungroup()

  if (add_zeros) {

    df_out <- df_out |>
      tidyr::pivot_wider(names_from = ntax, values_from = ncov, values_fill = list(percentCover = 0)) |>
      tidyr::pivot_longer(cols = !c0[!c0 %in% c(ntax, ncov)], names_to = ntax, values_to = ncov)
  }

  # returns
  return(df_out)

}
#' @rdname helpers
#' @export
date_match <- function(x, y, p) {
  # finds the closest date match in dhp for a vector of plant collection dates

  # xx is data_df$endDate, yy is NEON1::dhp
  # requires grouped dataframe
  stopifnot(length(p) == 1)

  if (p %in% unique(y$plot)) {

    # get only dates within the active group
    y0 <- y[which(y$plot == p), 'date']
    # apply along the length of the data dates and find match within dhp
    z <- sapply(x, \(xx) y0[which.min(abs(y0 - xx))])
    # reconvert to date
    z <- as.Date(z, origin = '1970-01-01')

  } else {

    z <- rep(NA, length(x))

  }

  return(z)

}
#' @rdname helpers
#' @export
sph <- function(vec) {

  v <- sapply(strsplit(vec, '\\+/-'), \(xx) xx[1])
  v <- ifelse(v == 'NULL', NA, v)
  v <- as.numeric(v)

  e <- sapply(strsplit(vec, '\\+/-'), \(xx) xx[2])
  e <- as.numeric(e)

  data.frame(val = v, err = e)

}
#' @rdname helpers
#' @export
order_to_class <- function() {
  # https://en.wikipedia.org/wiki/List_of_plant_orders
  o_t_c <- data.frame(
    order = c(
      'Apiales', 'Aquifoliales', 'Arecales', 'Asparagales', 'Asterales',
      'Caryophyllales', 'Celastrales', 'Cornales', 'Cyatheales', 'Ericales',
      'eupolypod II', 'Fabales', 'Gentianales', 'Geraniales', 'Gleicheniales',
      'Hymenophyllales', 'Hypnales', 'Lamiales', 'Liliales', 'Lycopodiales',
      'Malpighiales', 'Malvales', 'Marattiales', 'Myrtales', 'Ophioglossales',
      'Pandanales', 'Piperales', 'Poales', 'Polypodiales', 'Polypodiales-eupolypod I',
      'Polypodiales-eupolypod II', 'Rosales', 'Santalales', 'Sapindales'
    ),
    superorder = c(
      'Asteranae', 'Asteranae', 'Commelinids', 'Orchidanae', 'Asteranae',
      'Caryophyllanae', 'Rosanae', 'Cornanae', 'Cyatheales', 'Ericanae',
      'Polypodiales', 'Rosanae', 'Lamianae', 'Myrtanae', 'Gleicheniales',
      'Hymenophyllales', 'Hypnanae', 'Lamianae', 'Lilianae', 'Lycopodiales',
      'Rosanae', 'Myrtanae', 'Marattiales', 'Myrtanae', 'Ophioglossales',
      'Pandananae', 'Piperales', 'Commelinids', 'Polypodiales', 'Polypodiales',
      'Polypodiales', 'Rosanae', 'Santalanae ', 'Myrtanae'
    ),
    subclass = c(
      'Asteridae', 'Asteridae', 'Liliidae', 'Liliidae', 'Asteridae',
      'Asteridae', 'Rosidae', 'Asteridae', 'Polypodiidae', 'Asteridae',
      'Polypodiidae', 'Rosidae', 'Asteridae', 'Rosidae', 'Polypodiidae',
      'Polypodiidae', 'Bryidae', 'Asteridae', 'Liliidae', 'Lycopodiidae',
      'Rosidae', 'Rosidae', 'Marattiopsida', 'Rosidae', 'Ophioglossidae',
      'Liliidae', 'Magnoliidae', 'Liliidae', 'Polypodiidae', 'Polypodiidae',
      'Polypodiidae', 'Rosidae', 'Asteridae', 'Rosidae'
    ),
    class = c(
      'Magnoliopsida', 'Magnoliopsida', 'Magnoliopsida', 'Magnoliopsida', 'Magnoliopsida',
      'Magnoliopsida', 'Magnoliopsida', 'Magnoliopsida', 'Polypodiopsida', 'Magnoliopsida',
      'Polypodiopsida', 'Magnoliopsida', 'Magnoliopsida', 'Magnoliopsida', 'Polypodiopsida',
      'Polypodiopsida', 'Bryopsida', 'Magnoliopsida', 'Magnoliopsida', 'Lycopodiopsida',
      'Magnoliopsida', 'Magnoliopsida', 'Polypodiopsida', 'Magnoliopsida', 'Polypodiopsida',
      'Magnoliopsida', 'Magnoliopsida', 'Magnoliopsida', 'Polypodiopsida', 'Polypodiopsida',
      'Polypodiopsida', 'Magnoliopsida', 'Magnoliopsida', 'Magnoliopsida'
    )
  )

  return(o_t_c)
}
#' @rdname helpers
#' @export
ess_as_df <- function(ess, param) {
  require(Hmsc)
  require(coda)

  edf <- data.frame(term = names(ess), ess = as.numeric(ess))
  edf$param <- sapply(strsplit(names(ess), ','), \(xx) xx[1])
  edf$spp <- sapply(strsplit(names(ess), ','), \(xx) xx[2])
  edf$spp <- sapply(strsplit(edf$spp, ' '), \(xx) xx[2])

  if (param == 'beta') {
    edf$param <- gsub('B\\[', '', edf$param)
    edf$param <- sapply(strsplit(edf$param, ' '), \(xx) xx[1])

    ret <- edf
  }

  return(ret)

}
#' @rdname helpers
#' @export
posterior_from_coda <- function(mcmc, term, prob, average = F, drop_ns = F) {
  require(coda)

  mb <- summary(mcmc[[term]], quantiles = NULL)[[1]]
  mb <- data.frame(term = row.names(mb), mean = mb[, 1], row.names = NULL)

  hp <- coda::HPDinterval(mcmc[[term]], prob = prob)

  for (i in seq_along(hp)) colnames(hp[[i]]) <- paste(colnames(hp[[i]]), i, sep = '_')
  hp <- do.call('cbind', hp)
  hp <- data.frame(term = row.names(hp), hp, row.names = NULL)

  hpd <- dplyr::left_join(mb, hp, by = 'term')

  if (average) {
    lower <- rowMeans(hpd[, grepl('lower', colnames(hpd))], na.rm = T)
    upper <- rowMeans(hpd[, grepl('upper', colnames(hpd))], na.rm = T)
    hpd <- data.frame(term = hpd$term, mean = hpd$mean, lower, upper)

    if (drop_ns) {
      hpd <- hpd[!(hpd$lower <= 0 & hpd$upper >= 0), ]
    }
  }

  spp <- sapply(strsplit(hpd$term, ' '), \(xx) xx[3])
  var <- sapply(strsplit(hpd$term, ' '), \(xx) xx[1])
  var <- sapply(strsplit(var, '\\['), \(xx) xx[2])

  hpd$spp <- spp
  hpd$var <- var

  return(hpd)

}
#' @rdname helpers
#' @export
NEON_div_full_400 <- function() {
  # gets all subplots necessary to pool presence/absence by the full 400-m quad

  # 31.1.10 + 31.1.1 + 34.4.10 + 34.4.1 + 32.2.10 + 32.2.1 + 32.4.10 + 32.4.1 +
  # 40.1.10 + 40.1.1 + 40.3.10 + 40.3.1 + 41.1.10 + 41.1.1 + 41.4.10 + 41.4.1
  f0 <- c(
    '31.1.10', '31.1.1','34.4.10', '34.4.1', '32.2.10', '32.2.1', '32.4.10', '32.4.1',
    '40.1.10', '40.1.1', '40.3.10', '40.3.1', '41.1.10', '41.1.1', '41.4.10', '41.4.1'
  )

  return(f0)

}
