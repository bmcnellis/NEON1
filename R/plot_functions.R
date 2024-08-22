#' @title Plot functions for NEON1 analysis
#'
#' @description NA
#'
#' @rdname plot_functions
#' @name plot_functions
NULL
#' @rdname plot_functions
#' @export
HighResTiff <- function(plot_obj, file, width_in, height_in, resolution_dpi, ...) {

  if (inherits(plot_obj, 'ggplot')) {
    tiff(filename = file, width = width_in, height = height_in, units = 'in', res = resolution_dpi)
    print(plot_obj)
    dev.off()
  } else if (inherits(plot_obj, 'list')) {
    tiff(filename = file, width = width_in, height = height_in, units = 'in', res = resolution_dpi)
    Multiplot(plotlist = plot_obj, ...)
    dev.off()
  }
  invisible()

}
#' @rdname plot_functions
#' @export
lm_eqn <- function(df0, xvar, yvar){

  df <- data.frame(x = df0[[xvar]], y = df0[[yvar]])

  m <- lm(y ~ x, df)
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq))
}
#' @rdname plot_functions
#' @export
Multiplot <- function(..., plotlist = NULL, file, cols = 1, layout = NULL) {
  # BEM note: this is pretty much verbatim copied from the R cookbook
  #
  # Multiple plot function
  #
  # ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
  # - cols:   Number of columns in layout
  # - layout: A matrix specifying the layout. If present, 'cols' is ignored.
  #
  # If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
  # then plot 1 will go in the upper left, 2 will go in the upper right, and
  # 3 will go all the way across the bottom.
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])
    #plot(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
      #plot(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
#' @rdname plot_functions
#' @export
NMDS_screeplot <- function(mat, kmax = 10) {
  # see Figure 16.3 in:
  #     McCune, B., Grace, J.B., 2002. Analysis of Ecological Communities.
  #     MjM Software Design, Gleneden Beach, OR.

  require(vegan)
  require(ggplot2)

  stress <- sapply(seq(kmax), \(xx) {
    vegan::metaMDS(mat, k = xx, trace = 0)$stress
  })

  d0 <- data.frame(Dimension = seq(kmax), Stress = stress)

  p0 <- ggplot(data = d0, aes(x = Dimension, y = Stress)) +
    geom_line() +
    geom_point() +

    theme_bw() +
    theme(
      axis.text.x = element_text(color = 'black'),
      axis.text.y = element_text(color = 'black')
    )

  return(p0)

}
#' @rdname plot_functions
#' @export
c25 <- function(option = 1) {
  c0 <- c(
    "dodgerblue2", "#E31A1C",   "green4",   "#6A3D9A",     "#FF7F00",
    "black",       "gold1",     "skyblue2", "#FB9A99",     "palegreen2",
    "#CAB2D6",     "#FDBF6F",   "gray70",   "khaki2",      "maroon",
    "orchid1",     "deeppink1", "blue1",    "steelblue4",  "darkturquoise",
    "green1",      "yellow4",   "yellow3",  "darkorange4", "brown"
  )

  if (option == 1) {
    return(c0)
  } else if (option == 2) {
    # blues/greens
    return(c0[c(1, 3, 8, 10, 18, 19, 20, 21, 22)])
  } else {
    stop('bad option')
  }
}
plot_beta <- function(hM, post) {

  require(ggplot2)

  pp <- ifelse(abs(post$support) < 0.9, NA, post$mean)
  pp <- pp[, colSums(is.na(pp)) < 7]
  pp <- data.frame(hM$covNames, pp)
  pp <- tidyr::pivot_longer(pp, cols = 2:ncol(pp), values_to = 'effect', names_to = 'spp')
  pp <- pp[complete.cases(pp), ]
  pp$spp <- sapply(strsplit(pp$spp, '_'), \(xx) xx[1])
  pp <- dplyr::group_by(pp, hM.covNames, spp)
  pp <- dplyr::summarise(pp, effect = mean(effect, na.rm = T), .groups = 'drop')
  pp <- pp[which(!pp$hM.covNames %in% c('(Intercept)', 'elev:TWI', 'time_since_fence')), ]

  pp$hM.covNames <- ifelse(pp$hM.covNames == 'age_median', 'Flow age', pp$hM.covNames)
  pp$hM.covNames <- ifelse(pp$hM.covNames == 'elev', 'Elevation', pp$hM.covNames)
  pp$hM.covNames <- ifelse(pp$hM.covNames == 'pai', 'Light stress', pp$hM.covNames)
  pp$hM.covNames <- ifelse(pp$hM.covNames == 'time_since_fence', 'Pig fence age', pp$hM.covNames)
  pp$hM.covNames <- ifelse(pp$hM.covNames == 'TWI', 'Topographic run-in', pp$hM.covNames)

  ggplot(data = pp, aes(x = spp, y = effect)) +
    geom_point(size = 2.5) +
    facet_wrap(~ hM.covNames, scales = 'free_x') +

    theme_bw() +
    theme(
      axis.text.x = element_text(color = 'black', angle = 45, hjust = 1, size = 6),
      axis.text.y = element_text(color = 'black')
    ) +
    labs(y = 'Relative effect on plant presence', x = '')

}
#' @rdname plot_functions
#' @export
plot_vp <- function(hM, VP) {
  require(Hmsc)
  require(ggplot2)

  labs0 <- c(
    '(sample)', 'Flow age', 'Forest type', 'Grazing disturbance', 'Elevation',
    'Harvest disturbance', 'PAI', 'Ungulate disturbance', 'Years ungulate-free'
  )
  vals0 <- NEON1::c25(2)

  vpd <- VP$vals
  vpd <- data.frame(var = row.names(vpd), vpd)
  vpd <- tidyr::pivot_longer(vpd, cols = -1, names_to = 'spp', values_to = 'vp')
  vpd$vp <- round(vpd$vp * 100, 4)
  vpd$var <- ifelse(vpd$var == 'Random: plotDate', 'AAA_plotDate', vpd$var)

  plot0 <- ggplot(data = vpd, aes(x = spp, y = vp, fill = var)) +
    geom_bar(position = 'stack', stat = 'identity') +

    scale_fill_manual(values = vals0, labels = labs0, name = 'Variable') +

    theme_bw() +
    theme(
      axis.text.x = element_text(color = 'black', angle = 45, hjust = 1),
      axis.text.y = element_text(color = 'black')
    ) +
    labs(x = '', y = 'Proportion of explained variance')

  return(plot0)

}
