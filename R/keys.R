#' @title Key functions for NEON data
#'
#' @description Sometimes you need to key abbreviations and such!
#'
#' @rdname keys
#' @name keys
NULL
#' @rdname keys
#' @export
key_field_domain_id <- function(id) {

  k0 <- data.frame(
    id = seq(20),
    name = c(
      'Northeast', 'Mid-Atlantic', 'Southeast', 'Atlantic Neotropical',
      'Great Lakes', 'Prairie Peninsula', 'Appalachians & Cumberland Plateau', 'Ozarks Complex',
      'Northern Plains', 'Central Plains', 'Southern Plains',
      'Northern Rockies', 'Southern Rockies & Colorado Plateau',
      'Desert Southwest', 'Great Basin',
      'Pacific Northwest', 'Pacific Southwest',
      'Tundra', 'Taiga', 'Pacific Tropical'
    ),
    row.names = NULL
  )

  if (missing(id)) {
    return(k0)
  } else {
    if (!all(id %in% seq(20))) {
      stop('id(s) must be a domain number 1-20')
    } else {
      return(k0[id, 'name'])
    }
  }
}
#' @rdname keys
#' @export
key_age_range <- function(age_range, current_year = 2023) {

  # removing box2 in outline shapefile changed levels
  age_char <- c(
    #"A.D. 1984",
    "A.D. 1942", #"A.D. 1852", "A.D. 1410-1460",
    "200-750 yr", #"400-750 yr",
    "750-1,500 yr",
    "1,500-3,000 yr", "3,000-5,000 yr", "5,000-11,000 yr", "11,000-30,000 yr"
  )
  age_key <- data.frame(age_char, age_num = c(
    current_year - 1984, current_year - 1942, current_year - 1852, #current_year - ((1460 - 1410) / 2 + 1410),
    ((750 - 200) / 2) + 200, #((750 - 400) / 2) + 400,
    ((1500 - 750) / 2) + 750,
    ((3000 - 1500) / 2) + 1500, ((5000 - 3000) / 2) + 3000, ((11000 - 5000) / 2) + 5000, ((30000 - 11000) / 2) + 11000
  ))
  age_key$age_num <- round(age_key$age_num)

  stopifnot(is.character(age_range), all(age_range %in% age_char))

  age_num <- age_key$age_num[na.omit(match(age_range, age_key$age_char))]

  return(age_num)

}
#' @rdname keys
#' @export
relevel_age <- function(age_range) {

  age_fact <- ordered(age_range, levels = c(
    #"A.D. 1984",
    "A.D. 1942", #"A.D. 1852", "A.D. 1410-1460",
    "200-750 yr", #"400-750 yr",
    "750-1,500 yr",
    "1,500-3,000 yr", "3,000-5,000 yr", "5,000-11,000 yr", "11,000-30,000 yr"
  ))


}
#' @rdname keys
#' @export
key_plots_by_type <- function(plot_id) {

  plot_key <- data.frame(
    id <- seq(56),
    veg <- c(
      # in groups of 10
      rep(TRUE, 10),
      '', '', '', '', '', '', '', '', '', '',
      '', '', '', '', '', '', '', '', '', '',
      '', '', '', '', '', '', '', '', '', '',
      '', '', '', '', '', '', '', '', '', '',
      '', '', '', '', '', ''
    ),
    bird <- c(
      rep(TRUE, 10),
    ),
    location <- c(
      rep('distributed', 10),

    )
  )
}
#' @rdname keys
#' @export
key_spp <- function(code_vec) {
  spp_key[match(code_vec, spp_key$taxonID), 'scientificName']
}
