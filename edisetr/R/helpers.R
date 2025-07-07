
#' @export
get_base_change <- function(site_id) {
  stopifnot(is.character(site_id))

  str_extract(site_id, '[ACGT]_[ACGT]:.') %>%
    str_replace('_', '>') %>%
    str_replace(':', '(') %>%
    str_c(')')
}

#' @export
get_chrom <- function(site_id) {
  stopifnot(is.character(site_id))
  str_extract(site_id, '^[^_]+')
}

#' @export
get_pos <- function(site_id) {
  stopifnot(is.character(site_id))
  str_extract(site_id, '(?<=^[^_]{1,50}_)[0-9]+') %>% as.integer()
}

#' @export
cluster_check <- function(cluster, stop = TRUE) {
  if (! "multidplyr_cluster" %in% class(cluster)) {
    if (stop) {
      stop("Cluster is not a multidplyr cluster")
    } else {
      return(FALSE)
    }
  }
  if (!all(purrr::map_lgl(cluster, ~ .$is_alive()))) {
    if (stop) {
      stop("Cluster workers are not all alive")
    } else {
      return(FALSE)
    }
  }
  return(TRUE)
}
