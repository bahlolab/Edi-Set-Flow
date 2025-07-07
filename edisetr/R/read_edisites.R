
edisites_spec <- readr::cols(
  ID               = col_character(),
  ADAR             = col_logical(),
  PASS_FREQ        = col_double(),
  REDI_ACC         = col_character(),
  REDI_REP_TYPE    = col_character(),
  REDI_REP_ID      = col_character(),
  REP_ID           = col_character(),
  Consequence      = col_character(),
  IMPACT           = col_character(),
  SYMBOL           = col_character(),
  Gene             = col_character(),
  Feature_type     = col_character(),
  Feature          = col_character(),
  BIOTYPE          = col_character(),
  EXON             = col_character(),
  INTRON           = col_character(),
  cDNA_position    = col_integer(),
  CDS_position     = col_integer(),
  Protein_position = col_integer(),
  Amino_acids      = col_character(),
  Codons           = col_character(),
  STRAND           = col_integer(),
  .default         = col_integer()
)

#' @export
read_edisites <-  function(
    sites_files,
    cluster         = NULL,
    min_med_dp      = 5L,
    min_med_vaf     = 0.001,
    max_med_vaf     = 0.999,
    adar_only       = TRUE,
    sample_group    = NULL,
    grp_min_med_dp  = 0,
    grp_min_med_vaf = 0,
    grp_max_med_vaf = 1,
    drop_failed     = FALSE,
    snp_detect      = TRUE,
    snp_detect_vaf  = c(0.1, 0.30),
    snp_detect_cor  = 0.85,
    snp_detect_n    = 3,
    snp_detect_p    = 0.70,
    snp_detect_dp   = 10,
    .n_max          = Inf
)
{

  if (is_scalar_character(sites_files) && str_detect(sites_files, '\\*')) {
    sites_files <- Sys.glob(sites_files)
  }

  if (is.null(cluster)) {

    sites_full <- read_tsv(
      sites_files,
      na = '.',
      col_types = edisites_spec,
      n_max = .n_max
    )

    counts <-
      sites_full %>%
      select(site_id = ID, starts_with(c('NALT_', 'NREF_'))) %>%
      pivot_longer(
        cols      = -site_id,
        names_to  = c(".value", "sample_id"),
        names_sep = "(?<=(NALT|NREF))_"
      ) %>%
      rename(n_alt = NALT, n_ref = NREF) %>%
      mutate(
        n_alt = replace_na(n_alt, 0),
        n_ref = replace_na(n_ref, 0),
        depth = n_alt + n_ref,
        vaf   = n_alt / depth
      ) %>%
      arrange(get_chrom(site_id), get_pos(site_id), site_id)


    smry <-
      counts %>%
      group_by(site_id) %>%
      summarise(
        med_depth = median(depth),
        med_vaf   = median(vaf, na.rm = TRUE),
        is_snp   = snp_detect(
          vaf, depth,
          snp_detect_vaf = snp_detect_vaf,
          snp_detect_cor = snp_detect_cor,
          snp_detect_n = snp_detect_n,
          snp_detect_p = snp_detect_p,
          snp_detect_dp = snp_detect_dp),
        .groups   = 'drop'
      ) %>%
      mutate(
        reason =
          str_c(
            if_else(is_snp & snp_detect    , 'snp'        , '.'),
            if_else(med_depth < min_med_dp , 'min_med_dp' , '.'),
            if_else(med_vaf   < min_med_vaf, 'min_med_vaf', '.'),
            if_else(med_vaf   > max_med_vaf, 'max_med_vaf', '.'),
            sep = ';'
          ) %>%
          str_remove_all(';?\\.') %>%
          str_remove('^;') %>%
          (\(x) replace(x, x=="", NA_character_)),
        pass = is.na(reason),
      ) %>%
      select(site_id, pass, reason, med_depth, med_vaf) %>%
      inner_join(
        sites_full  %>%
          select(-starts_with(c('NALT_', 'NREF_'))) %>%
          rename(site_id = ID, is_adar = ADAR) %>%
          janitor::clean_names() %>%
          mutate(
            is_adar = replace_na(is_adar, FALSE),
            rep_type = case_when(
              str_starts(rep_id, 'Alu') ~ 'ALU',
              !is.na(rep_id)            ~ 'REP',
              TRUE                      ~ 'NONREP'
            ),
            .after = rep_id
          ),
        by = join_by(site_id)
      ) %>%
      mutate(
        reason = case_when(
          !adar_only | is_adar ~ reason,
          is.na(reason)        ~ 'non_adar',
          TRUE                 ~ str_c(reason, ';non_adar')
        ),
        pass = is.na(reason)
      )

    # remove failing sites from counts
    counts_pass <-
      semi_join(counts, filter(smry,  pass), by = 'site_id')

    counts_fail <-
      semi_join(counts, filter(smry, !pass), by = 'site_id')

    if (!is.null(sample_group)) {
      grp_smry <-
        counts %>%
        left_join(sample_group, by = join_by(sample_id)) %>%
        group_by(site_id, group) %>%
        summarise(
          med_depth = median(depth),
          med_vaf   = median(vaf, na.rm = TRUE),
          .groups   = 'drop'
        ) %>%
        mutate(
          pass =
            med_depth >= grp_min_med_dp &
            med_vaf   >= grp_min_med_vaf &
            med_vaf   <= grp_max_med_vaf
        ) %>%
        # mutate(
        #   reason =
        #     str_c(
        #       if_else(med_depth < grp_min_med_dp , 'grp_min_med_dp' , '.'),
        #       if_else(med_vaf   < grp_min_med_vaf, 'grp_min_med_vaf', '.'),
        #       if_else(med_vaf   > grp_max_med_vaf, 'grp_max_med_vaf', '.'),
        #       sep = ';'
        #     ) %>%
        #     str_remove_all(';?\\.') %>%
        #     str_remove('^;') %>%
        #     (\(x) replace(x, x=="", NA_character_)),
        #   pass = is.na(reason),
        # ) %>%
        group_by(site_id) %>%
        mutate(n_pass = sum(pass)) %>%
        ungroup()

      counts_fail <-
        counts_pass %>%
        left_join(sample_group, by = join_by(sample_id)) %>%
        semi_join(
          grp_smry %>% filter(!pass | n_pass < 2),
          by = c('site_id', 'group')
        ) %>%
        select(-group) %>%
        bind_rows(counts_fail)

      counts_pass <-
        counts_pass %>%
        anti_join(
          counts_fail,
          by = c('site_id', 'sample_id')
        )

      smry <-
        smry %>%
        left_join(
          grp_smry %>%
            select(site_id, group, pass) %>%
            chop(group) %>%
            mutate(
              group = map_chr(group, ~ str_c(sort(.), collapse=';'))
            ) %>%
            mutate(pass = if_else(pass, 'groups_pass', 'groups_fail')) %>%
            pivot_wider(names_from = pass, values_from = group) %>%
            bind_rows(tibble(groups_pass = character(), groups_fail = character())) %>%
            select(groups_pass, groups_fail, everything()),
          by = 'site_id'
        )

    }
    sample_summary <-
      counts_pass %>%
      group_by(sample_id) %>%
      summarise(
        med_depth = median(depth),
        med_vaf   = median(vaf[depth >= min_med_dp], na.rm = TRUE)
      )

    ret <-list(
      summary = smry,
      sample_summary = sample_summary,
      counts  = counts_pass
    )
    if (!drop_failed) {
      ret$counts_fail = counts_fail
    }
    return(ret)


  } else {
    # multidply cluster workflow, just calls this function recursively
    cluster_check(cluster)
    multidplyr::cluster_library(cluster, "edisetr")
    multidplyr::cluster_assign_partition(cluster, sites_files = sites_files)

    # create a temporary id, blocking others with file
    tmp <- tempfile('')
    file.create(tmp)
    id <- basename(tmp)
    multidplyr::cluster_assign(cluster, id = id)

    multidplyr::cluster_send(
      cluster,
      {
        res <-
          read_edisites(
            sites_files,
            min_med_dp      = !! min_med_dp,
            min_med_vaf     = !! min_med_vaf,
            max_med_vaf     = !! max_med_vaf,
            sample_group    = !! sample_group,
            grp_min_med_dp  = !! grp_min_med_dp,
            grp_min_med_vaf = !! grp_min_med_vaf,
            grp_max_med_vaf = !! grp_max_med_vaf,
            drop_failed     = !! drop_failed,
            snp_detect      = !! snp_detect,
            snp_detect_vaf  = !! snp_detect_vaf,
            snp_detect_cor  = !! snp_detect_cor,
            snp_detect_n    = !! snp_detect_n,
            snp_detect_p    = !! snp_detect_p,
            snp_detect_dp   = !! snp_detect_dp,
            .n_max          = !! .n_max
          )

        assign(str_c('counts_', id), res$counts)
        summary <- res$summary
        sample_summary <- res$sample_summary
        if (!(!!drop_failed)) {
          counts_fail <- res$counts_fail
        }

        rm(res)
        gc()
      }
    )

    ret <- list(
      summary =
        multidplyr::party_df(cluster, 'summary') %>%
        collect(),
      sample_summary =
        multidplyr::party_df(cluster, 'sample_summary') %>%
        collect() %>%
        group_by(sample_id) %>%
        summarise(
          med_depth = median(round(med_depth)),
          med_vaf   = median(med_vaf)
        ),
      counts  = multidplyr::party_df(cluster, str_c('counts_', id))
    )
    if (!drop_failed) {
      ret$counts_fail = collect(multidplyr::party_df(cluster, 'counts_fail'))
    }
    return(ret)

  }
}

snp_detect <- function(
    vaf,
    depth,
    snp_detect_vaf,
    snp_detect_cor,
    snp_detect_n,
    snp_detect_p,
    snp_detect_dp
)
{
  vaf <- vaf[depth > snp_detect_dp]
  gt  <- rep(NA_integer_, length(vaf))
  gt[vaf <= snp_detect_vaf[1]  ] <- 0L
  gt[vaf >= 1-snp_detect_vaf[1]] <- 2L
  gt[vaf >= snp_detect_vaf[2] &
     vaf <= 1-snp_detect_vaf[2]] <- 1L


  which_gt <- !is.na(gt)
  if (sum(which_gt) / length(gt) > snp_detect_p) {
    cnt <- table(gt)
    if (length(cnt) == 3 & min(cnt) > snp_detect_n) {
      if (cor(gt[which_gt], vaf[which_gt]) > snp_detect_cor) {
        return(TRUE)
      }
    }
  }
  return(FALSE)
}


