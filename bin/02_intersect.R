#!/usr/bin/env Rscript

require(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

name <- args[1]
siteStats_rediJOIN_fn <- args[2]
boundingENSG_fn <- args[3]
RM_Repeats_intersect_fn <- args[4]
genFeatures_intersect_fn <- args[5]
all_siteStats_fn <- args[6]
sample_thresh <- as.integer(args[7])
replace_chr <- as.logical(args[8])


# res_other <- readRDS(res_other_fn)
siteStats_rediJOIN <- readRDS(siteStats_rediJOIN_fn)
boundingENSG <- readRDS(boundingENSG_fn)
RM_Repeats_intersect <- readRDS(RM_Repeats_intersect_fn)
genFeatures_intersect <- readRDS(genFeatures_intersect_fn)
all_siteStats <- readRDS(all_siteStats_fn)

# Tag cognate strand ------------------------------------------------------

tag_cognate_p1 <-
  siteStats_rediJOIN %>%
  mutate(
    basechange = str_replace_all(basechange, "->", "_")
  ) %>%
  select(1:basechange, type) %>%
  rename(REDI = type) %>%
  left_join(
    boundingENSG,
    by = "siteID"
  ) %>%
  group_by(siteID) %>%
  mutate(n_oLap = n()) %>%
  ungroup() %>%
  # filter(is.na(strand)) %>% select(siteID) %>% distinct() %>% dim()
  select(
    1:10,
    n_oLap,
    everything()
  ) %>%
  mutate(strandMatch = case_when(
    is.na(strand) ~ NA_character_,
    basechange == "A_G" & n_oLap == 1 & strand == "+" ~ "single_cogn_Watson",
    basechange == "A_G" & n_oLap == 1 & strand == "-" ~ "single_cogn_Crick",
    TRUE ~ "other"
  )) %>%
  mutate(REDI = ifelse(is.na(REDI), "uncatalog", REDI))

tag_cognate_p2 <-
  tag_cognate_p1 %>%
  filter(strandMatch == "other" | is.na(strandMatch)) %>%
  group_by(siteID) %>%
  mutate(site_status = case_when(
    REDI == "uncatalog" & is.na(strandMatch) ~ "drop",
    REDI != "uncatalog" & strandMatch == "other" & basechange == "A_G" & all(strand == "+") ~ "multi_cogn_Watson",
    REDI != "uncatalog" & strandMatch == "other" & basechange == "A_G" & all(strand == "-") ~ "multi_cogn_Crick",
    REDI != "uncatalog" & strandMatch == "other" & basechange == "A_G" & any(strand == "+") ~ "multi_cogn_mix",
    # REDI!='uncatalog' & strandMatch=='other' & basechange=="T_C" &  any(strand=="-")  ~ 'multi_cogn_mix',
    REDI != "uncatalog" & strandMatch == "other" & basechange == "A_G" & !any(strand == "+") ~ "multi_noncogn",
    # REDI!='uncatalog' & strandMatch=='other' & basechange=="T_C" & !any(strand=="-")  ~ 'multi_noncogn',
    is.na(strandMatch) ~ "no_ENSG",
    TRUE ~ "oLap_mix"
  )) %>%
  ungroup()

sites_tagged <-
  rbind(
    tag_cognate_p1 %>% filter(!is.na(strandMatch), strandMatch != "other") %>%
      mutate(site_status = strandMatch),
    tag_cognate_p2
  )

sites_tagged_context <-
  sites_tagged %>%
  # select(siteID, dbSNP_status, REDI, basechange,strand, strandMatch,site_status,contains('gene')) %>% distinct() %>%
  select(siteID, REDI, basechange, strand, strandMatch, site_status, contains("gene")) %>%
  distinct() %>%
  mutate(site_context = case_when(
    str_detect(site_status, "single") ~ "cognate",
    site_status %in% c("multi_cogn_Watson", "multi_cogn_Crick") ~ "cognate",
    site_status %in% c("oLap_mix", "multi_cogn_mix", "no_ENSG") ~ "ambiguous", # NB 'no_ENSG' denotes sites in RADAR db, but without ENSG bounding genes
    site_status %in% c("multi_noncogn", "non_cognate") ~ "non_cognate",
    site_status %in% c("drop") ~ "drop"
  )) %>%
  filter(site_context != "drop")

saveRDS(
  sites_tagged_context,
  str_c(name, ".sites_tagged_context.rds")
)

# Drop non-cognate sites --------------------------------------------------

siteStats_anno <-
  sites_tagged_context %>%
  filter(
    site_status != "drop",
    site_status != "non_cognate",
    site_status != "multi_noncogn"
  ) %>%
  mutate(tag = ifelse(REDI == "uncatalog" & site_status == "no_ENSG", "drop", "retain")) %>%
  filter(tag != "drop") %>%
  rename(
    ENSGID = gene_id,
    ENSGbioType = gene_biotype
  )


siteStats_anno_RM <-
  siteStats_anno %>%
  left_join(
    select(RM_Repeats_intersect, 1, Repeat_id, RM_strand, RM_repType),
    by = "siteID"
  )

# Report genic features ---------------------------------------------------

genFeatures_spread <-
  genFeatures_intersect %>%
  semi_join(siteStats_anno_RM, by = "siteID") %>%
  mutate(dummy = 1) %>%
  spread(feature, dummy) %>%
  bind_rows(
    tibble(
      exon = numeric(),
      five_prime_utr = numeric(),
      gene = numeric(),
      start_codon = numeric(),
      stop_codon = numeric(),
      three_prime_utr = numeric()
    )
  ) %>% 
  mutate(
    intronic = ifelse(gene == 1 & is.na(exon) & is.na(five_prime_utr) & is.na(stop_codon) & is.na(three_prime_utr), 1, NA),
    exonic   = ifelse(is.na(intronic), 1, NA)
  ) %>%
  select(-exon) %>%
  select(siteID, five_prime_utr:exonic, -gene) %>%
  distinct() %>% # count(siteID) #
  gather(key, value, five_prime_utr:exonic) %>%
  na.omit() %>%
  distinct() %>%
  pivot_wider(
    names_from = key,
    values_from = value
  )

genFeatures_Nest <-
  genFeatures_spread %>%
  gather(key, value, -siteID) %>%
  na.omit() %>%
  select(-value) %>%
  # nest(-siteID, key) #will crash
  group_by(siteID) %>%
  nest("genicFeature" = key) %>%
  ungroup()

genFeatures_deDup <-
  genFeatures_spread %>%
  gather(key, value, -c(siteID)) %>%
  na.omit() %>%
  select(-value) %>%
  group_by(siteID) %>%
  mutate(n_Obs = n()) %>%
  ungroup() %>%
  # Drops 'exonic' counts that are redundant (i.e. with 3'utr / 5'utr)
  mutate(tag3 = ifelse(key == "exonic" & n_Obs > 1, "drop", "retain")) %>%
  # group_by(siteID) %>% filter(any(tag3=='drop')) %>% arrange(siteID)
  # NB a site may be intronic and exonic, but only the latter is also redundant w 3' utr etc.
  filter(tag3 == "retain") %>%
  group_by(siteID) %>%
  mutate(n_Obs = n()) %>%
  ungroup() %>%
  mutate(genicFeature_summ = ifelse(n_Obs > 1, paste0("multiple:_", n_Obs), key)) %>%
  # filter(genicFeature_summ=='multiple:_3') %>% arrange(siteID)
  select(siteID, genicFeature_summ) %>%
  left_join(genFeatures_Nest, by = "siteID")


# Finalize site annotation ------------------------------------------------

# 1. retain high-confidence RNA editing sites
# 2. join siteStats_anno_RM with genic features
siteStats_filt <-
  siteStats_anno_RM %>% # redundant on multi-cognate sites
  mutate(tag2 = case_when( # dbSNP_status=='commonSNP' ~ 'drop',
    site_context == "cognate" ~ "retain",
    site_context == "ambiguous" & !is.na(RM_repType) ~ "retain",
    site_context == "ambiguous" & REDI == "uncatalogue" & is.na(RM_repType) ~ "drop",
    TRUE ~ "drop"
  )) %>%
  filter(tag2 == "retain") %>%
  left_join(
    genFeatures_deDup %>%
      select(1, 2, genicFeature) %>%
      distinct(),
    by = c("siteID")
  ) %>%
  # left_join(anno %>% select(gene_id, description), by=c('ENSGID'='gene_id')) %>% #NR
  left_join(
    all_siteStats %>%
      mutate(basechange = str_replace_all(basechange, "->", "_")),
    by = c("siteID", "basechange")
  ) %>%
  mutate(RM = ifelse(RM_repType == "Alu", "ALU", "REP")) %>% # NAs are preserved
  mutate(
    REDI_RM = case_when(
      REDI != "uncatalog" ~ REDI, # prefer REDI annotation
      !is.na(RM) ~ paste0(RM, "_novel"),
      is.na(RM) ~ "NONREP_novel"
    )
  )

saveRDS(
  siteStats_filt,
  str_c(name, ".siteStats_filt.rds")
)

filt_forSamDepth <-
  siteStats_filt %>%
  filter(nSamples >= sample_thresh) %>%
  # filter siteID by prevalence. single donor; 3 tissues (subset from 5)
  select(siteID) %>%
  distinct() %>%
  separate(siteID, into = c("chr", "pos"), convert = TRUE, sep = "_") %>%
  mutate(pos = as.numeric(pos), start = pos - 1) %>%
  select(chr, start, pos) %>% 
  mutate(chr = `if`(replace_chr, str_c('chr', chr), chr)) %>% 
  arrange_all()

write_tsv(
  filt_forSamDepth,
  col_names = FALSE,
  file = str_c(name, ".sites_filt.bed.gz")
)
