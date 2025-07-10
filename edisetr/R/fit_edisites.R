
#' @export
fit_edisites <- function(
    counts_df,
    covariates_df,
    fixed_effects = character(),
    model = c('quasibinomial', 'binomial', 'linear', 'arcsine'),
    fdr_bh     = TRUE,
    fdr_storey = FALSE,
    fdr_ash    = FALSE,
    fdr_emp    = FALSE,
    permute_group = FALSE
)
{

  model     <- match.arg(model)

  stopifnot(
    is.data.frame(counts_df) || inherits(counts_df, 'multidplyr_party_df'),
    is.data.frame(covariates_df),
    all(c("site_id", "sample_id", "n_alt", "n_ref") %in% colnames(head(counts_df))),
    all(fixed_effects %in% colnames(covariates_df)),
    is.character(model) && length(model) > 0,
    rlang::is_bool(fdr_bh),
    rlang::is_bool(fdr_storey),
    rlang::is_bool(fdr_ash),
    rlang::is_bool(fdr_emp),
    rlang::is_bool(permute_group)
  )

  fixed_effects <- unique(c('group', fixed_effects))

  if (model == 'binomial') {
    family  <- binomial(link = "logit")
    formula <- reformulate(
      response = "cbind(n_alt, n_ref)",
      termlabels = fixed_effects
    )
  } else if (model == 'quasibinomial') {
    family  <- quasibinomial(link = "logit")
    formula <- reformulate(
      response = "cbind(n_alt, n_ref)",
      termlabels = fixed_effects
    )
  } else if (model == 'linear') {
    family  <- gaussian(link = 'identity')
    formula <- reformulate(
      response = "n_alt / (n_alt + n_ref)",
      termlabels = fixed_effects
    )
  } else if (model == 'arcsine') {
    family  <- gaussian(link = 'identity')
    formula <- reformulate(
      response = "asin(sqrt(n_alt / (n_alt + n_ref)))",
      termlabels = fixed_effects
    )
  }
  if (inherits(counts_df, 'multidplyr_party_df')) {
    cluster_assign(
      counts_df$cluster,
      formula = formula,
      family  =family
    )
  }

  results <-
    counts_df %>%
    filter(depth > 0) %>%
    inner_join(covariates_df, by = 'sample_id', copy = TRUE) %>%
    group_by(site_id) %>%
    filter(n_distinct(group) > 1) %>%
    (function(x) {
      if (permute_group) {
        mutate(x, group = sample(group))
      } else {
        x
      }
    }) %>%
    summarise(
      result = edisetr:::fit_glm(
        data    = pick(n_alt, n_ref, all_of(!!fixed_effects)),
        formula = formula,
        family  = family,
        .with_null = fdr_emp
      ),
    ) %>%
    ungroup() %>%
    collect()

  get_table <- function(table) {
    results %>%
      mutate(x = result[[table]]) %>%
      select(site_id, x) %>%
      unnest(x)
  }

  summary   <- get_table('summary')
  margins   <- get_table('margins')

  if (model == 'arcsine') {
    margins <- inverse_arcsine_margins(margins)
  }

  contrasts <-
    get_table('contrasts') %>%
    group_by(contrast) %>%
    calc_fdr(bh = fdr_bh, storey = fdr_storey, ash = fdr_ash, emp = fdr_emp) %>%
    ungroup()

  anova <-
    get_table('anova') %>%
    filter(term != '<none>') %>%
    group_by(term) %>%
    calc_fdr(bh = fdr_bh, storey = fdr_storey, ash = fdr_ash, emp = fdr_emp) %>%
    ungroup()

  rm(results); gc()

  if (inherits(counts_df, 'multidplyr_party_df')) {
    invisible(cluster_call(counts_df$cluster, {gc()}))
  }

  return(
    list(
      summary   = summary,
      margins   = margins,
      contrasts = contrasts,
      anova     = anova
    )
  )
}

make_clean_names <- memoise::memoise(janitor::make_clean_names)

#' @importFrom broom tidy
#' @importFrom janitor clean_names
fit_glm <- function(data, formula, family, .with_null = FALSE, .do_null = FALSE) {

  # make sure no invariant factors present in formula
  for (v in all.vars(formula[[3]])) {
    col <- data[[v]]
    if (!is.numeric(col)) {
      if (length(unique(na.omit(col))) < 2) {
        formula <- update(formula, paste(". ~ . -", v))
      }
    }
  }

  if (.do_null) {
    data$group <- sample(data$group)
  }

  fit <- glm(formula = formula, data = data, family = family)

  smry <- tidy(fit) %>% rename_with(make_clean_names)

  contrasts <-
    pairwise_glm_contrasts(fit, group_var = 'group') %>%
    mutate(dispersion = sum(residuals(fit, type = "pearson")^2) / df.residual(fit))

  margins <- group_emm(fit, data, group='group')

  if (family$family == 'binomial') {
    anova <-
      drop1(fit, ~ ., test="Chisq") %>%
      as.data.frame() %>%
      rownames_to_column('term') %>%
      rename_with(make_clean_names) %>%
      rename(p_value = pr_chi)
  } else {
    anova <-
      drop1(fit, ~ ., test="F") %>%
      as.data.frame() %>%
      rownames_to_column('term') %>%
      rename_with(make_clean_names) %>%
      rename(p_value = pr_f)
  }

  if (.with_null) {
    null_fit <- fit_glm(data, formula, family, .with_null = FALSE, .do_null = TRUE)
    contrasts$null_p_value <- null_fit$contrasts[[1]]$p_value
    anova$null_p_value     <- null_fit$anova[[1]]$p_value
  }

  result <-
    tibble(
      summary   = list(smry),
      contrasts = list(contrasts),
      margins   = list(margins),
      anova     = list(anova)
    )

  return(result)
}



calc_quasi <- function(X, name = 'quasibinomial', return_all = TRUE, clamp = TRUE) {

  Y <-
    filter(X, model == 'binomial') %>%
    mutate(
      model = name,
      dispersion_ = `if`(clamp, pmax(dispersion, 1), dispersion),
      std_error   = sqrt(dispersion_) * std_error,
      statistic   = estimate / std_error,
      p_value     = 2 * exp(pt(-abs(statistic), df = df_residual, log.p = TRUE))
    ) %>%
    select(-dispersion_)

  if (return_all) {
    return(arrange_all(bind_rows(X, Y)))
  }
  return(Y)
}

calc_fdr <- function(data, bh=FALSE, storey=FALSE, ash=FALSE, emp=FALSE) {
  if (bh) {
    data <-
      mutate(
        data,
        q_value_bh = p.adjust(p_value, method = 'BH')
      )
  }
  if (storey) {
    data <-
      mutate(
        data,
        q_value_storey = qvalue::qvalue(p_value)[['qvalues']]
      )
  }
  if (ash) {
    data <-
      mutate(
        data,
        q_value_ash   = ashr::ash(estimate, std_error, df = floor(median(df_residual)))$result$qvalue,
      )
  }
  if (emp) {
    data <-
      mutate(
        data,
        q_value_emp   = empirical_qvalue(p_obs = p_value, p_null = null_p_value)
      )
  }
  return(data)
}

pairwise_glm_contrasts <- function(fit, group_var) {
  # Get the levels and numeric index pairs
  levels_group <- fit$xlevels[[group_var]]
  nlev <- length(levels_group)
  comb <- utils::combn(nlev, 2)                     # 2×M matrix of row-indices
  i1  <- comb[1,]
  i2  <- comb[2,]

  #  Build the 'newdata' grid (one row per level)
  #  and force other factors to carry all their levels
  mf <- model.frame(fit)
  newdata <- as.data.frame(matrix(nrow = nlev, ncol = 0))
  newdata[[group_var]] <- factor(levels_group, levels = levels_group)

  # other factor covariates at their reference
  other_factors <- setdiff(names(fit$xlevels), group_var)
  for (v in other_factors) {
    newdata[[v]] <- factor(fit$xlevels[[v]][1],
                           levels = fit$xlevels[[v]])
  }
  # numeric covariates at their mean
  numeric_covs <- setdiff(names(mf),
                          c(group_var, names(fit$xlevels)))
  for (v in numeric_covs) {
    newdata[[v]] <- mean(mf[[v]], na.rm = TRUE)
  }

  # Get design matrix, coefs, vcov
  X        <- model.matrix(delete.response(terms(fit)), newdata)
  beta     <- coef(fit)
  V        <- vcov(fit)
  df_resid <- df.residual(fit)
  fam      <- family(fit)$family

  # Build one big contrast matrix (M × p)
  Cmat     <- X[i1, , drop = FALSE] - X[i2, , drop = FALSE]

  # Vectorized estimates and SEs
  est      <- as.vector(Cmat %*% beta)           # M estimates
  Vc       <- Cmat %*% V                         # M×p  matrix
  var      <- rowSums(Vc * Cmat)                 # M variances
  se       <- sqrt(var)                          # M SEs
  stat     <- est / se                           # Wald zs or ts

  # p-values on the correct scale
  if (fam %in% c("quasibinomial", "quasipoisson")) {
    pval <- 2 * exp(pt(-abs(stat), df = df_resid, log.p = TRUE))
  } else {
    pval <- 2 * exp(pnorm(-abs(stat), log.p = TRUE))
  }

  # Assemble result in one data.frame
  tibble(
    contrast   = paste(levels_group[i1], "-", levels_group[i2]),
    estimate   = est,
    std_error  = se,
    statistic  = stat,
    df_residual= df_resid,
    p_value    = pval
  )
}

group_emm <- function(fit, data, group="group", alpha = 0.05) {
  # Identify model predictors
  mf    <- model.frame(fit, data = data)
  preds <- names(mf)[-1]
  if (!(group %in% preds)) {
    stop("`group` must be one of the model predictors.")
  }

  # Build reference grid
  grids <- lapply(preds, function(v) {
    col <- data[[v]]
    if (is.numeric(col)) {
      mean(col, na.rm = TRUE)
    } else if (is.character(col)) {
      unique(na.omit(col))
    } else if (is.factor(col)) {
      levels(col)
    } else {
      stop("Don't know what to do with ", class(col), " column")
    }
  })
  names(grids) <- preds
  refg <- do.call(expand.grid, c(grids, list(stringsAsFactors = FALSE)))

  # Count observations per cell
  facs     <- preds[!sapply(data[preds], is.numeric)]
  cts      <- table(do.call(paste, c(data[facs], sep = "\r")))
  key_refg <- do.call(paste, c(refg[facs], sep = "\r"))
  refg$count <- as.integer(cts)[ match(key_refg, names(cts)) ]
  refg$count[is.na(refg$count)] <- 0L

  # Normalize weights within each group level
  refg$weight <- ave(refg$count, refg[[group]], FUN = function(x) {
    if (sum(x) == 0) rep(0, length(x)) else x / sum(x)
  })

  # Design matrix, coefs, vcov
  X    <- model.matrix(delete.response(terms(fit)), data = refg)
  bhat <- coef(fit)
  V    <- vcov(fit)

  # Compute EMM + CI for each level
  levels_g <- unique(refg[[group]])
  out <- lapply(levels_g, function(g) {
    sel      <- refg[[group]] == g
    w        <- refg$weight[sel]
    Xg       <- X[sel, , drop = FALSE]
    Lbar     <- colSums(Xg * w)
    eta      <- sum(Lbar * bhat)
    se_link  <- sqrt(as.numeric(Lbar %*% V %*% Lbar))
    df       <- df.residual(fit)
    crit     <- if (inherits(fit, "glm")) qnorm(1 - alpha/2) else qt(1 - alpha/2, df)
    ci_link  <- eta + c(-1, 1) * crit * se_link
    invlink  <- fit$family$linkinv
    mu_eta   <- fit$family$mu.eta(eta)
    se_resp  <- se_link * mu_eta

    data.frame(
      group     = g,
      estimate  = invlink(eta),
      std_error = se_resp,
      ci_lower  = invlink(ci_link[1]),
      ci_upper  = invlink(ci_link[2]),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, out)
}

inverse_arcsine_margins <- function(df) {
  # Preserve original transformed values
  orig_est   <- df$estimate
  orig_se    <- df$std_error
  orig_lower <- df$ci_lower
  orig_upper <- df$ci_upper

  # Overwrite with back‐transformed rates
  df$estimate   <- sin(orig_est)^2
  df$std_error  <- orig_se * sin(2 * orig_est)
  df$ci_lower   <- sin(orig_lower)^2
  df$ci_upper   <- sin(orig_upper)^2

  df
}

empirical_qvalue <- function(p_obs, p_null) {
  # sort observed
  o_obs    <- order(p_obs)
  p_sorted <- p_obs[o_obs]
  n_obs    <- length(p_sorted)

  # sort nulls
  p_null_sorted <- sort(p_null)
  n_null        <- length(p_null_sorted)

  # raw empirical q
  k     <- findInterval(p_sorted, p_null_sorted, rightmost.closed = TRUE)
  q_raw <- (k / n_null) / (seq_len(n_obs) / n_obs)

  # BH-style monotonic smoothing
  q_smooth <- rev(cummin(rev(q_raw)))

  # floor the *smoothed* q at the observed p
  q_bounded <- pmax(q_smooth, p_sorted)

  # optional: re-smooth to ensure monotonicity is preserved
  q_final <- rev(cummin(rev(q_bounded)))

  # cap at 1 and restore original order
  q_final <- pmin(1, q_final)
  out     <- numeric(n_obs)
  out[o_obs] <- q_final

  out
}