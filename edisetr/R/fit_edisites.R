
#' Fit per-site differential RNA editing models
#'
#' Fits a generalised linear model per editing site to test for differences in
#' editing rate across experimental groups. Each site is modelled independently
#' using `group` plus any `fixed_effects` as predictors, and the response is
#' derived from the per-sample alternate/reference allele counts. When
#' `counts_df` is a \pkg{multidplyr} party data frame the per-site fits are
#' distributed across the cluster's workers.
#'
#' @param counts_df A data frame (or \pkg{multidplyr} `party_df`) of per-sample
#'   allele counts, one row per site/sample, with at least the columns
#'   `site_id`, `sample_id`, `n_alt` (alternate/edited allele count) and
#'   `n_ref` (reference allele count). A `depth` column (`n_alt + n_ref`) is
#'   used to drop zero-depth observations. An optional `weights` column is
#'   passed through as GLM observation weights (ignored for the `*-ebbr` models,
#'   which derive their own weights).
#' @param covariates_df A data frame of per-sample covariates, one row per
#'   `sample_id`. Must contain `sample_id`, a `group` column, and any columns
#'   named in `fixed_effects`. Joined to `counts_df` by `sample_id`.
#' @param fixed_effects Character vector of additional covariate columns (from
#'   `covariates_df`) to include as fixed effects alongside `group`. Defaults to
#'   none.
#' @param model The model family / response transform, one of:
#'   \describe{
#'     \item{`quasibinomial`}{Quasibinomial logit GLM on `cbind(n_alt, n_ref)`
#'       (accounts for overdispersion).}
#'     \item{`binomial`}{Binomial logit GLM on `cbind(n_alt, n_ref)`.}
#'     \item{`linear`}{Gaussian (identity) model on the raw editing proportion
#'       `n_alt / (n_alt + n_ref)`.}
#'     \item{`arcsine`}{Gaussian model on the arcsine-square-root transformed
#'       proportion; margins are back-transformed to the rate scale.}
#'     \item{`linear-ebbr`}{As `linear`, but the response is the per-site
#'       empirical-Bayes Beta-binomial posterior mean (see
#'       `inverse_variance_weighting`). Requires the optional \pkg{ebbr}
#'       package.}
#'     \item{`arcsine-ebbr`}{As `arcsine`, but on the empirical-Bayes posterior
#'       mean. Requires the optional \pkg{ebbr} package.}
#'   }
#' @param inverse_variance_weighting Logical; for the `*-ebbr` models, weight
#'   each observation by the inverse of its Beta-binomial posterior variance so
#'   noisy low-depth samples contribute less. Default `TRUE`. Has no effect for
#'   the non-`ebbr` models.
#' @param weight_cap_q Numeric quantile in (0, 1) used to cap the per-site
#'   inverse-variance weights at `median(w) + mad(w) * qnorm(weight_cap_q)`,
#'   preventing a single very-high-depth sample from dominating a site's fit.
#'   Set to `NULL` or `NA` to disable capping. Default `0.95`. Only relevant
#'   when `inverse_variance_weighting = TRUE`.
#' @param fdr_bh Logical; add Benjamini-Hochberg `q_value_bh` to the contrasts
#'   and ANOVA tables. Default `TRUE`.
#' @param fdr_storey Logical; add Storey `q_value_storey` (requires
#'   \pkg{qvalue}). Default `FALSE`.
#' @param fdr_ash Logical; add adaptive-shrinkage `q_value_ash` (requires
#'   \pkg{ashr}). Default `FALSE`.
#' @param fdr_emp Logical; add an empirical `q_value_emp` computed against a
#'   null obtained by permuting `group` within each site. Doubles the fitting
#'   work. Default `FALSE`.
#' @param permute_group Logical; permute the `group` labels within each site
#'   before fitting, yielding a global null distribution (for diagnostics).
#'   Default `FALSE`.
#'
#' @return A named list of tidy data frames, each keyed by `site_id`:
#'   \describe{
#'     \item{`summary`}{Per-term model coefficient estimates.}
#'     \item{`margins`}{Estimated marginal means per `group` level (on the
#'       response scale; back-transformed for the arcsine models), with
#'       confidence intervals.}
#'     \item{`contrasts`}{Pairwise between-group contrasts with standard errors,
#'       test statistics, p-values and requested FDR columns.}
#'     \item{`anova`}{Per-term drop-one tests with p-values and requested FDR
#'       columns.}
#'   }
#'
#' @seealso [read_edisites()] for producing `counts_df`.
#' @export
fit_edisites <- function(
    counts_df,
    covariates_df,
    fixed_effects = character(),
    model = c('quasibinomial', 'binomial', 'linear', 'arcsine',
              'linear-ebbr', 'arcsine-ebbr'),
    inverse_variance_weighting = TRUE,
    weight_cap_q               = 0.95,
    fdr_bh     = TRUE,
    fdr_storey = FALSE,
    fdr_ash    = FALSE,
    fdr_emp    = FALSE,
    permute_group = FALSE
)
{

  model     <- match.arg(model)

  # The '*-ebbr' models rely on the optional 'ebbr' package (GitHub-only), so
  # fail early with install instructions rather than deep inside a cluster call.
  if (model %in% c('linear-ebbr', 'arcsine-ebbr') &&
      !requireNamespace('ebbr', quietly = TRUE)) {
    stop("model = '", model, "' requires the 'ebbr' package, which is not installed.\n",
         "Install it with: remotes::install_github('dgrtwo/ebbr')", call. = FALSE)
  }

  stopifnot(
    is.data.frame(counts_df) || inherits(counts_df, 'multidplyr_party_df')
  )
  counts_cols <- colnames(head(counts_df))

  stopifnot(
    is.data.frame(covariates_df),
    all(c("site_id", "sample_id", "n_alt", "n_ref") %in% counts_cols),
    all(fixed_effects %in% colnames(covariates_df)),
    is.character(model) && length(model) > 0,
    rlang::is_bool(inverse_variance_weighting),
    rlang::is_bool(fdr_bh),
    rlang::is_bool(fdr_storey),
    rlang::is_bool(fdr_ash),
    rlang::is_bool(fdr_emp),
    rlang::is_bool(permute_group)
  )

  if ('weights' %in% counts_cols) {
    message('fitting GLM using weights, remove weights column to disable')
  }

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
  } else if (model == 'linear-ebbr') {
    # linear model on the empirical-Bayes posterior mean (see eb_shrink_site)
    family  <- gaussian(link = 'identity')
    formula <- reformulate(
      response = "post_mean",
      termlabels = fixed_effects
    )
  } else if (model == 'arcsine-ebbr') {
    # arcsine model on the empirical-Bayes posterior mean (see eb_shrink_site)
    family  <- gaussian(link = 'identity')
    formula <- reformulate(
      response = "asin(sqrt(post_mean))",
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
        data    = pick(n_alt, n_ref, all_of(!!fixed_effects), any_of('weights')),
        formula = formula,
        family  = family,
        .with_null    = !!fdr_emp,
        .ebbr         = !!(model %in% c('linear-ebbr', 'arcsine-ebbr')),
        .iv_weighting = !!inverse_variance_weighting,
        .cap_q        = !!weight_cap_q
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

  if (model %in% c('arcsine', 'arcsine-ebbr')) {
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
fit_glm <- function(data, formula, family, .with_null = FALSE, .do_null = FALSE,
                    .ebbr = FALSE, .iv_weighting = TRUE, .cap_q = 0.95) {

  # empirical-Bayes shrinkage for the '*-ebbr' models: replace raw per-sample
  # rates with the Beta-binomial posterior mean and (optionally) precision
  # weights, both derived per site. Done before the invariant-factor check so
  # the post_mean response column exists when glm() evaluates the formula.
  if (.ebbr) {
    eb <- eb_shrink_site(data$n_alt, data$n_ref, .iv_weighting, .cap_q)
    data$post_mean <- eb$post_mean
    data$weights   <- eb$weights
  }

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
  if (! 'weights' %in% colnames(data)) {
    data$weights <- 1
  }

  fit <- glm(formula = formula, data = data, family = family,
             weights = weights)

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

# Empirical-Bayes shrinkage of a single site's per-sample editing rates.
# Fits a Beta(alpha0, beta0) prior across this site's samples (borrowing strength
# across samples measured at the same site) via ebbr, then returns each sample's
# Beta-binomial posterior:
#   alpha_i = n_alt_i + alpha0,  beta_i = n_ref_i + beta0
#   post_mean_i = alpha_i / (alpha_i + beta_i)
#   post_var_i  = alpha_i * beta_i / ((alpha_i + beta_i)^2 * (alpha_i + beta_i + 1))
# Precision weights (1 / post_var) shrink noisy low-depth samples; they are
# optionally capped per site (at median + MAD * qnorm(cap_q)) so one very-high-
# depth sample cannot dominate the fit. cap_q = NULL/NA disables the cap.
eb_shrink_site <- function(n_alt, n_ref, iv_weighting = TRUE, cap_q = 0.95) {
  n_tot <- n_alt + n_ref
  prior <- tryCatch(
    ebbr::ebb_fit_prior(tibble(x = n_alt[n_tot > 0], n = n_tot[n_tot > 0]), x, n),
    error = function(e) NULL
  )
  if (is.null(prior)) {
    # degenerate site (prior did not converge): fall back to raw rate, no weighting
    post_mean <- ifelse(n_tot > 0, n_alt / n_tot, 0)
    return(tibble(post_mean = post_mean, weights = rep(1, length(n_alt))))
  }
  a0    <- prior$parameters$alpha
  b0    <- prior$parameters$beta
  alpha <- n_alt + a0
  beta  <- n_ref + b0
  post_mean <- alpha / (alpha + beta)
  post_var  <- (alpha * beta) / ((alpha + beta)^2 * (alpha + beta + 1))
  if (iv_weighting) {
    w <- 1 / post_var
    if (!is.null(cap_q) && !is.na(cap_q)) {
      w <- pmin(w, median(w) + mad(w) * qnorm(cap_q))
    }
  } else {
    w <- rep(1, length(n_alt))
  }
  tibble(post_mean = post_mean, weights = w)
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
  preds <- names(mf)[-1] %>% setdiff("(weights)")
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
