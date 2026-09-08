################################################################################
## This file contains functions to perform analyses ############################
## Author: Malou Magnani and Carolien C.H.M. Maas ##############################
################################################################################

## ============================================================================
## KFRE prediction equations
## ============================================================================

## Compute the linear predictor (PI), which is the same for all time horizons.
## Coefficients and centering constants are the published 4-variable KFRE
## coefficients (age, sex, eGFR, urine albumin-to-creatinine ratio); the
## constants being subtracted (7.036, 0.5642, 7.222, 5.137) are the means
## used to center each predictor in the original derivation cohort.
PI_KFRE <- function(age, male, egfr, alb) {
  # age term: age (in decades) centered on its cohort mean
  PI <- -0.2201 *
    (age / 10 - 7.036) +
    # sex term: male is a 0/1 indicator, centered on cohort mean proportion male
    0.2467 * (male - 0.5642) -
    # eGFR term: eGFR (in 5 mL/min/1.73m^2 units) centered on its cohort mean
    0.5567 * (egfr / 5 - 7.222) +
    # albuminuria term: log(albumin) centered on its cohort mean
    0.4510 * (log(alb) - 5.137)
  
  return(PI)
}

## Compute 2y KFRE risk (P) from the linear predictor, using the 2-year
## baseline survival probability (0.9832) from the original KFRE derivation.
KFRE_risk_2y <- function(PI) {
  # 1 - S0(t)^exp(PI): standard Cox-model risk transform of PI into a
  # predicted probability of kidney failure by 2 years
  risk_2y <- 1 - (0.9832 ** exp(PI))
  
  return(risk_2y)
}

## Compute 5y KFRE risk (P) from the linear predictor, using the 5-year
## baseline survival probability (0.9365) from the original KFRE derivation.
KFRE_risk_5y <- function(PI) {
  # same risk transform as KFRE_risk_2y(), with the 5-year baseline survival
  risk_5y <- 1 - (0.9365 ** exp(PI))
  
  return(risk_5y)
}

## ============================================================================
## Discrimination
## ============================================================================

## Core call-site(s) for point estimates from riskRegression::Score().
## AUC and Brier/scaled-Brier (IPA) share ONE call with `cens.method =
## "ipcw"` (inverse probability of censoring weighting); calibration
## intercept/slope gets its own separate call with `cens.method = "pseudo"`
## (jackknife pseudo-values), which it needs for its pseudo-observations.
## These can't be combined into one call: Score() only accepts one
## `cens.method` per call, AND IPA specifically only resolves under "ipcw"
## -- under "pseudo" the null/reference model IPA is scaled against isn't
## computed, so it silently comes back NA even though the raw Brier score
## itself is fine. So this makes up to two Score() calls: one ipcw call
## (only when "auc" and/or "brier" is requested), one pseudo call (only
## when `calibration` is TRUE). Either call is skipped entirely when not
## needed -- an AUC-only request would make just one Score() call, though in
## practice every caller in this file (including the sensitivity analyses,
## via evaluate_predictions()/validate_predictions()) asks for all three.
##
## Every function that needs AUC, calibration, and/or Brier point estimates
## (`validate_predictions()` below, which `evaluate_predictions()` and, in
## turn, the sensitivity-analysis functions further down all build on) is a
## thin wrapper around this -- so this is the only place
## `riskRegression::Score()` is called for point estimates in this file. The
## only other place it's called is inside `bootstrap_brier_oe()` below,
## because a bootstrap CI needs its own resampling loop and can't share this
## call site (its `cens.method` is "ipcw" too, matching the Brier/IPA point
## estimate here -- for the same reason: IPA needs "ipcw" to resolve).
##
## `data` is a data frame (typically `cohort`, or a subset/recoded copy of
## it built by the caller) rather than pre-extracted time/status vectors, so
## this works both for the main results table and for sensitivity analyses
## that subset and/or recode the outcome first -- the caller does that
## subsetting/recoding on `data`'s columns before calling. time_auc/
## status_auc/time_restricted/status_restricted are each derived from
## `data` exactly once, right where they're needed below (AUC's inside the
## `need_auc` block, restricted's inside the `need_brier`/`calibration`
## block), instead of being computed upfront by every caller.
score_predictions <- function(pred_risks,
                              data,
                              horizon,
                              cause = 1,
                              conf.int = TRUE,
                              metrics = c("auc", "brier"),
                              calibration = TRUE) {
  # Hist() is exported by prodlim, but riskRegression's formula interface
  # expects it unqualified (Hist(...) ~ 1) -- bind it into scope explicitly
  # rather than relying on prodlim being attached via library()
  Hist <- prodlim::Hist
  
  results <- list()
  
  need_auc <- "auc" %in% metrics
  need_brier <- "brier" %in% metrics
  
  # --- Discrimination (AUC), via IPCW censoring handling, on the UNCAPPED
  # time/status. Kept as its own Score() call (rather than shared with
  # Brier below) specifically so it can use a different time basis than
  # Brier/calibration: AUC uses the UNCAPPED time_auc/status_auc (i.e.
  # time_to_event_inf/outcome_inf) to avoid tying the risk set exactly at
  # `times = horizon * 365.25` below (which happens when using the
  # horizon-capped time_to_event_<horizon>y, since it's truncated at that
  # same value for anyone reaching administrative end-of-follow-up -- this
  # starves the AUC control definition of subjects genuinely observed past
  # the horizon). Derived from `data` here, right before use, rather than
  # upfront by the caller. ---
  if (need_auc) {
    time_auc <- data[["time_to_event_inf"]]
    status_auc <- data[["outcome_inf"]]
    # riskRegression::Score() wants time-to-event/outcome in a data frame
    # with fixed column names ("time_to_event", "outcome"), so wrap here.
    score_data_auc <- data.frame(time_to_event = time_auc, outcome = status_auc)
    
    Score_auc <- riskRegression::Score(
      list(pred_risks),
      formula = Hist(time_to_event, outcome) ~ 1,
      cens.method = "ipcw",
      data = score_data_auc,
      times = horizon * 365.25,
      # define event of interest (was `outcome = 1`, which is not a valid
      # riskRegression::Score() argument and silently did nothing -- the
      # correct argument is `cause`)
      cause = cause,
      conf.int = conf.int,
      metrics = "auc"
    )
    
    # Score()'s AUC table has one row per model; "numeric" is riskRegression's
    # internal name for a plain numeric risk-prediction model (as opposed to
    # a fitted regression model), which is what `list(pred_risks)` becomes
    auc_model <- Score_auc$AUC$score[model == "numeric"]
    results$AUC <- auc_model$AUC
    results$AUC_CI <- if (conf.int) {
      c(Lower = auc_model$lower, Upper = auc_model$upper)
    } else {
      NULL
    }
  }
  
  # --- Brier score / scaled Brier (IPA), via IPCW censoring handling, on the
  # horizon-capped time/status (kept restricted, unlike AUC, per explicit
  # instruction). Both come from ONE ipcw call: empirically, Score()'s raw
  # Brier number is IDENTICAL whether requested under cens.method = "ipcw"
  # or "pseudo" (verified directly -- both are equally affected by the tie
  # at exactly `horizon * 365.25` when time is horizon-capped, unlike
  # calibration below, which genuinely is tie-invariant). So there's no
  # benefit to computing Brier from the pseudo/calibration call -- it's
  # just as biased there, and splitting it out only adds a redundant Score()
  # call. IPA requires "ipcw" specifically (it needs a properly-computed
  # null/reference model as its denominator, which doesn't get computed
  # under cens.method = "pseudo" -- Score() silently returns NA for IPA
  # there even though the raw Brier score itself resolves).
  #
  # conf.int is hardcoded FALSE here (NOT the `conf.int` argument) -- verified
  # directly that requesting Score()'s own analytic CI together with
  # summary = "ipa" makes the IPA POINT ESTIMATE itself come back NA (while
  # se/lower/upper get populated for Brier, not IPA -- an riskRegression
  # quirk, unrelated to censoring method or horizon-capping). This analytic
  # CI is never used anyway: `evaluate_predictions()` always overwrites
  # Brier_CI/Scaled_Brier_CI with the percentile bootstrap from
  # `bootstrap_brier_oe()` when `B > 0` -- so there's no cost to forcing it
  # off here. ---
  #
  # time_restricted/status_restricted (the horizon-capped columns) are
  # derived from `data` here, right before their first use -- shared by
  # both this Brier block and the calibration block below, so computed once
  # rather than by every caller.
  if (need_brier || calibration) {
    time_restricted <- data[[paste0("time_to_event_", horizon, "y")]]
    status_restricted <- data[[paste0("outcome_", horizon, "y")]]
    score_data_restricted <- data.frame(time_to_event = time_restricted, outcome = status_restricted)
  }
  
  if (need_brier) {
    Score_brier <- riskRegression::Score(
      list(pred_risks),
      formula = Hist(time_to_event, outcome) ~ 1,
      cens.method = "ipcw",
      data = score_data_restricted,
      times = horizon * 365.25,
      cause = cause,
      conf.int = FALSE,
      metrics = "brier",
      summary = "ipa"
    )
    
    # same "numeric" model-row convention as the AUC block above
    brier_model <- Score_brier$Brier$score[model == "numeric"]
    results$Brier_Score <- brier_model$Brier
    results$Scaled_Brier_Score <- brier_model$IPA
  }
  
  # --- Calibration intercept & slope, via pseudo-value censoring handling,
  # using unsmoothed pseudo-observations, on the horizon-capped time/status.
  # Unlike Brier above, this genuinely IS tie-invariant (verified directly):
  # the pseudo-value regression that drives Intercept/Slope is identical
  # whether time is horizon-capped, capped-but-tie-broken, or fully
  # uncapped -- only Score()'s reported Brier/IPA numbers are tie-sensitive,
  # not the pseudo-values themselves. ---
  if (calibration) {
    Score_cb <- riskRegression::Score(
      list(pred_risks),
      formula = Hist(time_to_event, outcome) ~ 1,
      cens.method = "pseudo",
      data = score_data_restricted,
      times = horizon * 365.25,
      cause = cause,
      conf.int = conf.int,
      plots = "calibration",
      # riskRegression's calibration plotframe is only populated when
      # metrics includes "brier" -- request it here purely to get that
      # plotframe; its Brier/IPA numbers are not used (see above)
      metrics = "brier"
    )
    
    # one pseudo-observation per patient: an unbiased estimate of their
    # individual outcome probability at `horizon`, used as the "observed"
    # side of the calibration regression below
    pseudos <- data.frame(Score_cb$Calibration$plotframe)
    
    # add the cloglog risk estimates (complementary log-log transform of the
    # predicted risk) -- calibration intercept/slope are estimated on this
    # scale, which is the natural link function for survival/CIF models
    pseudos$cll_pred <- log(-log(1 - pseudos$risk))
    
    # fit model for calibration intercept: regress the pseudo-observations on
    # an offset of the cloglog predicted risk (offset = coefficient fixed at
    # 1); the estimated mean/intercept term is the calibration-in-the-large
    fit_cal_int <- geepack::geese(
      pseudovalue ~ offset(cll_pred),
      data = pseudos,
      id = riskRegression_ID,
      scale.fix = TRUE,
      family = gaussian,
      mean.link = "cloglog",
      corstr = "independence",
      jack = TRUE
    )
    Intercept <- summary(fit_cal_int)$mean$estimate
    Intercept_SE <- summary(fit_cal_int)$mean$san.se
    
    # fit model for calibration slope: same idea, but now cll_pred is also
    # included as a covariate (not just an offset), so its coefficient
    # estimates the DEVIATION of the true slope from 1 (perfect calibration).
    # offset(cll_pred) centers the slope to zero, therefore later we need to add one
    fit_cal_slope <- geepack::geese(
      pseudovalue ~ offset(cll_pred) + cll_pred,
      data = pseudos,
      id = riskRegression_ID,
      scale.fix = TRUE,
      family = gaussian,
      mean.link = "cloglog",
      corstr = "independence",
      jack = TRUE
    )
    Slope_est <- summary(fit_cal_slope)$mean["cll_pred", ]$estimate
    Slope_SE <- summary(fit_cal_slope)$mean["cll_pred", ]$san.se
    
    results$pseudos <- pseudos
    results$Intercept <- Intercept
    # 95% CI via a normal approximation (estimate +/- 1.96 * robust SE)
    results$Intercept_CI <- c(
      Lower = Intercept - qnorm(0.975) * Intercept_SE,
      Upper = Intercept + qnorm(0.975) * Intercept_SE
    )
    # add 1 back, since the offset centered the slope coefficient at zero
    # (i.e. Slope_est = 0 means perfect calibration, slope = 1)
    results$Slope <- 1 + Slope_est
    results$Slope_CI <- c(
      Lower = 1 + (Slope_est - qnorm(0.975) * Slope_SE),
      Upper = 1 + (Slope_est + qnorm(0.975) * Slope_SE)
    )
  }
  
  return(results)
}

## ============================================================================
## Calibration, Brier score, and discrimination (single riskRegression::Score() call)
## ============================================================================

## Compute discrimination (AUC), calibration intercept & slope (with 95%
## CIs), and the Brier score / scaled Brier (IPA) point estimates, via
## `score_predictions()` above. This is a thin wrapper that just extracts
## `time`/`status` from `data`.
##
## AUC uses the UNCAPPED `time_to_event_inf`/`outcome_inf` columns; Brier and
## calibration use the horizon-capped `time_to_event_<horizon>y`/
## `outcome_<horizon>y` ones. The capped columns are truncated at exactly
## `horizon * 365.25` days for anyone reaching administrative
## end-of-follow-up, which is the same value passed as `times` to Score()
## below -- so most of the risk set is tied exactly at the evaluation
## horizon. That tie starves the AUC control definition of subjects
## genuinely observed past the horizon, so AUC is computed on `_inf` instead
## (uncapped follow-up: earliest of event/death/admin-censoring/emigration
## only), which removes the tie while `times = horizon * 365.25` still
## evaluates at the same horizon. Brier is left on the capped columns
## (matching how it's specified/reported elsewhere) even though it shares
## AUC's IPCW estimator and is subject to the same tie-driven bias -- see
## `score_predictions()` for that caveat. Calibration is unaffected either
## way (its cens.method = "pseudo" call doesn't do IPCW weighting), so it's
## also left on the capped columns.
validate_predictions <- function(pred_risks, data, horizon, conf.int = TRUE) {
  # delegate to the shared Score() call site, asking for everything
  # (AUC + calibration + Brier/scaled Brier) in one go; score_predictions()
  # derives the uncapped time/status (AUC) and horizon-capped time/status
  # (Brier, calibration) from `data` itself
  score_predictions(
    pred_risks = pred_risks,
    data = data,
    horizon = horizon,
    cause = 1,
    conf.int = conf.int,
    metrics = c("auc", "brier"),
    calibration = TRUE
  )
}

## Compute bootstrap CIs for the Brier score, scaled Brier (IPA), and O/E
## ratio TOGETHER, from ONE shared resampling loop -- rather than two
## separate boot::boot() calls (previously one here, one inside oe_ratio()).
## Those two calls only happened to draw the same patient indices because
## they were seeded and sized identically; nothing actually tied them
## together, so if B or the sample size ever differed between them, the
## Brier and O/E bootstraps would silently stop lining up. Sharing one loop
## here removes that fragility, and also fixes a real limitation the old
## O/E bootstrap had: since each replicate here has the resampled patients'
## own time/status data (not just their resampled risk predictions), the
## CIF can be refit on each resample via cmprsk::cuminc() -- rather than
## reusing the ORIGINAL, unresampled CIF for every replicate (the old
## behavior), which only captured resampling variability on the "expected"
## (mean predicted risk) side and likely made the O/E CI too narrow.
bootstrap_brier_oe <- function(pred_risks,
                               data,
                               horizon,
                               B = 500,
                               seed = 123) {
  # extract outcome variables according to horizon, and store them under
  # fixed names ("time_to_event"/"outcome") so the resampled data still has
  # the column names the Score() formula below expects -- horizon-capped
  # columns, matching the Brier/O-E point estimates elsewhere in this file.
  data$time_to_event <- data[[paste0("time_to_event_", horizon, "y")]]
  data$outcome <- data[[paste0("outcome_", horizon, "y")]]
  
  Hist <- prodlim::Hist
  
  # one bootstrap replicate: resample patients (with replacement) via
  # `indices`, then recompute Brier/IPA (via Score()) AND O/E (via a
  # freshly-refit CIF) on that SAME resampled sample
  boot_func <- function(data, indices, risk_var) {
    resampled_data <- data[indices, ]
    resampled_risk <- risk_var[indices]
    
    boot_brier <- riskRegression::Score(
      list(resampled_risk),
      formula = Hist(time_to_event, outcome) ~ 1,
      cens.method = "ipcw",
      data = resampled_data,
      conf.int = FALSE,
      times = horizon * 365.25,
      cause = 1,
      metrics = "brier",
      summary = "ipa"
    )
    numeric_rows <- boot_brier$Brier$score$model == "numeric"
    brier <- as.numeric(boot_brier$Brier$score$Brier[numeric_rows])
    ipa <- as.numeric(boot_brier$Brier$score$IPA[numeric_rows])
    
    # O/E ratio on this SAME resample: expected = mean of the resampled
    # predicted risks; observed = cumulative incidence from a CIF refit on
    # the resampled patients' own time/outcome data (see comment above this
    # function for why this is refit rather than reused).
    boot_cif <- cmprsk::cuminc(
      ftime = resampled_data$time_to_event,
      fstatus = resampled_data$outcome,
      cencode = 0
    )
    boot_observed <- cmprsk::timepoints(boot_cif, horizon * 365.25)
    oe <- boot_observed$est[1] / mean(resampled_risk)
    
    # boot::boot() wants a single named numeric vector back per replicate;
    # Brier, IPA, and O/E CIs are all derived from these same replicates
    c(brier = brier, ipa = ipa, oe = oe)
  }
  
  set.seed(seed)
  # runs boot_func() B times, each time with a fresh bootstrap resample (a
  # new `indices` vector) of the original data -- the SAME resample drives
  # Brier, IPA, and O/E every replicate, since all three are computed
  # inside the one boot_func() call above.
  boot_results <- boot::boot(
    data = data,
    statistic = boot_func,
    R = B,
    risk_var = pred_risks
  )
  
  # percentile bootstrap CI: 2.5th/97.5th percentiles of the B replicate
  # values (index = 1 -> Brier, 2 -> IPA, 3 -> O/E, matching the column
  # order boot_func() returns above)
  extract_ci <- function(index) {
    ci <- boot::boot.ci(boot_results, type = "perc", index = index)
    if (is.null(ci))
      NULL
    else
      c(Lower = ci$percent[4], Upper = ci$percent[5])
  }
  
  list(
    Brier_CI = extract_ci(1),
    Scaled_Brier_CI = extract_ci(2),
    O_E_CI = extract_ci(3)
  )
}

## Convenience wrapper: discrimination (AUC), calibration intercept/slope
## (+CI), Brier score, scaled Brier (IPA), and O/E ratio -- point estimates
## for ALL validation measures, plus bootstrap CIs for Brier/scaled-Brier/O-E
## (computed together, from one shared resampling loop -- see
## `bootstrap_brier_oe()`). This replaces the previous separate calls to
## `cal_int_slope()`, `brier_scores()`, a standalone timeROC AUC computation,
## and a separate `oe_ratio()` call.
##
## `CIF` is the cumulative incidence function object needed for the O/E
## ratio's point-estimate "observed" side (see `oe_ratio()` below) -- it
## isn't produced by Score(), so it still has to be passed in separately;
## everything else comes from the one Score() call inside
## `validate_predictions()`, plus the one bootstrap loop in
## `bootstrap_brier_oe()`.
evaluate_predictions <- function(pred_risks,
                                 data,
                                 horizon,
                                 CIF,
                                 B = 500,
                                 seed = 123) {
  # step 1: point estimates for AUC, calibration intercept/slope, Brier, and
  # scaled Brier, all from one Score() call (also includes the analytic
  # calibration/AUC CIs, since those don't need bootstrapping)
  results <- validate_predictions(pred_risks, data, horizon, conf.int = TRUE)
  
  # step 2: O/E ratio point estimate (needs `CIF`, which Score() doesn't
  # produce, so this is always its own call)
  results$O_E_Ratio <- oe_ratio(pred_risks = pred_risks, CIF = CIF, horizon = horizon)
  
  # step 3: bootstrap CIs for Brier, scaled Brier, and O/E TOGETHER, from one
  # shared resampling loop -- see `bootstrap_brier_oe()` for why this
  # replaces what used to be two independent bootstraps.
  if (B > 0) {
    boot_ci <- bootstrap_brier_oe(pred_risks, data, horizon, B = B, seed = seed)
    results$Brier_CI <- boot_ci$Brier_CI
    results$Scaled_Brier_CI <- boot_ci$Scaled_Brier_CI
    results$O_E_CI <- boot_ci$O_E_CI
  }
  
  return(results)
}

## ============================================================================
## O/E ratio
## ============================================================================

## Point estimate for the O/E ratio:
##   Expected = mean predicted risk across patients
##   Observed = cumulative incidence (from the fitted CIF) at `horizon`
## The bootstrap CI for O/E is no longer computed here -- see
## `bootstrap_brier_oe()` above, which computes it together with Brier/IPA
## from one shared resampling loop, refitting the CIF on each resample
## rather than reusing this function's original, unresampled `CIF` (which is
## what a standalone O/E bootstrap here would otherwise be limited to).
oe_ratio <- function(pred_risks, CIF, horizon) {
  expected <- mean(pred_risks)
  observed <- cmprsk::timepoints(CIF, horizon * 365.25)
  
  return(observed$est[1] / expected)
}

## ============================================================================
## Reclassification tables
## ============================================================================

## Create a reclassification table comparing a reference risk category
## (already present in `data`, e.g. `cat_risk_2y_ckd_epi_2009_cr`) against a
## new set of `risk_categories`, for a given set of category `levels`.
## This replaces the previous `generate_reclass_table_2y()` and
## `generate_reclass_table_5y()`, which were identical apart from the
## reference column name and category labels.
generate_reclass_table <- function(data,
                                   reference_col,
                                   # name (string) of the reference risk-category column in `data`
                                   risk_categories,
                                   # new risk categories to compare against the reference
                                   levels,
                                   # category labels, e.g. c(">40%", "<=40%") or c("<3%", "3% - 5%", ">5%")
                                   title) {
  # Ensure all levels are present in the data (so empty categories still
  # show up as a row/column of zeros, rather than being silently dropped)
  reference <- factor(data[[reference_col]], levels = levels)
  risk_categories <- factor(risk_categories, levels = levels)
  
  # Create a table of reclassification: rows = reference-equation category,
  # columns = new-equation category. Off-diagonal cells are patients whose
  # risk category changed between the two equations.
  table_raw <- table(reference, risk_categories)
  # Convert to percentages
  table_percent <- prop.table(table_raw) * 100
  # Format table to include both counts and percentages
  table_final <- paste0(table_raw, " (", round(table_percent, 1), "%)")
  # Convert to matrix for correct printing
  dim(table_final) <- dim(table_raw)
  dimnames(table_final) <- dimnames(table_raw)
  # Return the formatted table using kable for markdown output
  formatted_table <- knitr::kable(table_final, align = "c", caption = title) |>
    kableExtra::kable_styling(
      bootstrap_options = c("striped", "hover", "condensed", "responsive"),
      full_width = FALSE
    ) |>
    kableExtra::row_spec(0, bold = TRUE, background = "#D3D3D3") # Header row formatting
  
  return(list(table = table_final, formatted_table = formatted_table))
}

## ============================================================================
## Validation pipeline steps
## ============================================================================
## Everything below wraps one step of the validation pipeline (one figure or
## table) as a function: it takes the cohort + a few settings and does the
## whole step, including writing its output file(s). "5_ Validation.R" is now
## just a short, ordered sequence of calls to these functions, so almost all
## of the actual logic lives here and can be tested / edited / reused in one
## place instead of being spread across a 1400+ line script.
## ============================================================================

## ---- Step: main performance measures table ---------------------------------

## Return the permutation of `cols` (a vector of "<Measure>_<horizon>y[_
## <suffix>]" column names) that groups them by horizon first, then by
## subgroup suffix (e.g. "_omit1"/"_omit2", "_subgroup1".."_subgroup3" --
## "" for the main results table, which has no subgroup), then by measure
## (in the order given by `measure_order`). Each measure's point estimate
## and its 95% CI now live together in ONE cell (via `fmt_est_ci()`), so
## there's no separate "_CI" column left to interleave -- `is_ci` below is
## kept only as a harmless no-op tie-breaker in case a caller still passes
## a legacy "..._CI" column name (e.g. from code outside this file that
## hasn't been updated to the combined-cell format yet). E.g. applying this
## to c(Int_2y, Int_5y, Slope_2y, Slope_5y) yields Int_2y, Slope_2y, Int_5y,
## Slope_5y, rather than grouped by measure across all horizons first -- and
## for a sensitivity table's Int_2y_omit1, Int_2y_omit2, Slope_2y_omit1,
## Slope_2y_omit2, ... it keeps each suffix's own measures grouped together
## (rather than every suffix's Int first, then every suffix's Slope).
## Returns an integer index permutation (like `order()`), NOT the reordered
## names themselves -- used at `write_measures_xlsx()`'s call sites as
## `cols[order_measure_cols(colnames(...)[cols], ...)]`.
order_measure_cols <- function(cols, measure_order) {
  # horizon is the run of digits immediately before the "y" in each name
  # (e.g. "2" from "Int_2y_omit1"); measure is everything before
  # "_<horizon>y" (e.g. "Int"); suffix is whatever sits after "<horizon>y"
  # (e.g. "_omit1"; "" for the main table's plain "Int_2y")
  horizon <- as.numeric(sub("^.*?_([0-9]+)y.*$", "\\1", cols))
  measure <- sub("_[0-9]+y.*$", "", cols)
  is_ci <- grepl("_CI$", cols)
  suffix <- sub("_CI$", "", sub("^.*?_[0-9]+y", "", cols))
  measure_rank <- match(measure, measure_order)
  
  order(horizon, suffix, measure_rank, is_ci)
}

## Write a wide "measures" results table (one row per equation, plus one "N"
## row holding each horizon['s subgroup]'s sample size, columns like AUC_2y,
## Int_5y, Slope_5y, Brier_2y, ... -- each equation's cell already holding
## "estimate (lower; upper)" via `fmt_est_ci()`, not split across separate
## estimate/CI columns) to an Excel file split across topic-based tabs
## instead of one flat sheet: "Calibration" (intercept + slope), "OE" (O/E
## ratio, its own tab even though it's also a calibration measure), "AUC",
## and "Brier" (Brier + scaled Brier). Every tab is ordered horizon-major
## (and, for sensitivity tables, subgroup-major within horizon) via
## `order_measure_cols()` above, so e.g. Calibration reads Int_2y, Slope_2y,
## Int_5y, Slope_5y, ... instead of every horizon's Int first, then every
## horizon's Slope.
##
## Sample size is shown as a SECOND HEADER ROW ("N = 26505"), shown ONCE per
## horizon[/subgroup] block and merged across that block's columns (e.g.
## "Int_2y"/"Slope_2y" share one merged "N = 27125" cell, since both are
## computed on the same patients) -- matched via the
## "N_<horizon>y[_<suffix>]" columns in the results table's "N" row (written
## by `compute_performance_measures()`/`run_sensitivity_core()`), rather
## than N appearing as its own data column. This applies to every tab,
## including OE (O/E shares the exact same "N" row as every other measure --
## there's only ever one N per horizon/subgroup for the whole table, not a
## separate one per measure). Tables/tabs with no "N" row at all just get a
## blank second header row. Used by every function in this file that
## exports a measures table (`compute_performance_measures()`,
## `run_sensitivity_outcome_comparison()`, `run_sensitivity_egfr_subgroup()`),
## so this is the single place that layout is defined.
write_measures_xlsx <- function(results_df, file, rowNames = TRUE) {
  ordered_cols <- function(regex, measure_order) {
    cols <- grep(regex, colnames(results_df))
    cols[order_measure_cols(colnames(results_df)[cols], measure_order)]
  }
  
  calibration_cols <- ordered_cols("^(Int|Slope)_", c("Int", "Slope"))
  auc_cols <- ordered_cols("^AUC_", "AUC")
  brier_cols <- ordered_cols("^(Brier|Scaled_Brier)_", c("Brier", "Scaled_Brier"))
  
  col_groups <- list(
    Calibration = calibration_cols,
    OE          = grep("^OE_", colnames(results_df)),
    AUC         = auc_cols,
    Brier       = brier_cols
  )
  
  has_n_row <- "N" %in% rownames(results_df)
  data_rows <- setdiff(rownames(results_df), "N")
  
  # for a measure column like "Int_2y_omit1", find the value of its
  # matching "N_2y_omit1" column in the "N" row -- reuses the exact same
  # horizon/suffix parsing order_measure_cols() uses, so a column always
  # matches the N_ column for its own horizon and subgroup, never a
  # different one. Returns just the bare number (e.g. "27125"), not a
  # "N = ..." label -- the row this feeds into already has "N" as its own
  # row name, so repeating "N = " inside every cell would be redundant.
  n_value_for_col <- function(col) {
    if (!has_n_row)
      return("")
    horizon <- sub("^.*?_([0-9]+)y.*$", "\\1", col)
    suffix  <- sub("^.*?_[0-9]+y", "", col)
    n_col <- paste0("N_", horizon, "y", suffix)
    if (n_col %in% colnames(results_df) &&
        !is.na(results_df["N", n_col])) {
      as.character(results_df["N", n_col])
    } else {
      ""
    }
  }
  
  wb <- openxlsx::createWorkbook()
  for (sheet_name in names(col_groups)) {
    cols <- col_groups[[sheet_name]]
    # skip a tab entirely if this results_df has none of its columns (e.g.
    # a results_df built with B = 0 has no Brier/O-E CIs but still has the
    # point-estimate columns, so only a table missing a measure category
    # altogether -- not just its CI -- skips that tab)
    if (length(cols) == 0)
      next
    openxlsx::addWorksheet(wb, sheet_name)
    
    col_names <- colnames(results_df)[cols]
    # blank corner cell(s) above the row-name column, if row names are shown
    corner <- if (rowNames)
      ""
    else
      character(0)
    
    # row 1: the column titles themselves
    openxlsx::writeData(
      wb,
      sheet_name,
      x = matrix(c(corner, col_names), nrow = 1),
      startRow = 1,
      colNames = FALSE
    )
    
    # row 2: an actual "N" row (row name "N" in the leftmost column, bare
    # sample-size numbers in the measure columns) -- shown ONCE per
    # horizon[/subgroup] block rather than under every individual measure
    # column, e.g. Int_2y and Slope_2y share the same N (both computed on
    # the same patients for horizon 2), so "27125" only needs to appear
    # once, not repeated under both. This applies uniformly to every tab,
    # including OE. Blank only when there's genuinely no matching N_ value
    # (e.g. a table built without an "N" row at all).
    block_key <- function(col) {
      horizon <- sub("^.*?_([0-9]+)y.*$", "\\1", col)
      suffix  <- sub("^.*?_[0-9]+y", "", col)
      paste0(horizon, suffix)
    }
    keys <- vapply(col_names, block_key, character(1))
    values <- vapply(col_names, n_value_for_col, character(1))
    # TRUE at the first column of each consecutive run of the same block
    # key -- relies on `cols` already being horizon-major/suffix-major
    # (via order_measure_cols() above), so one block's columns are always
    # adjacent, never interleaved with another block's
    is_block_start <- c(TRUE, keys[-1] != keys[-length(keys)])
    n_row <- ifelse(is_block_start, values, "")
    n_row_label <- if (rowNames) "N" else character(0)
    openxlsx::writeData(
      wb,
      sheet_name,
      x = matrix(c(n_row_label, n_row), nrow = 1),
      startRow = 2,
      colNames = FALSE
    )
    
    # data rows (excluding the "N" row, which is written separately above),
    # starting right after the two header rows -- row names are built into
    # the matrix explicitly (rather than via writeData's rowNames = TRUE)
    # so there's no dependence on how that argument interacts with
    # colNames = FALSE
    body <- as.matrix(results_df[data_rows, cols, drop = FALSE])
    body_mat <- if (rowNames)
      cbind(data_rows, body)
    else
      body
    openxlsx::writeData(
      wb,
      sheet_name,
      x = body_mat,
      startRow = 3,
      colNames = FALSE
    )
  }
  openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
  
  return(invisible(wb))
}

## Format a point estimate together with its 95% CI as one string, e.g.
## fmt_est_ci(0.9284, c(Lower = 0.9123, Upper = 0.9451)) -> "0.928 (0.912;
## 0.945)" -- used everywhere an estimate + CI used to be written into two
## separate "<Measure>_<horizon>y[...]" / "<Measure>_<horizon>y[...]_CI"
## columns, so the pair now lands in a single Excel cell instead. `ci` is
## the c(Lower=, Upper=) vector `score_predictions()`/`evaluate_predictions()`
## return; `scale` lets the scaled-Brier's percentage values (already *100)
## share this helper without a separate code path.
fmt_est_ci <- function(est,
                       ci,
                       digits = 3,
                       scale = 1) {
  fmt <- paste0("%.", digits, "f")
  sprintf(paste0(fmt, " (", fmt, "; ", fmt, ")"),
          est * scale,
          ci["Lower"] * scale,
          ci["Upper"] * scale)
}

## Column names for one measurement "block" (N, AUC, calibration intercept &
## slope, O/E ratio, Brier, scaled Brier) for a given horizon, optionally
## with a subgroup `suffix` appended for sensitivity tables (e.g. "omit1",
## "subgroup2" -> "AUC_2y_omit1"). `suffix = NULL`/`""` (the default) is the
## main results table's plain columns (e.g. "AUC_2y", not "AUC_2y_" with a
## trailing underscore).
##
## This is the SINGLE source of truth for this column layout, shared by
## `compute_performance_measures()` (main table, no suffix) below and
## `sa_measure_colnames()` (sensitivity tables, with suffix) further down --
## so the two can never again drift out of sync the way `5__Validation.R`'s
## hand-maintained `measures` vector once did: it kept separate "..._CI"
## columns around long after `fmt_est_ci()` replaced them with one combined
## "estimate (CI)" cell, leaving those columns silently blank in every
## exported Excel file. Generating the column list here, from the same
## place that writes into it, makes that class of bug structurally
## impossible rather than something to remember to keep in sync by hand.
measure_colnames <- function(horizon, suffix = NULL) {
  tag <- if (is.null(suffix) || suffix == "")
    ""
  else
    paste0("_", suffix)
  c(
    paste0("N_", horizon, "y", tag),
    paste0("AUC_", horizon, "y", tag),
    paste0("Int_", horizon, "y", tag),
    paste0("Slope_", horizon, "y", tag),
    paste0("OE_", horizon, "y", tag),
    paste0("Brier_", horizon, "y", tag),
    paste0("Scaled_Brier_", horizon, "y", tag)
  )
}

## Discrimination (AUC), calibration (intercept & slope), Brier / scaled
## Brier, and O/E ratio for every equation and horizon -> writes an Excel
## table and returns it (invisibly). Column names come from
## `measure_colnames()` above -- not from a caller-supplied argument -- so
## this table's layout can't drift out of sync with what this function
## actually writes into it.
compute_performance_measures <- function(cohort,
                                         horizons,
                                         equations,
                                         model_names,
                                         B = 500,
                                         file = "Measures.xlsx") {
  # empty results table: one row per equation (+ one "N" row for sample
  # size per horizon, same convention as run_sensitivity_core() below), one
  # column per measure (populated below, then relabeled to `model_names` at
  # the end)
  col_names <- unlist(lapply(horizons, measure_colnames))
  results_df <- data.frame(matrix(nrow = length(equations) + 1, ncol = length(col_names)))
  rownames(results_df) <- c("N", equations)
  colnames(results_df) <- col_names
  
  # loop over every (horizon, equation) combination -- e.g. 2 horizons x 5
  # equations = 10 calls to evaluate_predictions() total
  for (horizon in horizons) {
    for (equation in equations) {
      # e.g. horizon = 2, equation = "ckd_epi_2009_cr" -> cohort$risk_2y_ckd_epi_2009_cr
      pred_risks <- cohort[[paste0("risk_", horizon, "y_", equation)]]
      # the cumulative incidence function object for this horizon (e.g.
      # `cin_2y`), expected to already exist in the environment (loaded
      # alongside `cohort` from cohort_predictions.RData)
      CIF <- get(paste0("cin_", horizon, "y"))
      
      # all point estimates + CIs for this equation/horizon, in one call
      performance <- evaluate_predictions(
        pred_risks = pred_risks,
        data = cohort,
        horizon = horizon,
        CIF = CIF,
        B = B
      )
      
      # record the horizon's sample size once (same across equations, since
      # every equation is scored on the same `cohort`) -- its own dedicated
      # "N_<horizon>y" column, in the pre-declared "N" row, so
      # write_measures_xlsx() can show it under that horizon's column titles
      if (equation == equations[1]) {
        results_df["N", paste0("N_", horizon, "y")] <- nrow(cohort)
      }
      
      # write each measure + its 95% CI into ONE "<Measure>_<horizon>y" cell
      # via fmt_est_ci() (e.g. "0.928 [0.912; 0.945]"), instead of a plain
      # estimate column plus a separate "..._CI" column
      results_df[equation, paste0("AUC_", horizon, "y")] <-
        fmt_est_ci(performance$AUC, performance$AUC_CI)
      results_df[equation, paste0("Int_", horizon, "y")] <-
        fmt_est_ci(performance$Intercept, performance$Intercept_CI)
      results_df[equation, paste0("Slope_", horizon, "y")] <-
        fmt_est_ci(performance$Slope, performance$Slope_CI)
      results_df[equation, paste0("OE_", horizon, "y")] <-
        sprintf("%.3f", performance$O_E_Ratio)
      
      # Brier/scaled-Brier/O-E CIs are only available when bootstrapping was
      # requested (B = 0 skips the bootstrap CIs for a quick/no-CI run) --
      # fall back to a plain estimate (no brackets) when there's no CI to
      # combine it with
      if (B > 0) {
        results_df[equation, paste0("Brier_", horizon, "y")] <-
          fmt_est_ci(performance$Brier_Score, performance$Brier_CI)
        results_df[equation, paste0("Scaled_Brier_", horizon, "y")] <-
          fmt_est_ci(
            performance$Scaled_Brier_Score,
            performance$Scaled_Brier_CI,
            digits = 1,
            scale = 100
          )
        results_df[equation, paste0("OE_", horizon, "y")] <-
          fmt_est_ci(performance$O_E_Ratio, performance$O_E_CI)
      } else {
        results_df[equation, paste0("Brier_", horizon, "y")] <-
          sprintf("%.3f", performance$Brier_Score)
        results_df[equation, paste0("Scaled_Brier_", horizon, "y")] <-
          sprintf("%.1f", performance$Scaled_Brier_Score * 100)
      }
    }
  }
  
  # swap the internal equation codes for their display names before export
  # (the "N" row keeps its own name -- see write_measures_xlsx())
  rownames(results_df) <- c("N", model_names)
  write_measures_xlsx(results_df, file = file)
  
  return(invisible(results_df))
}

## ---- Step: calibration plots (+ risk histogram insets) ---------------------

## Predicted-risk histogram for one horizon/trim combination, for a single
## equation (used as an inset under the combined calibration plot).
build_risk_histogram <- function(cohort,
                                 horizon,
                                 trim,
                                 max_5y,
                                 color = "darkorange4",
                                 xlab = "CKD-EPIcr 2009") {
  var <- paste0("risk_", horizon, "y_ckd_epi_2009_cr")
  # narrower bins for the trimmed (zoomed-in) view, since it covers a much
  # smaller x-axis range and needs finer resolution to still look like a
  # histogram rather than 1-2 wide bars. Untrimmed binwidth is the same for
  # both horizons; only the trimmed binwidth differs by horizon.
  binwidth <- if (!trim) {
    0.025
  } else if (horizon == 2) {
    0.005
  } else {
    0.01
  }
  
  hist <- ggplot2::ggplot(cohort, ggplot2::aes(x = .data[[var]])) +
    ggplot2::geom_histogram(
      binwidth = binwidth,
      boundary = 0,
      fill = color,
      color = "white"
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::labs(
      title = NULL,
      y = "Count",
      x = paste("Predicted Risk using", xlab)
    ) +
    ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      plot.background = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = 14),
      axis.title = ggplot2::element_text(size = 14)
    )
  
  if (trim) {
    # zoom the x-axis to the 95th percentile of 5-year risk (`max_5y`),
    # so the trimmed view isn't dominated by a long tail of high-risk outliers
    hist <- hist + ggplot2::coord_cartesian(xlim = c(0, max_5y))
  }
  
  return(hist)
}

## Combined calibration plot (all equations, one panel) for both horizons,
## trimmed and untrimmed, with a predicted-risk histogram strip underneath.
## Writes both PNGs.
plot_calibration_curves <- function(cohort,
                                    horizons,
                                    equations,
                                    model_names,
                                    colors_palette,
                                    max_5y,
                                    file_full,
                                    file_trimmed) {
  # the 4 histogram insets (2y/5y x trimmed/untrimmed) for the reference
  # equation (CKD-EPIcr 2009), stacked underneath the calibration curves below
  histogram_risk <- list(
    risk_2y_ckd_epi_2009_cr = build_risk_histogram(cohort, 2, FALSE, max_5y),
    risk_2y_ckd_epi_2009_cr_trimmed = build_risk_histogram(cohort, 2, TRUE, max_5y),
    risk_5y_ckd_epi_2009_cr = build_risk_histogram(cohort, 5, FALSE, max_5y),
    risk_5y_ckd_epi_2009_cr_trimmed = build_risk_histogram(cohort, 5, TRUE, max_5y)
  )
  
  # predicted risk for every equation, one list per horizon, in `equations` order
  pred_list_2y <- lapply(equations, function(eq)
    cohort[[paste0("risk_2y_", eq)]])
  pred_list_5y <- lapply(equations, function(eq)
    cohort[[paste0("risk_5y_", eq)]])
  
  # build one calibration-curve plot per (trim, horizon) combination -- 4
  # plots total (2 horizons x trimmed/untrimmed) -- each overlaying all
  # equations. Stored in a list keyed by trim/horizon (rather than assigned
  # to loose variables named "combined_cal_plot_5_trimmed" etc.), so they can
  # be looked up explicitly below instead of relying on those names existing
  # in scope. Loop order (trim outer, horizon inner) doesn't matter for the
  # result.
  combined_cal_plots_by <- list(full = list(), trimmed = list())
  
  for (trim in c(TRUE, FALSE)) {
    trim_key <- if (trim) "trimmed" else "full"
    for (horizon in horizons) {
      calibration_data <- data.frame()
      pred_list <- if (horizon == 2)
        pred_list_2y
      else
        pred_list_5y
      
      # for each equation: get its calibration pseudo-observations, smooth
      # them against predicted risk, and append to the shared long-format
      # data frame used for the multi-line plot below
      for (nr_equation in seq_along(pred_list)) {
        pseudos <- validate_predictions(
          pred_risks = pred_list[[nr_equation]],
          data = cohort,
          horizon = horizon,
          conf.int = TRUE
        )$pseudos |>
          dplyr::arrange(risk)
        
        # loess-smooth the raw (noisy) pseudo-observations against predicted
        # risk, to get a smooth calibration curve rather than a scatter
        smooth_pseudos <- predict(stats::loess(
          pseudovalue ~ risk,
          data = pseudos,
          degree = 1,
          span = 0.33
        ))
        
        temp_df <- data.frame(risk = pseudos$risk,
                              observed = smooth_pseudos,
                              model = model_names[nr_equation])
        calibration_data <- dplyr::bind_rows(calibration_data, temp_df)
      }
      
      # fix the legend/line order to match `model_names`, not alphabetical
      calibration_data$model <- factor(calibration_data$model, levels = model_names)
      
      # one line per equation: predicted risk (x) vs. smoothed observed risk (y)
      combined_cal_plot <- ggplot2::ggplot(calibration_data,
                                           ggplot2::aes(x = risk, y = observed, color = model)) +
        ggplot2::geom_line(linewidth = 0.75, alpha = 0.8) +
        ggplot2::scale_color_manual(values = colors_palette, breaks = model_names) +
        ggplot2::labs(
          title = paste0(horizon, "-year KFRE"),
          x = "Predicted Risks",
          y = "Observed Risks"
        ) +
        ggthemes::theme_clean() +
        ggplot2::theme(
          legend.title = ggplot2::element_blank(),
          legend.background = ggplot2::element_rect(colour = NA),
          legend.position = "bottom",
          plot.subtitle = ggplot2::element_text(size = 10),
          panel.border = ggplot2::element_blank(),
          plot.background = ggplot2::element_blank(),
          axis.text = ggplot2::element_text(size = 14),
          axis.title = ggplot2::element_text(size = 14)
        )
      
      # add the diagonal reference line (perfect calibration) and zoom the
      # axes. NOTE: as written, this produces the SAME xlim/ylim (0, max_5y)
      # for both the trimmed and untrimmed plots -- only the invisible
      # diagonal segment's endpoint differs (max_5y vs. 1), which has no
      # visible effect once cropped by coord_cartesian. If the untrimmed
      # plot is meant to show the full 0-1 risk scale (as the comment this
      # replaced suggested), the untrimmed branch's coord_cartesian needs
      # xlim/ylim = c(0, 1) instead -- flagging rather than changing this,
      # since it changes the figure's appearance.
      if (trim) {
        combined_cal_plot <- combined_cal_plot +
          ggplot2::annotate(
            "segment",
            x = 0,
            y = 0,
            xend = max_5y,
            yend = max_5y,
            linetype = "dashed",
            color = "gray40"
          ) +
          ggplot2::coord_cartesian(xlim = c(0, max_5y),
                                   ylim = c(0, max_5y))
      } else {
        combined_cal_plot <- combined_cal_plot +
          ggplot2::annotate(
            "segment",
            x = 0,
            y = 0,
            xend = 1,
            yend = 1,
            linetype = "dashed",
            color = "gray40"
          ) +
          ggplot2::coord_cartesian(xlim = c(0, max_5y),
                                   ylim = c(0, max_5y))
      }
      
      combined_cal_plots_by[[trim_key]][[as.character(horizon)]] <- combined_cal_plot
    }
  }
  
  # --- Untrimmed figure: 2y + 5y calibration curves, with histograms below ---
  combined_cal_plots <-
    (combined_cal_plots_by$full[["2"]] + combined_cal_plots_by$full[["5"]]) +
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = "bottom")
  combined_hist_plots <-
    histogram_risk$risk_2y_ckd_epi_2009_cr + histogram_risk$risk_5y_ckd_epi_2009_cr
  
  ggplot2::ggsave(
    filename = file_full,
    plot = (combined_cal_plots / combined_hist_plots) + patchwork::plot_layout(heights = c(3, 1)),
    width = 9,
    height = 6,
    dpi = 300
  )
  
  # --- Trimmed figure: same layout, using the trimmed/zoomed plots instead ---
  combined_cal_plots_trimmed <-
    (combined_cal_plots_by$trimmed[["2"]] + combined_cal_plots_by$trimmed[["5"]]) +
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = "bottom")
  combined_hist_plots_trimmed <-
    histogram_risk$risk_2y_ckd_epi_2009_cr_trimmed + histogram_risk$risk_5y_ckd_epi_2009_cr_trimmed
  
  ggplot2::ggsave(
    filename = file_trimmed,
    plot = (combined_cal_plots_trimmed / combined_hist_plots_trimmed) + patchwork::plot_layout(heights = c(3, 1)),
    width = 9,
    height = 6,
    dpi = 300
  )
  
  return(invisible(NULL))
}

## ---- Step: discrimination forest plot ---------------------------------------

## Forest plot of the time-dependent AUC (with 95% CI) per equation, for both
## horizons -- read straight from "Measures.xlsx" (produced by
## `compute_performance_measures()`), since that already has the AUC + CI
## columns this plot needs. Writes the PNG and returns the plot (invisibly).
plot_discrimination_forest <- function(measures_file, out_file) {
  # Measures.xlsx now stores each estimate together with its CI in ONE cell,
  # formatted as "0.928 (0.912; 0.945)" (via fmt_est_ci()) -- this pulls the
  # point estimate and the two CI bounds back out of that single string.
  parse_est_ci <- function(x) {
    est <- as.numeric(sub("\\s*\\(.*$", "", x))
    bounds <- sub("^[^(]*\\(", "", x)
    bounds <- gsub("\\)", "", bounds)
    parts <- strsplit(bounds, ";\\s*")
    data.frame(est = est,
               lo  = as.numeric(sapply(parts, `[`, 1)),
               hi  = as.numeric(sapply(parts, `[`, 2)))
  }
  
  # Measures.xlsx is now written with write_measures_xlsx() as separate
  # Calibration/AUC/Brier tabs (see that function above) -- read the "AUC"
  # sheet explicitly rather than relying on sheet order
  td_AUC <- openxlsx::read.xlsx(measures_file, sheet = "AUC")
  # the first (unnamed) column from writeData(rowNames = TRUE) holds the
  # model names -- give it an explicit name to work with
  colnames(td_AUC)[1] <- "Equation"
  # write_measures_xlsx() now writes a second header row ("N = ...") right
  # under the column titles, which read.xlsx() has no way to distinguish
  # from a data row -- drop it here (it isn't per-equation data anyway).
  td_AUC <- td_AUC[-1, ]
  
  # keep the original combined text as the plot label ("0.928 [0.912;
  # 0.945]"), and parse out the numeric estimate/lo/hi for plotting
  lab_2y <- td_AUC$AUC_2y
  lab_5y <- td_AUC$AUC_5y
  p_2y <- parse_est_ci(lab_2y)
  p_5y <- parse_est_ci(lab_5y)
  
  td_AUC <- td_AUC |>
    dplyr::mutate(
      AUC_2y = p_2y$est,
      AUC_5y = p_5y$est,
      lo_2y  = p_2y$lo,
      hi_2y  = p_2y$hi,
      lo_5y  = p_5y$lo,
      hi_5y  = p_5y$hi,
      lab_2y = lab_2y,
      lab_5y = lab_5y,
      row_num  = dplyr::row_number(),
      # alternate shading every other row, for readability
      shade    = row_num %% 2 == 0,
      # reverse the row order so the first equation in the table ends up at
      # the TOP of the plot (ggplot draws factor level 1 at the bottom)
      Equation = factor(Equation, levels = rev(Equation))
    )
  
  # reshape from one row per equation (wide: AUC_2y, AUC_5y side by side) to
  # one row per equation x horizon (long), so both horizons can share one
  # plot via facet_wrap() below
  td_AUC_long <- td_AUC |>
    tidyr::pivot_longer(
      cols = c(AUC_2y, AUC_5y),
      names_to = "horizon",
      values_to = "AUC"
    ) |>
    dplyr::mutate(
      lo    = dplyr::if_else(horizon == "AUC_2y", lo_2y, lo_5y),
      hi    = dplyr::if_else(horizon == "AUC_2y", hi_2y, hi_5y),
      label = dplyr::if_else(horizon == "AUC_2y", lab_2y, lab_5y),
      horizon = factor(
        horizon,
        levels = c("AUC_2y", "AUC_5y"),
        labels = c("2-year AUC (95% CI)", "5-year AUC (95% CI)")
      )
    )
  
  # --- Dynamic x-axis window, derived from the actual AUC/CI values ---
  # rather than numbers tuned for one particular dataset. AUC's theoretical
  # ceiling is 1.00 (also where the "perfect discrimination" reference line
  # sits), so only the lower bound needs to adapt to the data.
  step <- 0.02
  data_lo <- min(td_AUC_long$lo, na.rm = TRUE)
  # two grid steps of padding below the lowest CI lower bound (rounded onto
  # the `step` grid), rather than one -- with only one step, the leftmost
  # whisker/point can end up sitting right on top of the equation-name text
  # when its CI comes close to the axis edge
  x_lower <- floor(data_lo / step) * step - 2 * step
  x_upper <- 1.00
  
  x_limits <- c(x_lower, x_upper)
  x_breaks <- seq(x_lower, x_upper, by = step)
  
  # how far the panel extends beyond x_upper, to leave room for the
  # "AUC [CI]" text labels. Rather than guessing a fraction of the panel
  # width (which either left dead space or clipped the labels, depending on
  # how it was tuned -- the text has a fixed PHYSICAL width, which doesn't
  # scale neatly with the data range either way), this measures the actual
  # rendered width of the widest label via grid, then solves for the exact
  # `panel_max` that gives it just enough room, plus a small safety cushion.
  fig_width <- 6.4 # must match the `width` passed to ggsave() below
  
  # ggplot2's `size` aesthetic is in mm; convert to points the same way
  # ggplot2 does internally (1 pt = 25.4/72.27 mm) so the measurement uses
  # the same font size the plot will actually render at
  label_size_pt <- 2.8 * (72.27 / 25.4)
  widest_label <- td_AUC_long$label[which.max(nchar(td_AUC_long$label))]
  
  # measure in a throwaway null graphics device, so this works reliably
  # regardless of what device (if any) is already open in the session
  grDevices::pdf(NULL)
  label_width_in <- grid::convertWidth(grid::grobWidth(grid::textGrob(
    widest_label, gp = grid::gpar(fontsize = label_size_pt)
  )), "inches", valueOnly = TRUE)
  grDevices::dev.off()
  
  small_gap <- step * 0.25 # tiny buffer between the AUC = 1.00 line and where the label text starts
  x_right   <- x_upper + small_gap
  # 20% safety cushion, since this measurement won't perfectly match the
  # final saved plot's font rendering
  text_in <- label_width_in * 1.2
  # solve panel_max so the measured label width -- converted to data units,
  # which itself depends on panel_max, since a wider panel means more
  # inches per data unit -- fits exactly:
  #   text_in = (text_data_units / (panel_max - x_limits[1])) * fig_width
  panel_max <- (x_right - (text_in / fig_width) * x_limits[1]) / (1 - text_in / fig_width)
  
  x_left <- x_limits[1] # equation-name labels start exactly at the panel's left edge
  n_eq   <- nlevels(td_AUC_long$Equation)
  
  # build the alternating-row shading rectangles: one row per shaded
  # equation/horizon combination, spanning the full plot width
  shade_df <- td_AUC_long |>
    dplyr::filter(shade) |>
    dplyr::mutate(y_pos = as.integer(Equation),
                  ymin = y_pos - 0.5,
                  ymax = y_pos + 0.5) |>
    dplyr::distinct(horizon, Equation, .keep_all = TRUE)
  
  # forest plot: point + horizontal CI per equation, faceted by horizon
  AUC_plot <- ggplot2::ggplot(td_AUC_long,
                              ggplot2::aes(
                                y = Equation,
                                x = AUC,
                                xmin = lo,
                                xmax = hi
                              )) +
    # alternating row shading, drawn first so points/text render on top.
    # Uses the exact numeric panel bounds (x_limits[1] to panel_max) rather
    # than xmin = -Inf/xmax = Inf -- mixing scale_x_continuous(limits = ...)
    # (hard limits, drops out-of-range data) with coord_cartesian(xlim = ...)
    # (a visual zoom) means -Inf/Inf don't reliably resolve to the same
    # visible edges those two systems computed separately, which left the
    # shading not reaching all the way to the left edge behind the labels.
    ggplot2::geom_rect(
      data = shade_df,
      ggplot2::aes(ymin = ymin, ymax = ymax),
      xmin = x_limits[1],
      xmax = panel_max,
      inherit.aes = FALSE,
      fill = "lightgrey",
      color = NA
    ) +
    # dashed reference line at AUC = 1.00 (perfect discrimination)
    ggplot2::geom_vline(
      xintercept = x_limits[2],
      color = "black",
      linewidth = 0.4,
      linetype = "dashed"
    ) +
    # the CI whisker and point estimate for each equation
    ggplot2::geom_errorbar(
      width = 0.15,
      orientation = "y",
      linewidth = 0.4,
      color = "black"
    ) +
    ggplot2::geom_point(shape = 15,
                        size = 2,
                        color = "black") +
    # "AUC [CI]" text, to the right of the panel (allowed to overflow the
    # panel border because of coord_cartesian(clip = "off") below)
    ggplot2::geom_text(
      ggplot2::aes(x = x_right, label = label),
      hjust = 0,
      size = 2.8,
      color = "black"
    ) +
    # equation name, drawn as text INSIDE the panel at the far left (the
    # real y-axis labels are hidden below, so this replaces them -- doing it
    # this way lets the shading rectangles render behind the text)
    ggplot2::geom_text(
      ggplot2::aes(x = x_left, y = Equation, label = Equation),
      hjust = 0,
      size = 2.8,
      color = "black",
      inherit.aes = FALSE
    ) +
    ggplot2::scale_x_continuous(
      limits = c(x_limits[1], panel_max),
      # sprintf() forces two decimals on every break (e.g. "0.70", "1.00")
      # -- plain x_breaks[-1] would let R's default numeric formatting drop
      # trailing zeros (e.g. "0.7", "1")
      labels = c("", sprintf("%.2f", x_breaks[-1])),
      breaks = c(x_breaks[2], x_breaks[-1]),
      expand = c(0, 0)
    ) +
    # clip = "off" lets the equation-name and AUC-label text spill outside
    # the strict panel area; xlim extends a bit past x_upper to leave room
    # for the "AUC [CI]" labels on the right
    ggplot2::coord_cartesian(
      clip = "off",
      xlim = c(x_limits[1], panel_max),
      ylim = c(0.5, n_eq + 0.5)
    ) +
    ggplot2::facet_wrap( ~ horizon, ncol = 1) +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(
      axis.line.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(color = "black"),
      panel.grid = ggplot2::element_blank(),
      legend.position = "none",
      # margin set to zero -- all the room for the "AUC [CI]" text (drawn at
      # x_right, outside the panel via clip = "off") now comes from
      # `panel_max`/`label_fraction` above instead of a physical margin
      plot.margin = ggplot2::margin(0, 0, 0, 0),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(
        face = "bold",
        size = 10,
        hjust = 1,
        color = "black"
      )
    ) +
    # the real axis line/ticks are hidden above (axis.line.x = element_blank()),
    # so this manually draws a short x-axis segment under the bottom facet only
    ggplot2::geom_segment(
      data = data.frame(
        x = x_breaks[2],
        xend = x_upper,
        y = -0.1,
        yend = -0.1,
        horizon = levels(td_AUC_long$horizon)[2]
      ),
      ggplot2::aes(
        x = x,
        xend = xend,
        y = y,
        yend = yend
      ),
      inherit.aes = FALSE,
      color = "black",
      linewidth = 0.5
    )
  
  suppressMessages(ggplot2::ggsave(
    out_file,
    AUC_plot,
    width = fig_width,
    # same variable used to size panel_max above -- keep these in sync
    height = 4,
    dpi = 300
  ))
  
  return(invisible(AUC_plot))
}

## ---- Step: decision curve analysis (DCA) ------------------------------------

## Net benefit plot + table across equations, for both horizons. The raw DCA
## computation is slow, so it's cached to `cache_file` and skipped unless
## `new_DCA_run = TRUE`.
run_dca <- function(cohort,
                    horizons,
                    equations,
                    model_names,
                    colors_palette,
                    max_5y,
                    new_DCA_run = FALSE,
                    cache_file,
                    plot_file,
                    table_file) {
  # dcurves::dca() needs a factor outcome (rather than the 0/1/2 numeric
  # coding used elsewhere), labelled for its competing-risks handling
  cohort <- cohort |>
    dplyr::mutate(
      outcome_2y_dca = factor(
        outcome_2y,
        levels = 0:2,
        labels = c("censor", "KFRT", "Death")
      ),
      outcome_5y_dca = factor(
        outcome_5y,
        levels = 0:2,
        labels = c("censor", "KFRT", "Death")
      )
    )
  
  if (new_DCA_run) {
    # --- Slow path: recompute DCA from scratch for every equation/horizon ---
    for (horizon in horizons) {
      cat("Time horizon", horizon, "years\n")
      # max threshold probability to evaluate net benefit at (2-year risk is
      # generally higher, hence the wider range)
      max_tp <- if (horizon == 2)
        0.5
      else
        0.2
      
      for (nr_equation in seq_along(equations)) {
        cat("Equation", equations[nr_equation], "\n")
        # e.g. Surv(time = time_to_event_2y, event = outcome_2y_dca) ~ risk_2y_ckd_epi_2009_cr
        formula <- as.formula(
          paste0(
            "survival::Surv(time = time_to_event_",
            horizon,
            "y,\n",
            "        event = outcome_",
            horizon,
            "y_dca) ~ risk_",
            horizon,
            "y_",
            equations[nr_equation]
          )
        )
        
        # restrict to the eGFR range where a 2-year (resp. 5-year) prediction
        # is clinically relevant, matching how DCA is normally reported for
        # this kind of risk equation
        res_cohort <- if (horizon == 2) {
          cohort[cohort$ckd_epi_2021_cr >= 10 & cohort$ckd_epi_2021_cr < 30, ]
        } else {
          cohort[cohort$ckd_epi_2021_cr >= 30 & cohort$ckd_epi_2021_cr < 60, ]
        }
        
        # run DCA across a grid of threshold probabilities from 0 to max_tp
        dca <- dcurves::dca(
          formula,
          data = res_cohort,
          thresholds = seq(0, max_tp, 0.01),
          time = horizon * 365
        )
        assign(paste0("dca_", horizon, "y_", nr_equation), dca)
        
        table <- dplyr::as_tibble(dca$dca)
        assign(paste0("table_dca_", horizon, "y_", nr_equation),
               table)
      }
    }
    # cache all 20 DCA objects/tables (2 horizons x 5 equations x 2 kinds) to
    # disk, so future runs can skip straight to the fast path below
    save(
      dca_2y_1,
      dca_2y_2,
      dca_2y_3,
      dca_2y_4,
      dca_2y_5,
      dca_5y_1,
      dca_5y_2,
      dca_5y_3,
      dca_5y_4,
      dca_5y_5,
      table_dca_2y_1,
      table_dca_2y_2,
      table_dca_2y_3,
      table_dca_2y_4,
      table_dca_2y_5,
      table_dca_5y_1,
      table_dca_5y_2,
      table_dca_5y_3,
      table_dca_5y_4,
      table_dca_5y_5,
      file = cache_file
    )
  } else {
    # --- Fast path: load the previously cached DCA objects/tables ---
    load(cache_file)
  }
  
  # build one net-benefit plot per horizon, from the (freshly computed or
  # cached) per-equation DCA tables
  for (horizon in horizons) {
    # plot-specific axis limits and vertical reference lines (clinically
    # meaningful threshold probabilities), tuned per horizon
    if (horizon == 2) {
      max_tp <- max_5y
      vlines <- c(0.1, 0.4)
      min_nb <- -0.02
      max_nb <- 0.14
    } else {
      max_tp <- max_5y
      vlines <- c(0.03, 0.05)
      min_nb <- -0.005
      max_nb <- 0.016
    }
    
    # assemble a wide "net benefit" table: one row per threshold, one column
    # per equation, plus the "treat all" / "treat none" reference strategies
    for (nr_equation in seq_along(equations)) {
      table_dca <- get(paste0("table_dca_", horizon, "y_", nr_equation))
      
      # "treat all"/"treat none" don't depend on the equation, so only need
      # to be pulled out once (from the first equation's table)
      if (nr_equation == 1) {
        treat_none <- table_dca |> dplyr::filter(label == "Treat None") |> dplyr::select(threshold, net_benefit)
        treat_all <- table_dca |> dplyr::filter(label == "Treat All") |> dplyr::select(net_benefit)
        NB <- data.frame(
          thresholds = treat_none$threshold,
          treat_none = treat_none$net_benefit,
          treat_all = treat_all$net_benefit
        )
      }
      
      # net benefit for this specific equation, added as its own column
      NB_dca <- table_dca |>
        dplyr::filter(label == paste0("risk_", horizon, "y_", equations[nr_equation])) |>
        dplyr::select(net_benefit)
      NB[, equations[nr_equation]] <- NB_dca$net_benefit
    }
    # keep the full wide table around (used later to build Net Benefit Table.xlsx)
    assign(paste0("NB_", horizon, "y"), NB)
    
    # long format for ggplot: one row per (threshold, strategy) pair
    NB_long <- NB |> tidyr::pivot_longer(
      cols = -thresholds,
      names_to = "variable",
      values_to = "value"
    )
    
    # x-axis breaks: ggplot's default "nice" breaks, plus the clinically
    # meaningful vlines thresholds, so they always get their own tick mark
    max_x <- max(max_tp, max(vlines, na.rm = TRUE))
    default_breaks <- scales::extended_breaks()(c(0, max_x))
    custom_x_breaks <- sort(unique(c(default_breaks, vlines)))
    custom_y_breaks <- seq(min(min_nb), max_nb, ifelse(horizon == 2, 0.02, 0.005))
    
    # one line per strategy (treat all / treat none / each equation)
    dca_plot <- ggplot2::ggplot(NB_long,
                                ggplot2::aes(x = thresholds, y = value, color = variable)) +
      ggplot2::geom_line() +
      ggplot2::scale_x_continuous(breaks = custom_x_breaks, labels = scales::percent) +
      ggplot2::scale_y_continuous(breaks = custom_y_breaks) +
      ggplot2::coord_cartesian(ylim = c(min_nb, max_nb)) +
      ggplot2::labs(
        title = paste0(horizon, "-year KFRE CKD-EPI equations"),
        subtitle = paste0(
          "For patients with eGFRcr 2009 ",
          ifelse(horizon == 2, 10, 30),
          "\u2013",
          ifelse(horizon == 2, 29, 59),
          " mL/min/1.73m\u00b2 (N=",
          unique(get(paste0(
            "table_dca_", horizon, "y_1"
          ))$n),
          ")"
        ),
        x = "Threshold Probability",
        y = "Net Benefit"
      ) +
      ggplot2::theme(
        legend.title = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(size = 10),
        legend.margin = ggplot2::margin(
          t = 0,
          r = 0,
          b = 0,
          l = -20
        ),
        plot.subtitle = ggplot2::element_text(size = 10),
        plot.background = ggplot2::element_rect(fill = "white", color = NA),
        panel.background = ggplot2::element_rect(fill = "white", color = NA),
        panel.border = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor.x = ggplot2::element_blank(),
        panel.grid.minor.y = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_line(
          color = "black",
          linetype = "dotted",
          linewidth = 0.5
        ),
        axis.line.x = ggplot2::element_line(),
        axis.line.y = ggplot2::element_line(),
        axis.ticks = ggplot2::element_line(),
        plot.title = ggplot2::element_text(face = "bold"),
        axis.text = ggplot2::element_text(size = 14),
        axis.title = ggplot2::element_text(size = 14)
      ) +
      ggplot2::scale_color_manual(
        values = c("darkviolet", "black", colors_palette),
        labels = c("Treat All", "Treat None", model_names),
        breaks = c("treat_all", "treat_none", equations)
      ) +
      ggplot2::guides(color = ggplot2::guide_legend(ncol = 7, byrow = TRUE)) +
      ggplot2::geom_vline(
        xintercept = vlines[1],
        linetype = "dashed",
        color = "black"
      ) +
      ggplot2::geom_vline(
        xintercept = vlines[2],
        linetype = "dashed",
        color = "black"
      )
    
    assign(paste0("dca_plot_", horizon, "y"), dca_plot)
  }
  
  # combine the 2y and 5y net-benefit plots side by side, sharing one legend
  combined_dca_plot <- (dca_plot_2y + dca_plot_5y) +
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = "bottom")
  
  ggplot2::ggsave(
    filename = plot_file,
    plot = combined_dca_plot,
    width = 10,
    height = 4,
    dpi = 300
  )
  
  # net-benefit summary table: just the values at the clinically meaningful
  # threshold(s) picked out via `vlines` above (0.1/0.4 for 2y, 0.03/0.05 for 5y)
  NB_table <- rbind(data.frame(time = 2, NB_2y[NB_2y$thresholds %in% c(0.1, 0.4), ]), data.frame(time = 5, NB_5y[NB_5y$thresholds %in% c(0.03, 0.05), ]))
  openxlsx::write.xlsx(NB_table, rowNames = FALSE, file = table_file)
  
  return(invisible(NULL))
}

## ---- Step: predicted risk distributions & risk differences -----------------

## Long-format predicted risk (all equations) for one horizon, with equations
## ordered/labelled for the ridge plots (reverse of `equations`, so the
## reference equation ends up at the bottom).
build_risk_distribution_long <- function(cohort, horizon, equations, model_names) {
  # e.g. horizon = 2 -> c("risk_2y_ckd_epi_2009_cr", "risk_2y_ckd_epi_2021_cr", ...)
  cols <- paste0("risk_", horizon, "y_", equations)
  # reversed so the FIRST equation in `equations` ends up at the BOTTOM of
  # the ridge plots below (ggplot draws factor level 1 at the bottom)
  equation_order <- rev(cols)
  
  cohort |>
    dplyr::select(lopnr, dplyr::all_of(cols)) |>
    # reshape from one column per equation (wide) to one row per patient x
    # equation (long), which ggridges needs for the density-ridge plots
    tidyr::pivot_longer(
      cols = dplyr::all_of(cols),
      names_to = "Equation",
      values_to = "Risk"
    ) |>
    dplyr::mutate(Equation = factor(
      Equation,
      levels = equation_order,
      labels = rev(model_names)
    ))
}

## Build one ridge-density plot (0-100% + a zoomed inset) for a given
## long-format risk data frame.
build_risk_distribution_plot <- function(risk_long,
                                         horizon,
                                         zoom_limit,
                                         zoom_fill_limit,
                                         zoom_breaks,
                                         title_full,
                                         title_zoom,
                                         x_lab_full = "Predicted risk") {
  # the "full" ridge plot: density of predicted risk per equation, across the
  # entire 0-100% range
  plot_full <- ggplot2::ggplot(risk_long,
                               ggplot2::aes(
                                 x = Risk,
                                 y = Equation,
                                 fill = ggplot2::after_stat(x)
                               )) +
    ggridges::geom_density_ridges_gradient(scale = 1,
                                           rel_min_height = 0.000001,
                                           alpha = 0.8) +
    ggplot2::scale_fill_viridis_c(
      name = "Risk",
      option = "C",
      limits = c(0, 0.05),
      oob = scales::squish
    ) +
    ggplot2::scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
    ggplot2::scale_y_discrete(expand = c(0, 0)) +
    ggplot2::labs(title = title_full, x = x_lab_full, y = "Equation") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.spacing = ggplot2::unit(0.5, "lines"),
      axis.text.y = ggplot2::element_text(size = 10)
    )
  
  # the "zoomed" ridge plot: same data, but restricted to the low-risk range
  # (`zoom_limit`) where most of the cohort actually sits -- this is what
  # gets inset into the top-right corner of the full plot below
  plot_zoom <- ggplot2::ggplot(risk_long,
                               ggplot2::aes(
                                 x = Risk,
                                 y = Equation,
                                 fill = ggplot2::after_stat(x)
                               )) +
    ggridges::geom_density_ridges_gradient(scale = 0.9,
                                           rel_min_height = 0.000001,
                                           alpha = 0.8) +
    ggplot2::scale_fill_viridis_c(
      name = "Risk",
      option = "C",
      limits = c(0, zoom_fill_limit),
      oob = scales::squish
    ) +
    ggplot2::scale_x_continuous(limits = c(0, zoom_limit), breaks = zoom_breaks) +
    ggplot2::scale_y_discrete(expand = c(0, 0)) +
    ggplot2::labs(title = title_zoom, x = "Risk Score", y = "Equation") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = ggplot2::alpha("white", 0.8), color = NA),
      panel.background = ggplot2::element_rect(fill = ggplot2::alpha("white", 0.8), color = NA),
      legend.title = ggplot2::element_text(size = 8),
      legend.text = ggplot2::element_text(size = 6),
      legend.key.size = ggplot2::unit(0.5, "lines"),
      panel.spacing = ggplot2::unit(0.5, "lines"),
      axis.text.y = ggplot2::element_text(size = 10),
      axis.text.x = ggplot2::element_text(size = 8),
      axis.title.x = ggplot2::element_text(size = 8),
      axis.title.y = ggplot2::element_text(size = 8),
      plot.title = ggplot2::element_text(size = 8)
    )
  
  # overlay the zoomed plot as an inset in the top-right of the full plot
  plot_full +
    ggplot2::annotation_custom(
      ggplot2::ggplotGrob(plot_zoom),
      xmin = 0.4,
      xmax = 1.05,
      ymin = 3,
      ymax = 6
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = 8),
      legend.text = ggplot2::element_text(size = 6),
      plot.tag = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 10),
      axis.text = ggplot2::element_text(size = 10),
      axis.text.y = ggplot2::element_text(size = 10),
      plot.background = ggplot2::element_rect(fill = "white", color = "white"),
      panel.background = ggplot2::element_rect(fill = "white", color = "white")
    )
}

## Combined 2y + 5y predicted-risk distribution plot. Writes the PNG.
plot_risk_distributions <- function(risk_2y_long, risk_5y_long, out_file) {
  # each horizon gets its own zoom window/fill scale, tuned to where most of
  # that horizon's predicted risks actually fall (2y risks run higher than 5y)
  plot_2y <- build_risk_distribution_plot(
    risk_2y_long,
    2,
    zoom_limit = 0.005,
    zoom_fill_limit = 0.004,
    zoom_breaks = seq(0, 0.005, by = 0.001),
    title_full = "2-year KFRE risk 0-100%",
    title_zoom = "2-year KFRE risk 0-0.5%"
  )
  plot_5y <- build_risk_distribution_plot(
    risk_5y_long,
    5,
    zoom_limit = 0.01,
    zoom_fill_limit = 0.008,
    zoom_breaks = seq(0, 0.01, by = 0.002),
    title_full = "5-year KFRE risk 0-100%",
    title_zoom = "5y KFRE Risk Distributions by Equation 0-1%",
    x_lab_full = "Risk Score"
  )
  
  # stack 2y above 5y, sharing one legend
  combined_dist_plot <- plot_2y + plot_5y + patchwork::plot_layout(ncol = 1, guides = "collect")
  
  ggplot2::ggsave(
    filename = out_file,
    plot = combined_dist_plot,
    width = 15,
    height = 10,
    dpi = 300
  )
  
  return(invisible(NULL))
}

## Risk differences vs. a reference equation, for one horizon's long-format
## risk data frame.
build_risk_difference_plot <- function(risk_long,
                                       horizon,
                                       reference_label,
                                       title,
                                       x_lab) {
  # per patient, subtract the reference equation's predicted risk from every
  # other equation's; drop the (now all-zero) reference equation afterward
  risk_diff_long <- risk_long |>
    dplyr::group_by(lopnr) |>
    dplyr::mutate(risk_diff = Risk - Risk[Equation == reference_label]) |>
    dplyr::ungroup() |>
    dplyr::filter(Equation != reference_label)
  
  ggplot2::ggplot(
    risk_diff_long,
    ggplot2::aes(
      x = risk_diff,
      y = Equation,
      fill = ggplot2::after_stat(x)
    )
  ) +
    ggridges::geom_density_ridges_gradient(scale = 1,
                                           rel_min_height = 0.00001,
                                           alpha = 0.8) +
    ggplot2::scale_fill_viridis_c(
      name = "Risk Difference",
      option = "C",
      limits = c(-0.0025, 0.0025),
      oob = scales::squish
    ) +
    ggplot2::scale_x_continuous(limits = c(-0.01, 0.01),
                                breaks = seq(-0.01, 0.01, by = 0.001)) +
    ggplot2::scale_y_discrete(expand = c(0, 0)) +
    ggplot2::geom_vline(xintercept = 0,
                        linetype = "dashed",
                        color = "black") +
    ggplot2::labs(title = title, x = x_lab, y = "Equation") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.spacing = ggplot2::unit(0.5, "lines"),
      axis.text.y = ggplot2::element_text(size = 10)
    )
}

## Combined 2y + 5y risk-difference plot vs. a reference equation
## (default "CKD-EPIcr 2009"). Writes the PNG.
plot_risk_differences <- function(risk_2y_long,
                                  risk_5y_long,
                                  reference_label = "CKD-EPIcr 2009",
                                  out_file) {
  # same reference equation used for both horizons; titles/axis labels are
  # built here since they reference `reference_label` directly
  plot_2y <- build_risk_difference_plot(
    risk_2y_long,
    2,
    reference_label,
    title = paste("2y KFRE risk difference compared to", reference_label),
    x_lab = paste("Risk difference (Risk equation -", reference_label, ")")
  )
  plot_5y <- build_risk_difference_plot(
    risk_5y_long,
    5,
    reference_label,
    title = paste("5y KFRE Risk Differences Compared to", reference_label),
    x_lab = paste("Risk Difference (Risk Equation -", reference_label, ")")
  )
  
  # stack 2y above 5y, sharing one legend
  combined_plot_diff <- plot_2y + plot_5y + patchwork::plot_layout(ncol = 1, guides = "collect")
  
  ggplot2::ggsave(
    filename = out_file,
    plot = combined_plot_diff,
    width = 15,
    height = 10,
    dpi = 300
  )
  
  return(invisible(NULL))
}

## ---- Step: eGFR distributions -----------------------------------------------

## Distribution of eGFR (the predictor, not predicted risk) across equations.
## Writes the PNG.
plot_egfr_distributions <- function(cohort, equations, model_names, out_file) {
  # `equations` here are also the eGFR column names themselves (e.g.
  # "ckd_epi_2009_cr" is both the equation code and the eGFR column) -- one
  # row per patient x equation, in the order given by `equations`
  cohort_egfr_long <- cohort |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(equations),
      names_to = "Equation",
      values_to = "eGFR"
    ) |>
    dplyr::mutate(Equation = factor(Equation, levels = equations, labels = model_names))
  
  egfr_plot <- ggplot2::ggplot(cohort_egfr_long,
                               ggplot2::aes(
                                 x = eGFR,
                                 y = Equation,
                                 fill = ggplot2::after_stat(x)
                               )) +
    ggridges::geom_density_ridges_gradient(scale = 0.9,
                                           rel_min_height = 0.001,
                                           alpha = 0.8) +
    ggplot2::scale_fill_viridis_c(
      name = "eGFR",
      option = "C",
      limits = c(0, 60),
      oob = scales::squish
    ) +
    ggplot2::scale_x_continuous(limits = c(0, 60), breaks = seq(0, 60, by = 10)) +
    ggplot2::scale_y_discrete(limits = rev, expand = c(0, 0)) +
    ggplot2::labs(title = "eGFR distributions", x = "eGFR (mL/min/1.73m\u00b2)", y = "Equation") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.spacing = ggplot2::unit(1.5, "lines"),
      axis.text.y = ggplot2::element_text(size = 10)
    ) +
    ggplot2::theme(plot.margin = ggplot2::margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0
    ))
  
  ggplot2::ggsave(
    filename = out_file,
    plot = egfr_plot,
    width = 15,
    height = 10,
    dpi = 300
  )
  
  return(invisible(NULL))
}

## Summary statistics (n, min/max, 5th/95th percentile, median, and counts +
## percentages below 10 / at-or-above 60 mL/min/1.73m^2) for one eGFR
## equation's values. `x` is the eGFR vector for a single equation (e.g.
## `cohort$ckd_epi_2009_cr`); `name` labels the resulting row and defaults
## to `x`'s own variable name, so `egfr_summary(cohort$ckd_epi_2009_cr)`
## needs no extra argument. Missing values are dropped before summarizing.
## Returns a one-row data frame; row-bind across equations (e.g.
## `dplyr::bind_rows(lapply(equations, function(eq) egfr_summary(cohort[[eq]],
## eq)))`) for a table covering all of them.
egfr_summary <- function(x, name = deparse(substitute(x))) {
  x <- x[!is.na(x)]
  
  data.frame(
    equation      = name,
    n             = length(x),
    min           = round(min(x), 1),
    p5            = round(quantile(x, 0.05, names = FALSE), 1),
    median        = round(median(x), 1),
    p95           = round(quantile(x, 0.95, names = FALSE), 1),
    max           = round(max(x), 1),
    n_below_10    = sum(x < 10),
    pct_below_10  = round(100 * sum(x < 10) / length(x), 1),
    n_above_60    = sum(x >= 60),
    pct_above_60  = round(100 * sum(x >= 60) / length(x), 1)
  )
}

## Build one row per equation via `egfr_summary()` above (row-bound into a
## single table, one row per entry in `equations`, labelled by `equations`
## itself) and write it to an Excel file. Thin wrapper so callers don't need
## to inline the `dplyr::bind_rows(lapply(...))` + `openxlsx::write.xlsx()`
## themselves. Returns the combined table (invisibly).
write_egfr_summary_table <- function(cohort, equations, out_file) {
  egfr_table <- dplyr::bind_rows(lapply(equations, function(eq)
    egfr_summary(cohort[[eq]], name = eq)))
  
  openxlsx::write.xlsx(egfr_table, rowNames = FALSE, file = out_file)
  
  return(invisible(egfr_table))
}

## ---- Step: reclassification tables ------------------------------------------

## Reclassification of patients into risk categories, each equation vs. the
## reference (CKD-EPIcr 2009), for the full sample and split by outcome.
## Writes two Excel files (2-year and 5-year).
build_reclassification_tables <- function(cohort,
                                          equations,
                                          file_2y = "Reclassification 2-year predictions.xlsx",
                                          file_5y = "Reclassification 5-year predictions.xlsx") {
  # 2-year: binary <=40% vs. >40% risk; 5-year: three risk bands
  levels_2y <- c(">40%", "\u226440%")
  levels_5y <- c("<3%", "3% - 5%", ">5%")
  
  # restrict to the eGFR range where each horizon's prediction is clinically
  # relevant (same ranges used for the DCA restriction in run_dca())
  cohort_risk_2y_reclas <- cohort |> dplyr::filter(ckd_epi_2009_cr >= 10 &
                                                     ckd_epi_2009_cr <= 30)
  cohort_risk_5y_reclas <- cohort |> dplyr::filter(ckd_epi_2009_cr >= 30 &
                                                     ckd_epi_2009_cr <= 60)
  
  # categorize every equation's predicted risk into the levels above, adding
  # one "cat_risk_<horizon>y_<equation>" column per equation
  for (equation in equations) {
    risk_2y <- cohort_risk_2y_reclas[[paste0("risk_2y_", equation)]]
    cohort_risk_2y_reclas[[paste0("cat_risk_2y_", equation)]] <- ifelse(risk_2y > 0.40, ">40%", "\u226440%")
    
    risk_5y <- cohort_risk_5y_reclas[[paste0("risk_5y_", equation)]]
    cohort_risk_5y_reclas[[paste0("cat_risk_5y_", equation)]] <-
      ifelse(risk_5y < 0.03,
             "<3%",
             ifelse(risk_5y >= 0.03 &
                      risk_5y <= 0.05, "3% - 5%", ">5%"))
  }
  
  # split into "had the outcome" (cases, outcome == 1) vs. "did not" (control,
  # censored or competing event), so reclassification can be reported
  # separately for each group as well as for everyone combined
  cohort_risk_2y_reclas_cases <- subset(cohort_risk_2y_reclas, outcome_2y == 1)
  cohort_risk_2y_reclas_control <- subset(cohort_risk_2y_reclas, outcome_2y %in% c(0, 2))
  cohort_risk_5y_reclas_cases <- subset(cohort_risk_5y_reclas, outcome_5y == 1)
  cohort_risk_5y_reclas_control <- subset(cohort_risk_5y_reclas, outcome_5y %in% c(0, 2))
  
  table_2y <- c()
  table_5y <- c()
  
  # for every non-reference equation (equations[-1], since equations[1] IS
  # the reference), build 3 reclassification tables (all/cases/control) vs.
  # CKD-EPIcr 2009, and stack them into one combined table per horizon
  for (equation in equations[-1]) {
    table_2y_all <- generate_reclass_table(
      data = cohort_risk_2y_reclas,
      reference_col = "cat_risk_2y_ckd_epi_2009_cr",
      risk_categories = cohort_risk_2y_reclas[[paste0("cat_risk_2y_", equation)]],
      levels = levels_2y,
      title = paste("Reclassification Table all patients after 2 years", equation)
    )
    table_5y_all <- generate_reclass_table(
      data = cohort_risk_5y_reclas,
      reference_col = "cat_risk_5y_ckd_epi_2009_cr",
      risk_categories = cohort_risk_5y_reclas[[paste0("cat_risk_5y_", equation)]],
      levels = levels_5y,
      title = paste("Reclassification Table all patients after 5 years", equation)
    )
    
    table_2y_cases <- generate_reclass_table(
      data = cohort_risk_2y_reclas_cases,
      reference_col = "cat_risk_2y_ckd_epi_2009_cr",
      risk_categories = cohort_risk_2y_reclas_cases[[paste0("cat_risk_2y_", equation)]],
      levels = levels_2y,
      title = paste(
        "Reclassification Table patients with KFRT after 2 years",
        equation
      )
    )
    table_5y_cases <- generate_reclass_table(
      data = cohort_risk_5y_reclas_cases,
      reference_col = "cat_risk_5y_ckd_epi_2009_cr",
      risk_categories = cohort_risk_5y_reclas_cases[[paste0("cat_risk_5y_", equation)]],
      levels = levels_5y,
      title = paste(
        "Reclassification Table patients with KFRT after 5 years",
        equation
      )
    )
    
    table_2y_control <- generate_reclass_table(
      data = cohort_risk_2y_reclas_control,
      reference_col = "cat_risk_2y_ckd_epi_2009_cr",
      risk_categories = cohort_risk_2y_reclas_control[[paste0("cat_risk_2y_", equation)]],
      levels = levels_2y,
      title = paste(
        "Reclassification Table patients without KFRT after 2 years",
        equation
      )
    )
    table_5y_control <- generate_reclass_table(
      data = cohort_risk_5y_reclas_control,
      reference_col = "cat_risk_5y_ckd_epi_2009_cr",
      risk_categories = cohort_risk_5y_reclas_control[[paste0("cat_risk_5y_", equation)]],
      levels = levels_5y,
      title = paste(
        "Reclassification Table patients without KFRT after 5 years",
        equation
      )
    )
    
    # stack all/cases/control vertically (with a blank spacer row between
    # each), then append as a new block of columns for this equation --
    # so the final table reads left-to-right as one block per equation
    table_2y <- cbind(
      table_2y,
      rbind(
        table_2y_all$table,
        rep("", ncol(table_2y_all$table)),
        table_2y_cases$table,
        rep("", ncol(table_2y_all$table)),
        table_2y_control$table
      )
    )
    table_5y <- cbind(
      table_5y,
      rbind(
        table_5y_all$table,
        rep("", ncol(table_5y_all$table)),
        table_5y_cases$table,
        rep("", ncol(table_5y_all$table)),
        table_5y_control$table
      )
    )
  }
  
  # Prepend a single "Category" column, built directly from the known row
  # layout (all-patients block, blank spacer, cases block, blank spacer,
  # control block -- each block `length(levels)` rows tall) rather than
  # trying to carry a category column through generate_reclass_table() and
  # the rbind()/cbind() loop above. This is deterministic and easy to
  # verify by inspection: its length is guaranteed to equal nrow(table_2y).
  category_col_2y <- c(levels_2y, "", levels_2y, "", levels_2y)
  category_col_5y <- c(levels_5y, "", levels_5y, "", levels_5y)
  
  stopifnot(
    "category_col_2y doesn't match table_2y's row count" = length(category_col_2y) == nrow(table_2y),
    "category_col_5y doesn't match table_5y's row count" = length(category_col_5y) == nrow(table_5y)
  )
  
  table_2y <- cbind(Category = category_col_2y, table_2y)
  table_5y <- cbind(Category = category_col_5y, table_5y)
  
  openxlsx::write.xlsx(table_2y, file = file_2y, rowNames = FALSE)
  openxlsx::write.xlsx(table_5y, file = file_5y, rowNames = FALSE)
  
  return(invisible(list(
    table_2y = table_2y, table_5y = table_5y
  )))
}

## ---- Step: sensitivity analyses ---------------------------------------------

## Fit a subgroup-specific cumulative-incidence curve (KFRT vs. death, as
## competing risks), used by run_sensitivity_core() below to compute a
## subgroup's O/E ratio. compute_performance_measures() (the main table) can
## reuse the full-cohort `cin_2y`/`cin_5y` objects loaded from
## cohort_predictions.RData, but a sensitivity subgroup's O/E has to be
## measured against a CIF fit on that SAME restricted set of patients, so it
## is refit here from `data`'s horizon-capped time/status columns (i.e.
## `data` is the subgroup's own `sub_data`, already selected/recoded by the
## caller -- same convention as score_predictions()/evaluate_predictions()
## elsewhere in this file).
##
## Uses the same event coding as everywhere else in this file: fstatus = 0
## censored, 1 = KFRT (the cause of interest), 2 = death (competing risk).
## cmprsk::cuminc() with no `group` argument returns one curve per cause,
## named "1 1", "1 2", ... in increasing cause order, and oe_ratio() (further
## up) always reads $est[1] as the cause-of-interest estimate -- exactly as
## it does for the pre-fit `cin_2y`/`cin_5y` objects -- so the causes must
## sort with 1 (KFRT) first, matching how those objects are built.
fit_subgroup_cif <- function(data, horizon) {
  cmprsk::cuminc(
    ftime = data[[paste0("time_to_event_", horizon, "y")]],
    fstatus = data[[paste0("outcome_", horizon, "y")]],
    cencode = 0
  )
}

## Column names for one measurement "block" (AUC, calibration intercept &
## slope, O/E ratio, Brier, scaled Brier -- all with 95% CIs) for a given
## horizon and subgroup `suffix` (e.g. "omit1" or "subgroup2"). Used by
## run_sensitivity_core() below to keep every sensitivity analysis's output
## table laid out consistently.
##
## Thin wrapper around `measure_colnames()` (defined above, alongside
## `compute_performance_measures()`) -- that function is the single source
## of truth for this column layout, shared between the main results table
## and every sensitivity table, so the two can never drift out of sync.
##
## This mirrors compute_performance_measures()'s main table in full: O/E is
## included (using a per-subgroup CIF from fit_subgroup_cif() above, rather
## than the full-cohort cin_2y/cin_5y), and Brier/scaled-Brier get bootstrap
## CIs too, via run_sensitivity_core() calling evaluate_predictions()
## instead of score_predictions(). Note this makes the sensitivity analyses
## considerably slower (a shared bootstrap now runs for EVERY equation x
## horizon x subgroup combination -- see `bootstrap_brier_oe()`) -- lower
## `B` in run_sensitivity_core() if this becomes a bottleneck.
sa_measure_colnames <- function(horizon, suffix) {
  measure_colnames(horizon, suffix)
}

## Write one equation's AUC/calibration/O-E/Brier measures (for one horizon
## and subgroup `suffix`) into `results_df` at `row`, using the column names
## from `sa_measure_colnames()` above. `performance` is the list returned by
## `evaluate_predictions()` (AUC/calibration point estimates + analytic CIs,
## plus bootstrap CIs for O/E and Brier/scaled Brier). Returns the updated
## data frame (data frames aren't modified in place in R, so call sites need
## `results_df <- write_sa_measures(...)`).
write_sa_measures <- function(results_df,
                              row,
                              horizon,
                              suffix,
                              performance) {
  results_df[row, paste0("AUC_", horizon, "y_", suffix)] <-
    fmt_est_ci(performance$AUC, performance$AUC_CI)
  results_df[row, paste0("Int_", horizon, "y_", suffix)] <-
    fmt_est_ci(performance$Intercept, performance$Intercept_CI)
  results_df[row, paste0("Slope_", horizon, "y_", suffix)] <-
    fmt_est_ci(performance$Slope, performance$Slope_CI)
  results_df[row, paste0("OE_", horizon, "y_", suffix)] <-
    fmt_est_ci(performance$O_E_Ratio, performance$O_E_CI)
  results_df[row, paste0("Brier_", horizon, "y_", suffix)] <-
    fmt_est_ci(performance$Brier_Score, performance$Brier_CI)
  results_df[row, paste0("Scaled_Brier_", horizon, "y_", suffix)] <-
    fmt_est_ci(
      performance$Scaled_Brier_Score,
      performance$Scaled_Brier_CI,
      digits = 1,
      scale = 100
    )
  
  return(results_df)
}

## Shared engine behind every run_sensitivity_*() function below. Handles
## the part that was previously duplicated three times almost verbatim:
## building the results table's column layout, looping over horizon x
## subgroup x equation, fitting the subgroup's CIF and calling
## evaluate_predictions(), recording the subgroup's N once, writing each
## equation's measures, and saving the Excel file. Callers only need to say
## what the subgroups ARE.
##
## `subgroups` is a named list -- one entry per subgroup, name = the suffix
## used in column names (e.g. "omit1", "subgroup2", "egfr_range") -- where
## each entry is a FUNCTION OF `horizon` that returns list(select_patients
## = <logical vector over cohort's rows>, data = <data frame passed to
## evaluate_predictions(), usually cohort[select_patients, ], optionally with
## columns recoded>). Taking a function of `horizon` (rather than a fixed
## list) lets a subgroup's definition depend on the horizon (e.g. outcome
## comparison's "everyone except status X BY THIS HORIZON") while a
## subgroup that doesn't depend on horizon can just ignore the argument and
## return the same precomputed list every time (e.g. the eGFR ones below).
##
## `pred_risks` is read from `cohort` (not from a subgroup's own `data`)
## since the risk columns themselves are never recoded, only the outcome
## columns some subgroups touch.
##
## `B` is the number of bootstrap replicates used for the Brier/scaled-Brier
## and O/E confidence intervals (passed straight through to
## `evaluate_predictions()`/`bootstrap_brier_oe()`; `B = 0` skips the
## bootstrap entirely, same convention as `compute_performance_measures()`).
## Since this loop calls `evaluate_predictions()` once per equation x
## horizon x subgroup, and each call runs ONE shared B-replicate bootstrap
## (Brier, scaled Brier, and O/E all computed from the same resample -- see
## `bootstrap_brier_oe()`), the total number of bootstrap resamples run is
## `B * length(equations) * length(horizons) * length(subgroups)` -- lower
## `B` (or pass `B = 0` for point estimates
## only) if a sensitivity analysis becomes too slow.
run_sensitivity_core <- function(cohort,
                                 horizons,
                                 equations,
                                 model_names,
                                 out_file,
                                 subgroups,
                                 B = 500,
                                 seed = 123) {
  suffixes <- names(subgroups)
  
  # build the full column-name vector up front (one measurement block per
  # horizon x suffix combination) so the output table has a stable, known
  # column order regardless of loop iteration order
  col_names <- unlist(lapply(horizons, function(h) {
    unlist(lapply(suffixes, function(s)
      sa_measure_colnames(h, s)))
  }))
  
  # one extra row ("N") on top of one row per equation
  results_df <- data.frame(matrix(nrow = length(equations) + 1, ncol = length(col_names)))
  rownames(results_df) <- c("N", equations)
  colnames(results_df) <- col_names
  
  for (horizon in horizons) {
    for (suffix in suffixes) {
      # resolved once per (horizon, suffix) -- not once per equation, since
      # the subgroup itself doesn't depend on which equation is being
      # scored
      subgroup <- subgroups[[suffix]](horizon)
      select_patients <- subgroup$select_patients
      sub_data <- subgroup$data
      
      # subgroup-specific CIF (KFRT vs. death), refit on this subgroup's own
      # patients via fit_subgroup_cif() above, rather than reusing the
      # full-cohort cin_2y/cin_5y -- computed once per (horizon, suffix),
      # not once per equation, since it only depends on the outcome, not on
      # which equation is being scored.
      subgroup_CIF <- fit_subgroup_cif(sub_data, horizon)
      
      for (equation in equations) {
        pred_risks <- cohort[[paste0("risk_", horizon, "y_", equation)]][select_patients]
        
        # ALL measures (AUC, calibration intercept/slope, O/E ratio, Brier,
        # scaled Brier), with CIs -- via evaluate_predictions(), the same
        # function compute_performance_measures() uses for the main table.
        # time_auc/status_auc/time_restricted/status_restricted are derived
        # from `sub_data` inside score_predictions() (called by
        # validate_predictions() inside evaluate_predictions()) itself.
        performance <- evaluate_predictions(
          pred_risks = pred_risks,
          data = sub_data,
          horizon = horizon,
          CIF = subgroup_CIF,
          B = B,
          seed = seed
        )
        
        # record the subgroup size once per (horizon, suffix) -- it's the
        # same across equations, so only compute it on the first pass. Its
        # own dedicated "N_<horizon>y_<suffix>" column (from
        # sa_measure_colnames() above) rather than piggybacked onto
        # "AUC_<horizon>y_<suffix>", so it's not stuck inside the AUC tab
        # only -- write_measures_xlsx() surfaces N_ columns in every tab.
        if (equation == equations[1]) {
          results_df["N", paste0("N_", horizon, "y_", suffix)] <- sum(select_patients)
        }
        results_df <- write_sa_measures(results_df, equation, horizon, suffix, performance)
      }
    }
  }
  
  rownames(results_df) <- c("N", model_names)
  write_measures_xlsx(results_df, file = out_file)
  
  return(invisible(results_df))
}

## Discrimination (AUC), calibration (intercept & slope), and Brier / scaled
## Brier recomputed after omitting one of the two competing outcomes
## (status_num 1 = KFRT, 2 = death) at a time. Writes an Excel file.
run_sensitivity_outcome_comparison <- function(cohort,
                                               horizons,
                                               equations,
                                               model_names,
                                               out_file,
                                               B = 500,
                                               seed = 123) {
  # keep everyone EXCEPT those with the omitted outcome BY the horizon
  # (subgroup definition stays horizon-specific -- hence a function of
  # `horizon`, per run_sensitivity_core()'s `subgroups` contract). AUC uses
  # the uncapped _inf columns, Brier/calibration use the horizon-capped
  # columns -- see comment in validate_predictions() above.
  make_omit_subgroup <- function(status_num) {
    function(horizon) {
      select_patients <- cohort[[paste0("outcome_", horizon, "y")]] != status_num
      sub_data <- cohort[select_patients, ]
      
      # when omitting death (status_num == 1), the remaining "outcome 2"
      # values (death) can't occur among the selected patients, but recode
      # KFRT (already coded 1) defensively in case any 2s slipped through
      # -- applied to both time bases
      if (status_num == 1) {
        sub_data$outcome_inf[sub_data$outcome_inf == 2] <- 1
        sub_data[[paste0("outcome_", horizon, "y")]][sub_data[[paste0("outcome_", horizon, "y")]] == 2] <- 1
      }
      
      list(select_patients = select_patients, data = sub_data)
    }
  }
  
  # status_num 1 = drop patients who had KFRT (recompute measures on death
  # vs. censored only); status_num 2 = drop patients who died (recompute
  # measures on KFRT vs. censored only) -- i.e. omit each competing outcome
  # in turn
  subgroups <- list(omit1 = make_omit_subgroup(1), omit2 = make_omit_subgroup(2))
  
  run_sensitivity_core(cohort, horizons, equations, model_names, out_file, subgroups, B = B, seed = seed)
}

## Insert `suffix` right before the file extension of `file` (e.g.
## add_filename_suffix("SA_eGFR.xlsx", "ckd_epi_2009_cr") ->
## "SA_eGFR_ckd_epi_2009_cr.xlsx"). Used by run_sensitivity_egfr_subgroup()
## below to give each creatinine-eGFR comparator its own output file.
add_filename_suffix <- function(file, suffix) {
  sub("(\\.[^.]+)$", paste0("_", suffix, "\\1"), file)
}

## Discrimination (AUC), calibration (intercept & slope), and Brier / scaled
## Brier recomputed within subgroups defined by the relative difference
## between eGFRcys (2012) and a creatinine-based eGFR equation -- run once
## per comparator in `cr_equations` (default: the 2009 and 2021 CKD-EPI
## creatinine equations), since either one is a defensible choice for the
## "creatinine-based eGFR" side of the disagreement and the subgroup
## definition (and everything downstream) shifts depending on which is
## used. Each comparator gets its own Excel file, named via
## `add_filename_suffix()` above (e.g. "SA_eGFR.xlsx" ->
## "SA_eGFR_ckd_epi_2009_cr.xlsx", "SA_eGFR_ckd_epi_2021_cr.xlsx"). Returns
## a named list of each comparator's results data frame (invisibly), named
## by `cr_equations`.
run_sensitivity_egfr_subgroup <- function(cohort,
                                          horizons,
                                          equations,
                                          model_names,
                                          out_file,
                                          cr_equations = c("ckd_epi_2009_cr", "ckd_epi_2021_cr"),
                                          B = 500,
                                          seed = 123) {
  results_list <- list()
  
  for (cr_equation in cr_equations) {
    results_list[[cr_equation]] <- run_sensitivity_egfr_subgroup_one(
      cohort = cohort,
      horizons = horizons,
      equations = equations,
      model_names = model_names,
      cr_equation = cr_equation,
      out_file = add_filename_suffix(out_file, cr_equation),
      B = B,
      seed = seed
    )
  }
  
  return(invisible(results_list))
}

## One creatinine-eGFR comparator's worth of the eGFR-agreement subgroup
## sensitivity analysis -- the body of run_sensitivity_egfr_subgroup()
## above, factored out so it can be run once per entry in `cr_equations`.
## `cr_equation` is the column name of the creatinine-based eGFR to compare
## cystatin-C eGFR against (e.g. "ckd_epi_2009_cr").
run_sensitivity_egfr_subgroup_one <- function(cohort,
                                              horizons,
                                              equations,
                                              model_names,
                                              cr_equation,
                                              out_file,
                                              B = 500,
                                              seed = 123) {
  # relative difference between the cystatin-C-based eGFR and this
  # creatinine-based comparator -- large positive/negative values flag
  # patients where the two markers disagree substantially
  cohort$eGFRdiff <- (cohort$ckd_epi_2012_cys - cohort[[cr_equation]]) / cohort[[cr_equation]]
  
  # classify into 3 subgroups by that relative difference: cystatin-C eGFR
  # notably lower, roughly equal, or notably higher than creatinine eGFR
  # (using +/-20% as the cutoff)
  eGFRsubgroup_labels <- c("eGFRcys < eGFRcr", "eGFRcys ~ eGFRcr", "eGFRcys > eGFRcr")
  cohort$eGFR_subgroup <- factor(
    dplyr::case_when(cohort$eGFRdiff <= -0.2 ~ 1, cohort$eGFRdiff >= 0.2  ~ 3, TRUE                    ~ 2),
    levels = 1:3,
    labels = eGFRsubgroup_labels
  )
  
  # this subgroup definition doesn't depend on horizon, so each function
  # just ignores the `horizon` argument run_sensitivity_core() passes it
  # and returns the same precomputed selection every time
  make_egfr_subgroup <- function(label) {
    select_patients <- cohort$eGFR_subgroup == label
    sub_data <- cohort[select_patients, ]
    function(horizon)
      list(select_patients = select_patients, data = sub_data)
  }
  
  subgroups <- stats::setNames(lapply(eGFRsubgroup_labels, make_egfr_subgroup),
                               paste0("subgroup", seq_along(eGFRsubgroup_labels)))
  
  run_sensitivity_core(cohort, horizons, equations, model_names, out_file, subgroups, B = B, seed = seed)
}

## Discrimination (AUC), calibration (intercept & slope), and Brier / scaled
## Brier recomputed after restricting to patients whose eGFR falls within
## [egfr_min, egfr_max] on EVERY eGFR equation in `equations` at once (not
## just the one being evaluated) -- i.e. the same restricted patient set is
## used for every equation's measures, so equations are still compared on a
## common population, just one narrowed to the clinically relevant eGFR
## range where all equations agree on inclusion. A patient with a missing
## eGFR for any equation is excluded (treated as failing the range check
## for that equation, not silently passed through). Writes a single Excel
## file (a single subgroup, unlike run_sensitivity_egfr_subgroup() above).
run_sensitivity_egfr_range <- function(cohort,
                                       horizons,
                                       equations,
                                       model_names,
                                       out_file,
                                       egfr_min = 10,
                                       egfr_max = 60,
                                       B = 500,
                                       seed = 123) {
  # TRUE for a patient only if every equation's eGFR is non-missing and
  # within [egfr_min, egfr_max]; Reduce(`&`, ...) requires ALL equations to
  # pass, not just one. Doesn't depend on horizon, so computed once here
  # rather than inside the subgroup function below.
  in_range_by_equation <- lapply(equations, function(eq) {
    eGFR <- cohort[[eq]]
    ! is.na(eGFR) & eGFR >= egfr_min & eGFR <= egfr_max
  })
  select_patients <- Reduce(`&`, in_range_by_equation)
  sub_data <- cohort[select_patients, ]
  
  subgroups <- list(
    egfr_range = function(horizon)
      list(select_patients = select_patients, data = sub_data)
  )
  
  run_sensitivity_core(cohort, horizons, equations, model_names, out_file, subgroups, B = B, seed = seed)
}