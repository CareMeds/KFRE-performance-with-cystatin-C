################################################################################
## This file contains functions to perform analyses ############################
## Author: Malou Magnani #######################################################
################################################################################

## Compute the linear predictor (PI), which is the same for all time horizons
PI_KFRE <- function(age, male, egfr, alb) {
  PI <- -0.2201 *
    (age / 10 - 7.036) +
    0.2467 * (male - 0.5642) -
    0.5567 * (egfr / 5 - 7.222) +
    0.4510 * (log(alb) - 5.137)

  return(PI)
}

## Compute 2y KFRE risk (P)
KFRE_risk_2y <- function(PI) {
  risk_2y <- 1 - (0.9832**exp(PI))

  return(risk_2y)
}

# Compute 5y KFRE risk (P)
KFRE_risk_5y <- function(PI) {
  risk_5y <- 1 - (0.9365**exp(PI))

  return(risk_5y)
}

## Compute C-statistic and CIs with bootstrapping
c_statistic <- function(
  PI, # prognostic index
  data, # cohort data
  horizon, # time horizon
  bootstrap = FALSE, # to calculate CIs
  B = 500,
  seed = 123
) {
  # Calculate C-statistic
  discr <- Hmisc::rcorrcens(
    formula = formula(paste0(
      "survival::Surv(data$time_to_event_",
      horizon,
      "y_modified, data$outcome_",
      horizon,
      "y_modified) ~ I(-PI)"
    ))
  )
  c_stat <- discr[1, 1]

  # Initialize results with just the C-statistic
  results <- list(
    C_Statistic = c_stat
  )

  # Add bootstrapping if requested
  if (bootstrap) {
    # Define bootstrap function
    boot_func <- function(data, indices, risk_var) {
      # Create resampled dataset
      resampled_data <- data[indices, ]
      # Apply the risk variable to resampled data
      resampled_risk <- risk_var[indices]
      # Calculate C-statistic on resampled data
      boot_discr <- Hmisc::rcorrcens(
        formula = formula(paste0(
          "survival::Surv(resampled_data$time_to_event_",
          horizon,
          "y_modified, resampled_data$outcome_",
          horizon,
          "y_modified) ~ I(-resampled_risk)"
        ))
      )
      boot_c_stat <- boot_discr[1, 1]
      return(as.numeric(boot_c_stat))
    }
    # Set seed for reproducibility
    set.seed(seed)
    # Run bootstrap
    boot_results <- boot::boot(
      data = data,
      statistic = boot_func,
      R = B,
      risk_var = PI
    )
    # Calculate bootstrap confidence intervals
    c_stat_ci <- boot::boot.ci(boot_results, type = "perc")
    # Add only the CI to the results
    results$C_Statistic_CI <- if (!is.null(c_stat_ci)) {
      c(Lower = c_stat_ci$percent[4], Upper = c_stat_ci$percent[5])
    } else {
      NULL
    }
  }

  return(results)
}

## Compute O/E ratio with CIs using bootstrapping
oe_ratio <- function(
  pred_risks,
  CIF,
  horizon,
  bootstrap = FALSE,
  B = 500,
  seed = 123
) {
  # Calculate original O/E ratio
  expected <- mean(pred_risks)
  observed <- cmprsk::timepoints(CIF, horizon * 365.25)
  oe <- observed$est[1] / expected

  # Initialize results with just the O/E ratio
  results <- list(
    O_E_Ratio = oe
  )

  # Add bootstrapping if requested
  if (bootstrap) {
    # Define bootstrap function
    boot_func <- function(data, indices, pred_risks, cif_obj, time_horizon) {
      # Resample risk predictions
      resampled_risks <- pred_risks[indices]

      # Calculate expected risk in bootstrap sample
      boot_expected <- mean(resampled_risks)

      # For observed, we need to recalculate CIF with resampled data
      # This depends on how your CIF object is created
      # This is a placeholder - you'll need to adjust based on your data structure
      boot_observed <- cmprsk::timepoints(cif_obj, time_horizon)

      # Calculate O/E ratio
      boot_oe <- boot_observed$est[1] / boot_expected

      return(as.numeric(boot_oe))
    }

    # Set seed for reproducibility
    set.seed(seed)

    # Create a dummy dataset for bootstrapping
    # We need this because boot() requires a data argument
    dummy_data <- data.frame(id = 1:length(pred_risks))

    # Run bootstrap
    boot_results <- boot::boot(
      data = dummy_data,
      statistic = boot_func,
      R = B,
      pred_risks = pred_risks,
      cif_obj = CIF,
      time_horizon = horizon * 365.25
    )

    # Calculate bootstrap confidence intervals
    oe_ci <- boot::boot.ci(boot_results, type = "perc")

    # Add only the CI to the results
    results$O_E_CI <- if (!is.null(oe_ci)) {
      c(Lower = oe_ci$percent[4], Upper = oe_ci$percent[5])
    } else {
      NULL
    }
  }

  return(results)
}

## Compute calibration intercept & slope with 95% CIs
cal_int_slope <- function(pred_risks, data, horizon) {
  # extract outcome variables according to horizon
  data$time_to_event <- eval(parse(
    text = paste0("data$time_to_event_", horizon, "y")
  ))
  data$outcome <- eval(parse(text = paste0("data$outcome_", horizon, "y")))

  Score <- riskRegression::Score(
    # Give score list of predictions for all individuals
    # (This will be used for the predicted part of validation)
    list(pred_risks),
    # Define the model as subdistributional hazards model
    formula = Hist(time_to_event, outcome) ~ 1,
    # specify pseudovalues to be used
    cens.method = "pseudo",
    # define validation dataset with values for covariates,
    # this is to be used for observed probabilities
    data = data,
    # define prediction horizon
    times = horizon * 365.25,
    # define event of interest
    outcome = 1,
    conf.int = TRUE,
    # define validation methods
    plots = "calibration",
    summary = c("ipa"),
    metrics = "brier"
  )
  
  # extract unsmoothed pseudo-observations 
  pseudos <- data.frame(Score$Calibration$plotframe)
  
  # add the cloglog risk estimates
  pseudos$cll_pred <- log(-log(1 - pseudos$risk))
  
  # fit model for calibration intercept
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
  
  # fit model for calibration slope
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
  
  return(list(
    Intercept = Intercept,
    Intercept_CI = c(Intercept - qnorm(0.975) * Intercept_SE,
                     Intercept + qnorm(0.975) * Intercept_SE),
    Slope = 1 + Slope_est,
    Slope_CI = c(1 + (Slope_est - qnorm(0.975) * Slope_SE),
                 1 + (Slope_est + qnorm(0.975) * Slope_SE))
  ))
}

## Compute Scaled Brier & Brier score, with CIs using bootstrapping
brier_scores <- function(
  pred_risks,
  data,
  horizon,
  bootstrap = FALSE,
  B = 500,
  seed = 123
) {
  # extract outcome variables according to horizon
  data$time_to_event <- eval(parse(
    text = paste0("data$time_to_event_", horizon, "y")
  ))
  data$outcome <- eval(parse(text = paste0("data$outcome_", horizon, "y")))

  # Original Brier score calculation
  brier <- riskRegression::Score(
    list(pred_risks),
    formula = Hist(time_to_event, outcome) ~ 1,
    cens.method = "pseudo",
    data = data,
    times = horizon * 365.25,
    outcome = 1,
    conf.int = TRUE,
    metrics = "brier",
    summary = "ipa"
  )

  # Extract values for the numeric model only
  numeric_model <- brier$Brier$score[model == "numeric"]

  # Store original values
  brier_score <- numeric_model$Brier
  ipa_score <- numeric_model$IPA

  # Initialize results with just the scores
  results <- list(
    Brier_Score = brier_score,
    Scaled_Brier_Score = ipa_score
  )

  # Add bootstrapping if requested
  if (bootstrap) {
    # Define bootstrap function
    boot_func <- function(data, indices, risk_var) {
      resampled_data <- data[indices, ]
      resampled_risk <- risk_var[indices]

      boot_brier <- riskRegression::Score(
        list(resampled_risk),
        formula = Hist(time_to_event, outcome) ~ 1,
        cens.method = "pseudo",
        data = resampled_data,
        conf.int = TRUE,
        times = horizon * 365.25,
        outcome = 1,
        metrics = "brier",
        summary = "ipa"
      )

      ipa_value <- boot_brier$Brier$score$IPA[
        boot_brier$Brier$score$model == "numeric"
      ]
      brier_value <- boot_brier$Brier$score$Brier[
        boot_brier$Brier$score$model == "numeric"
      ]

      return(c(brier = as.numeric(brier_value), ipa = as.numeric(ipa_value)))
    }

    # Set seed for reproducibility
    set.seed(seed)

    # Run bootstrap
    boot_results <- boot::boot(
      data = data,
      statistic = boot_func,
      R = B,
      risk_var = pred_risks
    )

    # Calculate bootstrap confidence intervals
    brier_ci <- boot::boot.ci(boot_results, type = "perc", index = 1)
    ipa_ci <- boot::boot.ci(boot_results, type = "perc", index = 2)

    # Add only the CIs to the results
    results$Brier_CI <- if (!is.null(brier_ci)) {
      c(Lower = brier_ci$percent[4], Upper = brier_ci$percent[5])
    } else {
      NULL
    }

    results$Scaled_Brier_CI <- if (!is.null(ipa_ci)) {
      c(Lower = ipa_ci$percent[4], Upper = ipa_ci$percent[5])
    } else {
      NULL
    }
  }

  return(results)
}

# Create 2-year reclassification tables
generate_reclass_table_2y <- function(data, risk_categories, title) {
  # Ensure all levels are present in the data
  data$cat_risk_2y_ckd_epi_2009_cr <- factor(
    data$cat_risk_2y_ckd_epi_2009_cr,
    levels = c(">40%", "≤40%")
  )
  risk_categories <- factor(risk_categories, levels = c(">40%", "≤40%"))
  # Create a table of reclassification
  table_raw <- table(data$cat_risk_2y_ckd_epi_2009_cr, risk_categories)
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

# Create 5-year reclassification tables
generate_reclass_table_5y <- function(data, risk_categories, title) {
  # Ensure all levels are present in the data
  data$cat_risk_5y_ckd_epi_2009_cr <- factor(
    data$cat_risk_5y_ckd_epi_2009_cr,
    levels = c("<3%", "3% - 5%", ">5%")
  )
  risk_categories <- factor(
    risk_categories,
    levels = c("<3%", "3% - 5%", ">5%")
  )
  # Create a table of reclassification
  table_raw <- table(data$cat_risk_5y_ckd_epi_2009_cr, risk_categories)
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
