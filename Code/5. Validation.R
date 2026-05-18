################################################################################
## This script is intended to validate predictions #############################
## Author: Malou Magnani #######################################################
################################################################################
# remove history
rm(list = ls(all.names = TRUE))

# set seed for reproducibility
set.seed(27)

# set directory to save results
setwd("P:/SCREAM2/SCREAM2_Research/Malou Magnani/Final/Results new/")

# load data sets
load("P:/SCREAM2/SCREAM2_Research/Malou Magnani/Final/Data/cohort_predictions.RData")

# load functions
source("P:/SCREAM2/SCREAM2_Research/Malou Magnani/Final/Code/Functions for analyses.R")

# load libraries
library(survival)       # time-to-event analyses
library(patchwork)      # combine plots

################################################################################
### For each equation, calculate performance measures ##########################
################################################################################
horizons <- c(2, 5)
measures <- c(
  "AUC_2y",
  "AUC_2y_CI",
  "AUC_5y",
  "AUC_5y_CI",
  "Int_2y",
  "Int_2y_CI",
  "Int_5y",
  "Int_5y_CI",
  "Slope_2y",
  "Slope_2y_CI",
  "Slope_5y",
  "Slope_5y_CI",
  "OE_2y",
  "OE_2y_CI",
  "OE_5y",
  "OE_5y_CI",
  "Brier_2y",
  "Brier_2y_CI",
  "Brier_5y",
  "Brier_5y_CI",
  "Scaled_Brier_2y",
  "Scaled_Brier_2y_CI",
  "Scaled_Brier_5y",
  "Scaled_Brier_5y_CI"
)
equations <- c(
  "ckd_epi_2009_cr",
  "ckd_epi_2021_cr",
  "ckd_epi_2012_cys",
  "ckd_epi_2012_cr_cys",
  "ckd_epi_2021_cr_cys"
)
model_names <- c(
  "CKD-EPIcr 2009",
  "CKD-EPIcr 2021",
  "CKD-EPIcys 2012",
  "CKD-EPIcrcys 2012",
  "CKD-EPIcrcys 2021"
)
results_df <- data.frame(matrix(nrow = length(equations), ncol = length(measures)))
rownames(results_df) <- equations
colnames(results_df) <- measures
n_bootstraps <- 500
for (horizon in horizons) {
  for (equation in equations) {
    # Prognostic index
    PI <- eval(parse(text = paste0("cohort$PI_", equation)))
    time <- eval(parse(text = paste0(
      "cohort$time_to_event_", horizon, "y"
    )))
    status <- eval(parse(text = paste0("cohort$outcome_", horizon, "y")))
    
    # Time-dependent AUC by Blanche
    discrimination <- timeROC::timeROC(
      T = time,
      delta = status,
      marker = PI,
      times = horizon * 365,
      cause = 1,
      iid = FALSE
    )
    results_df[equation, paste0("AUC_", horizon, "y")] <-
      sprintf("%.3f", discrimination$AUC_1[2])
    
    # calibration intercept and slope
    pred_risks <- eval(parse(text = paste0(
      "cohort$risk_", horizon, "y_", equation
    )))
    int_slope <- cal_int_slope(pred_risks = pred_risks,
                               data = cohort,
                               horizon = horizon)
    results_df[equation, paste0("Int_", horizon, "y")] <-
      sprintf("%.3f", int_slope$Intercept)
    results_df[equation, paste0("Int_", horizon, "y_CI")] <-
      paste0(
        "[",
        sprintf("%.3f", int_slope$Intercept_CI[1]),
        "; ",
        sprintf("%.3f", int_slope$Intercept_CI[2]),
        "]"
      )
    results_df[equation, paste0("Slope_", horizon, "y")] <-
      sprintf("%.3f", int_slope$Slope)
    results_df[equation, paste0("Slope_", horizon, "y_CI")] <-
      paste0(
        "[",
        sprintf("%.3f", int_slope$Slope_CI[1]),
        "; ",
        sprintf("%.3f", int_slope$Slope_CI[2]),
        "]"
      )
    
    # O/E ratio
    CIF <- eval(parse(text = paste0("cin_", horizon, "y")))
    OE <- oe_ratio(
      pred_risks = pred_risks,
      CIF = CIF,
      horizon = horizon,
      B = n_bootstraps
    )
    results_df[equation, paste0("OE_", horizon, "y")] <-
      sprintf("%.3f", OE$O_E_Ratio)
    
    # brier score
    brier <- brier_scores(
      pred_risks = pred_risks,
      data = cohort,
      horizon = horizon,
      B = n_bootstraps
    )
    results_df[equation, paste0("Brier_", horizon, "y")] <-
      sprintf("%.3f", brier$Brier_Score)
    results_df[equation, paste0("Scaled_Brier_", horizon, "y")] <-
      sprintf("%.1f", brier$Scaled_Brier_Score * 100)
    
    # add CI from bootstrapping if turned on
    if (n_bootstraps > 0) {
      # discrimination
      # discrimination <- timeROC::timeROC(T=time,
      #                                    delta=status,
      #                                    marker=PI,
      #                                    times=horizon*365,
      #                                    cause=1,
      #                                    iid=TRUE)
      # discrimination_CI <- stats::confint(discrimination)
      # results_df[equation, paste0("AUC_", horizon, "y_CI")] <-
      #   paste0("[", sprintf("%.3f", discrimination_CI$CI_AUC_1[1]/100), "; ",
      #          sprintf("%.3f", discrimination_CI$CI_AUC_1[2]/100), "]")
      
      # O/E ratio
      results_df[equation, paste0("OE_", horizon, "y_CI")] <-
        paste0("[",
               sprintf("%.3f", OE$O_E_CI["Lower"]),
               "; ",
               sprintf("%.3f", OE$O_E_CI["Upper"]),
               "]")
      
      # Brier score
      results_df[equation, paste0("Brier_", horizon, "y_CI")] <-
        paste0(
          "[",
          sprintf("%.3f", brier$Brier_CI["Lower"]),
          "; ",
          sprintf("%.3f", brier$Brier_CI["Upper"]),
          "]"
        )
      results_df[equation, paste0("Scaled_Brier_", horizon, "y_CI")] <-
        paste0(
          "[",
          sprintf("%.1f", brier$Scaled_Brier_CI["Lower"] * 100),
          "; ",
          sprintf("%.1f", brier$Scaled_Brier_CI["Upper"] * 100),
          "]"
        )
    }
  }
}
# save to Excel
rownames(results_df) <- model_names
openxlsx::write.xlsx(results_df, rowNames = TRUE, file = "Measures.xlsx")

################################################################################
### Histograms of predicted risks
################################################################################
colors_palette <- c("darkorange4",
                    "darkred",
                    "darkorchid4",
                    "darkblue",
                    "darkgreen")
quantile_cut <- function(x)
  quantile(x, 0.95, na.rm = TRUE)
max_5y <- max(sapply(cohort[, "risk_5y_ckd_epi_2009_cr"], quantile_cut))

horizons_long <- rep(rep(horizons, each = length(equations)), 2)
trim_long <- rep(c(TRUE, FALSE), each = length(equations) * 2)

histogram_risk <- Map(function(equation,
                               horizon,
                               xlab,
                               trim = TRUE,
                               max_5y,
                               color) {
  equation <- paste0("risk_", horizon, "y_", equation)
  hist <- ggplot2::ggplot(cohort, ggplot2::aes(x = .data[[equation]])) +
    ggplot2::geom_histogram(
      binwidth = ifelse(
        horizon == 2 & !trim,
        0.025,
        ifelse(
          horizon == 2 & trim,
          0.005,
          ifelse(horizon == 5 & !trim, 0.025, ifelse(horizon == 5 &
                                                       trim, 0.01))
        )
      ),
      boundary = 0,
      fill = color,
      color = "white"
    ) +
    ggthemes::theme_clean() +
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
    hist <- hist + ggplot2::coord_cartesian(xlim = c(0, max_5y))
  }
  return(hist)
},
equations,
horizons_long,
model_names,
trim_long,
max_5y,
rep(colors_palette, length(horizons_long) / length(colors_palette)))

names(histogram_risk) <- paste0("risk_",
                                horizons_long,
                                "y_",
                                equations,
                                ifelse(trim_long, "_trimmed", ""))

################################################################################
### Combined calibration plot for all models for both horizons #################
################################################################################
pred_list_2y <- list(
  cohort$risk_2y_ckd_epi_2009_cr,
  cohort$risk_2y_ckd_epi_2021_cr,
  cohort$risk_2y_ckd_epi_2012_cys,
  cohort$risk_2y_ckd_epi_2012_cr_cys,
  cohort$risk_2y_ckd_epi_2021_cr_cys
)
pred_list_5y <- list(
  cohort$risk_5y_ckd_epi_2009_cr,
  cohort$risk_5y_ckd_epi_2021_cr,
  cohort$risk_5y_ckd_epi_2012_cys,
  cohort$risk_5y_ckd_epi_2012_cr_cys,
  cohort$risk_5y_ckd_epi_2021_cr_cys
)
for (trim in c(TRUE, FALSE)) {
  for (horizon in horizons) {
    # Initialize calibration data
    calibration_data <- data.frame()
    
    # Loop over models
    pred_list <- eval(parse(text = paste0("pred_list_", horizon, "y")))
    for (nr_equation in 1:length(pred_list)) {
      Score <- riskRegression::Score(
        list("model" = pred_list[[nr_equation]]),
        formula = eval(parse(
          text = paste0(
            "Hist(time_to_event_",
            horizon,
            "y, outcome_",
            horizon,
            "y) ~ 1"
          )
        )),
        cens.method = "pseudo",
        data = cohort,
        times = horizon * 365.25,
        outcome = 1,
        conf.int = TRUE,
        plots = "calibration"
      )
      
      # Extract and smooth pseudo-values
      pseudos <- data.frame(Score$Calibration$plotframe) |>
        dplyr::arrange(risk)
      smooth_pseudos <- predict(stats::loess(
        pseudovalue ~ risk,
        data = pseudos,
        degree = 1,
        span = 0.33
      ))
      
      # Store in a data frame
      temp_df <- data.frame(risk = pseudos$risk,
                            observed = smooth_pseudos,
                            model = model_names[nr_equation])
      
      # Append to calibration_data
      calibration_data <- dplyr::bind_rows(calibration_data, temp_df)
    }
    
    # reorder models
    calibration_data$model <- factor(calibration_data$model, levels = model_names)
    
    # create combined calibration plot
    combined_cal_plot <- ggplot2::ggplot(calibration_data,
                                         ggplot2::aes(x = risk,
                                                      y = observed, 
                                                      color = model)) +
      ggplot2::geom_line(linewidth = 0.75, alpha = 0.8) +
      ggplot2::scale_color_manual(values = colors_palette,
                                  breaks = model_names) +
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
        ggplot2::coord_cartesian(xlim = c(0, max_5y), ylim = c(0, max_5y))
    } else{
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
        ggplot2::coord_cartesian(xlim = c(0, max_5y), ylim = c(0, max_5y))
    }
    assign(paste0(
      "combined_cal_plot_",
      horizon,
      ifelse(trim, "_trimmed", "")
    ),
    combined_cal_plot)
  }
}

# position two plots side by side
combined_cal_plots <-
  (combined_cal_plot_2 + combined_cal_plot_5) +
  patchwork::plot_layout(guides = "collect") &
  ggplot2::theme(legend.position = "bottom")
combined_hist_plots <-
  histogram_risk$risk_2y_ckd_epi_2009_cr +
  histogram_risk$risk_5y_ckd_epi_2009_cr

# save plot
ggplot2::ggsave(
  filename = paste("Combined calibration plot.png"),
  plot = (combined_cal_plots / combined_hist_plots) +
    patchwork::plot_layout(heights = c(3, 1)),
  width = 9,
  height = 6,
  dpi = 300
)

# position two plots side by side
combined_cal_plots_trimmed <-
  (combined_cal_plot_2_trimmed + combined_cal_plot_5_trimmed) +
  patchwork::plot_layout(guides = "collect") &
  ggplot2::theme(legend.position = "bottom")
combined_hist_plots_trimmed <-
  histogram_risk$risk_2y_ckd_epi_2009_cr_trimmed +
  histogram_risk$risk_5y_ckd_epi_2009_cr_trimmed

# save plot
ggplot2::ggsave(
  filename = paste("Combined calibration plot trimmed.png"),
  plot = (combined_cal_plots_trimmed / combined_hist_plots_trimmed) +
    patchwork::plot_layout(heights = c(3, 1)),
  width = 9,
  height = 6,
  dpi = 300
)

################################################################################
### Distribution plots 2-years #################################################
################################################################################
equation_order_2y <- c(
  "risk_2y_ckd_epi_2021_cr_cys",
  "risk_2y_ckd_epi_2012_cr_cys",
  "risk_2y_ckd_epi_2012_cys",
  "risk_2y_ckd_epi_2021_cr",
  "risk_2y_ckd_epi_2009_cr"
)

# long format
cohort_risk_2y_long <- cohort |>
  dplyr::select(
    lopnr,
    risk_2y_ckd_epi_2009_cr,
    risk_2y_ckd_epi_2021_cr,
    risk_2y_ckd_epi_2012_cys,
    risk_2y_ckd_epi_2012_cr_cys,
    risk_2y_ckd_epi_2021_cr_cys
  ) |>
  tidyr::pivot_longer(cols = starts_with("risk_"),
                      names_to = "Equation",
                      values_to = "Risk") |>
  dplyr::mutate(
    Equation = factor(Equation, levels = equation_order_2y),
    Equation = dplyr::recode(
      Equation,
      "risk_2y_ckd_epi_2009_cr" = "CKD-EPIcr 2009",
      "risk_2y_ckd_epi_2021_cr" = "CKD-EPIcr 2021",
      "risk_2y_ckd_epi_2012_cys" = "CKD-EPIcys 2012",
      "risk_2y_ckd_epi_2012_cr_cys" = "CKD-EPIcrcys 2012",
      "risk_2y_ckd_epi_2021_cr_cys" = "CKD-EPIcrcys 2021"
    )
  )

# 0-100% density ridge plot
plot_KFRE_2y <- ggplot2::ggplot(cohort_risk_2y_long,
                                ggplot2::aes(x = Risk, y = Equation, fill = stat(x))) +
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
  ggplot2::labs(title = "2-year KFRE risk 0-100%",
                x = "Predicted risk",
                y = "Equation") +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    panel.spacing = ggplot2::unit(0.5, "lines"),
    axis.text.y = ggplot2::element_text(size = 10)
  )

# 0-0.5% Create density ridge plot
plot_KFRE_zoomed_00.5_2y <- ggplot2::ggplot(cohort_risk_2y_long,
                                            ggplot2::aes(x = Risk, y = Equation, fill = stat(x))) +
  ggridges::geom_density_ridges_gradient(scale = 0.9,
                                         rel_min_height = 0.000001,
                                         alpha = 0.8) +
  ggplot2::scale_fill_viridis_c(
    name = "Risk",
    option = "C",
    limits = c(0, 0.004),
    oob = scales::squish
  ) +
  ggplot2::scale_x_continuous(limits = c(0, 0.005),
                              breaks = seq(0, 0.005, by = 0.001)) +
  ggplot2::scale_y_discrete(expand = c(0, 0)) +
  ggplot2::labs(title = "2-year KFRE risk 0-0.5%", x = "Risk Score", y = "Equation") +
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

# Combine plots
plot_KFRE_combined_2y <- plot_KFRE_2y +
  ggplot2::annotation_custom(
    ggplot2::ggplotGrob(plot_KFRE_zoomed_00.5_2y),
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
    plot.background = ggplot2::element_rect(fill = "transparent", color = "transparent"),
    panel.background = ggplot2::element_rect(fill = "white", color = "transparent")
  )

################################################################################
### Distribution plots 2-years #################################################
################################################################################
equation_order_5y <- c(
  "risk_5y_ckd_epi_2021_cr_cys",
  "risk_5y_ckd_epi_2012_cr_cys",
  "risk_5y_ckd_epi_2012_cys",
  "risk_5y_ckd_epi_2021_cr",
  "risk_5y_ckd_epi_2009_cr"
)

# long format
cohort_risk_5y_long <- cohort |>
  dplyr::select(
    lopnr,
    risk_5y_ckd_epi_2009_cr,
    risk_5y_ckd_epi_2021_cr,
    risk_5y_ckd_epi_2012_cys,
    risk_5y_ckd_epi_2012_cr_cys,
    risk_5y_ckd_epi_2021_cr_cys
  ) |>
  tidyr::pivot_longer(cols = starts_with("risk_"),
                      names_to = "Equation",
                      values_to = "Risk") |>
  dplyr::mutate(
    Equation = factor(Equation, levels = equation_order_5y),
    Equation = dplyr::recode(
      Equation,
      "risk_5y_ckd_epi_2009_cr" = "CKD-EPIcr 2009",
      "risk_5y_ckd_epi_2021_cr" = "CKD-EPIcr 2021",
      "risk_5y_ckd_epi_2012_cys" = "CKD-EPIcys 2012",
      "risk_5y_ckd_epi_2012_cr_cys" = "CKD-EPIcrcys 2012",
      "risk_5y_ckd_epi_2021_cr_cys" = "CKD-EPIcrcys 2021"
    )
  )

# 0-100% density ridge plot
plot_KFRE_5y <- ggplot2::ggplot(cohort_risk_5y_long,
                                ggplot2::aes(x = Risk, y = Equation, fill = stat(x))) +
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
  ggplot2::labs(title = "5-year KFRE risk 0-100%", 
                x = "Risk Score", 
                y = "Equation") +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    panel.spacing = ggplot2::unit(0.5, "lines"),
    axis.text.y = ggplot2::element_text(size = 10)
  )

# 0-1% Create density ridge plot
plot_KFRE_zoomed_01_5y <- ggplot2::ggplot(cohort_risk_5y_long,
                                          ggplot2::aes(x = Risk, 
                                                       y = Equation, 
                                                       fill = stat(x))) +
  ggridges::geom_density_ridges_gradient(scale = 0.9,
                                         rel_min_height = 0.000001,
                                         alpha = 0.8) +
  ggplot2::scale_fill_viridis_c(
    name = "Risk",
    option = "C",
    limits = c(0, 0.008),
    oob = scales::squish
  ) +
  ggplot2::scale_x_continuous(limits = c(0, 0.01),
                              breaks = seq(0, 0.01, by = 0.002)) +
  ggplot2::scale_y_discrete(expand = c(0, 0)) +  # uniform spacing
  ggplot2::labs(title = "5y KFRE Risk Distributions by Equation 0-1%", 
                x = "Risk Score", 
                y = "Equation") +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    plot.background = ggplot2::element_rect(fill = ggplot2::alpha("white", 0.8), 
                                            color = NA),
    panel.background = ggplot2::element_rect(fill = ggplot2::alpha("white", 0.8), 
                                             color = NA),
    legend.title = ggplot2::element_text(size = 8),
    legend.text = ggplot2::element_text(size = 6),
    legend.key.size = ggplot2::unit(0.5, "lines"),
    panel.spacing = ggplot2::unit(0.5, "lines"),
    axis.text.y = ggplot2::element_text(size = 10),
    axis.text.x = ggplot2::element_text(size = 10),
    axis.title.x = ggplot2::element_text(size = 8),
    axis.title.y = ggplot2::element_text(size = 8),
    plot.title = ggplot2::element_text(size = 8)
  )

# Combine plots
plot_KFRE_combined_5y <- plot_KFRE_5y +
  ggplot2::annotation_custom(
    ggplot2::ggplotGrob(plot_KFRE_zoomed_01_5y),
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

# combine plots for 2 and 5 years prediction horizon
plot_KFRE_combined_2y <- plot_KFRE_combined_2y + ggplot2::theme(legend.position = "bottom")
plot_KFRE_combined_5y <- plot_KFRE_combined_5y + ggplot2::theme(legend.position = "bottom")
combined_dist_plot <- plot_KFRE_combined_2y + plot_KFRE_combined_5y +
  patchwork::plot_layout(ncol = 1, guides = "collect")
ggplot2::ggsave(
  filename = "Combined distribution plots.png",
  plot = combined_dist_plot,
  width = 15,
  height = 10,
  dpi = 300
)

################################################################################
### Risk differences ###########################################################
################################################################################
### 2-year KFRE
cohort_risk_diff_2y_long <- cohort_risk_2y_long |>
  dplyr::group_by(lopnr) |>
  dplyr::mutate(risk_diff = Risk - Risk[Equation == "CKD-EPIcr 2009"]) |>
  dplyr::ungroup() |>
  dplyr::filter(Equation != "CKD-EPIcr 2009")

# -0.5 to 0.5 %
plot_KFRE_diff_2y <- ggplot2::ggplot(cohort_risk_diff_2y_long,
                                     ggplot2::aes(x = risk_diff, 
                                                  y = Equation, 
                                                  fill = stat(x))) +
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
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "black"
  ) +
  ggplot2::labs(title = "2y KFRE risk difference compared to CKD-EPIcr 2009",
                x = "Risk difference (Risk equation - Risk CKD-EPIcr 2009)",
                y = "Equation") +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    panel.spacing = ggplot2::unit(0.5, "lines"),
    axis.text.y = ggplot2::element_text(size = 10)
  )

### 5-year KFRE
cohort_risk_diff_5y_long <- cohort_risk_5y_long |>
  dplyr::group_by(lopnr) |>
  dplyr::mutate(risk_diff = Risk - Risk[Equation == "CKD-EPIcr 2009"]) |>
  dplyr::ungroup() |>
  dplyr::filter(Equation != "CKD-EPIcr 2009")

# -1 to 1 %
plot_KFRE_diff_5y <- ggplot2::ggplot(cohort_risk_diff_5y_long,
                                     ggplot2::aes(x = risk_diff, 
                                                  y = Equation, 
                                                  fill = stat(x))) +
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
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "black"
  ) +
  ggplot2::labs(title = "5y KFRE Risk Differences Compared to CKD-EPIcr 2009", 
                x = "Risk Difference (Risk Equation - Risk CKD-EPIcr 2009)",
                y = "Equation") +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    panel.spacing = ggplot2::unit(0.5, "lines"),
    axis.text.y = ggplot2::element_text(size = 10)
  )

# combine plots
combined_plot_diff <- plot_KFRE_diff_2y + plot_KFRE_diff_5y +
  plot_layout(ncol = 1, guides = "collect")
ggplot2::ggsave(
  filename = "Combined risk difference.png",
  plot = combined_plot_diff,
  width = 15,
  height = 10,
  dpi = 300
)

################################################################################
### eGFR distributions #########################################################
################################################################################
# long format
cohort_egfr_long <- cohort |>
  tidyr::pivot_longer(
    cols = c(
      ckd_epi_2021_cr_cys,
      ckd_epi_2012_cr_cys,
      ckd_epi_2012_cys,
      ckd_epi_2021_cr,
      ckd_epi_2009_cr
    ),
    names_to = "Equation",
    values_to = "eGFR"
  ) |>
  dplyr::mutate(
    Equation = factor(Equation, levels = equations),
    Equation = dplyr::recode(
      Equation,
      "ckd_epi_2009_cr" = "CKD-EPIcr 2009",
      "ckd_epi_2021_cr" = "CKD-EPIcr 2021",
      "ckd_epi_2012_cys" = "CKD-EPIcys 2012",
      "ckd_epi_2012_cr_cys" = "CKD-EPIcrcys 2012",
      "ckd_epi_2021_cr_cys" = "CKD-EPIcrcys 2021"
    )
  )

# create plot
egfr_plot <- ggplot2::ggplot(cohort_egfr_long,
                             ggplot2::aes(x = eGFR, y = Equation, fill = stat(x))) +
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
  ggplot2::scale_y_discrete(expand = c(0, 0)) +
  ggplot2::labs(title = "eGFR distributions", 
                x = "eGFR (mL/min/1.73m²)",
                y = "Equation") +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    panel.spacing = ggplot2::unit(1.5, "lines"),
    axis.text.y = ggplot2::element_text(size = 10)
  ) +
  ggplot2::theme(plot.margin = ggplot2::margin(20, 10, 20, 10))
ggplot2::ggsave(
  filename = "eGFR distribution.png",
  plot = egfr_plot,
  width = 15,
  height = 10,
  dpi = 300
)

################################################################################
### Reclassification tables ####################################################
################################################################################
# 2y eGFR <30
cohort_risk_2y_reclas <- cohort |>
  dplyr::filter(ckd_epi_2009_cr >= 10 & ckd_epi_2009_cr <= 30)

# 5y eGFR 30-60
cohort_risk_5y_reclas <- cohort |>
  dplyr::filter(ckd_epi_2009_cr >= 30 & ckd_epi_2009_cr <= 60)

# categorize 2-year predicted risk below and above 40%
# categorize 5-year predicted risk in three categories
for (equation in equations) {
  risk <- eval(parse(text = paste0(
    "cohort_risk_2y_reclas$risk_2y_", equation
  )))
  cohort_risk_2y_reclas[, paste0("cat_risk_2y_", equation)] <-
    ifelse(risk > 0.40, ">40%", "≤40%")
  
  risk <- eval(parse(text = paste0(
    "cohort_risk_5y_reclas$risk_5y_", equation
  )))
  cohort_risk_5y_reclas[, paste0("cat_risk_5y_", equation)] <-
    ifelse(risk < 0.03,
           "<3%",
           ifelse(risk >= 0.03 & risk <= 0.05, "3% - 5%", ">5%"))
}

# Split data based on outcome
cohort_risk_2y_reclas_cases <- subset(cohort_risk_2y_reclas, outcome_2y == 1)
cohort_risk_2y_reclas_control <- subset(cohort_risk_2y_reclas, outcome_2y %in% c(0, 2))
cohort_risk_5y_reclas_cases <- subset(cohort_risk_5y_reclas, outcome_5y == 1)
cohort_risk_5y_reclas_control <- subset(cohort_risk_5y_reclas, outcome_5y %in% c(0, 2))

# Make two reclassification tables for 2-year and 5-year predictions
table_2y <- c()
table_5y <- c()
for (equation in equations[-1]) {
  # for all patients
  table_2y_all <- generate_reclass_table_2y(
    data = cohort_risk_2y_reclas,
    risk_categories = eval(parse(
      text = paste0("cohort_risk_2y_reclas$cat_risk_2y_", equation)
    )),
    title = paste("Reclassification Table all patients after 2 years", equation)
  )
  table_5y_all <- generate_reclass_table_5y(
    data = cohort_risk_5y_reclas,
    risk_categories = eval(parse(
      text = paste0("cohort_risk_5y_reclas$cat_risk_5y_", equation)
    )),
    title = paste("Reclassification Table all patients after 5 years", equation)
  )
  
  # for patients with the outcome
  table_2y_cases <- generate_reclass_table_2y(
    data = cohort_risk_2y_reclas_cases,
    risk_categories = eval(parse(
      text = paste0("cohort_risk_2y_reclas_cases$cat_risk_2y_", equation)
    )),
    title = paste(
      "Reclassification Table patients with KFRT after 2 years",
      equation
    )
  )
  table_5y_cases <- generate_reclass_table_5y(
    data = cohort_risk_5y_reclas_cases,
    risk_categories = eval(parse(
      text = paste0("cohort_risk_5y_reclas_cases$cat_risk_5y_", equation)
    )),
    title = paste(
      "Reclassification Table patients with KFRT after 5 years",
      equation
    )
  )
  
  # for patients without the outcome
  table_2y_control <- generate_reclass_table_2y(
    data = cohort_risk_2y_reclas_control,
    risk_categories = eval(parse(
      text = paste0("cohort_risk_2y_reclas_control$cat_risk_2y_", equation)
    )),
    title = paste(
      "Reclassification Table patients without KFRT after 2 years",
      equation
    )
  )
  table_5y_control <- generate_reclass_table_5y(
    data = cohort_risk_5y_reclas_control,
    risk_categories = eval(parse(
      text = paste0("cohort_risk_5y_reclas_control$cat_risk_5y_", equation)
    )),
    title = paste(
      "Reclassification Table patients without KFRT after 5 years",
      equation
    )
  )
  
  # final table
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
# write to Excel
openxlsx::write.xlsx(table_2y, file = "Reclassification 2-year predictions.xlsx")
openxlsx::write.xlsx(table_5y, file = "Reclassification 5-year predictions.xlsx")

################################################################################
### DCA analysis
################################################################################
cohort <-
  cohort |>
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

new_DCA_run <- FALSE
if (new_DCA_run) {
  for (horizon in horizons) {
    cat("Time horizon", horizon, "years\n")
    if (horizon == 2) {
      max_tp <- 0.5
    } else{
      max_tp <- 0.2
    }
    
    for (nr_equation in 1:5) {
      cat("Equation", equations[nr_equation], "\n")
      # Create formula dynamically
      formula <- as.formula(
        paste0(
          "survival::Surv(time = time_to_event_",
          horizon,
          "y,
        event = outcome_",
          horizon,
          "y_dca) ~ risk_",
          horizon,
          "y_",
          equations[nr_equation]
        )
      )
      
      # restrict cohort
      if (horizon == 2) {
        res_cohort <- cohort[cohort$ckd_epi_2021_cr >= 10 &
                               cohort$ckd_epi_2021_cr < 30, ]
      } else if (horizon == 5) {
        res_cohort <- cohort[cohort$ckd_epi_2021_cr >= 30 &
                               cohort$ckd_epi_2021_cr < 60, ]
      }
      
      # Run DCA
      dca <- dcurves::dca(
        formula,
        data = res_cohort,
        thresholds = seq(0, max_tp, 0.01),
        time = horizon * 365
      )
      assign(paste0("dca_", horizon, "y_", nr_equation), dca)
      
      # Extract table
      table <- dplyr::as_tibble(dca$dca)
      assign(paste0("table_dca_", horizon, "y_", nr_equation), table)
    }
  }
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
    file = "DCA_analysis.Rdata"
  )
} else{
  load("DCA_analysis.Rdata")
}

for (horizon in horizons) {
  # set maximum threshold, vertical lines, and minimum net benefit on x-axis
  if (horizon == 2) {
    max_tp <- max_5y
    vlines <- c(0.1, 0.4)
    min_nb <- -0.02
    max_nb <- 0.14
  } else{
    max_tp <- max_5y
    vlines <- c(0.03, 0.05)
    min_nb <- -0.005
    max_nb <- 0.016
  }
  
  for (nr_equation in 1:5) {
    # extract table
    table_dca <- eval(parse(text = paste0(
      "table_dca_", horizon, "y_", nr_equation
    )))
    
    if (nr_equation == 1) {
      # extract net benefit for treat none
      treat_none <- table_dca |>
        dplyr::filter(label == "Treat None") |>
        dplyr::select(threshold, net_benefit)
      
      # extract net benefit for treat all
      treat_all <- table_dca |>
        dplyr::filter(label == "Treat All") |>
        dplyr::select(net_benefit)
      
      # attach to dt
      NB <- data.frame(
        thresholds = treat_none$threshold,
        treat_none = treat_none$net_benefit,
        treat_all = treat_all$net_benefit
      )
    }
    
    # extract net benefit of model
    NB_dca <- table_dca |>
      dplyr::filter(label == paste0("risk_", horizon, "y_", equations[nr_equation])) |>
      dplyr::select(net_benefit)
    
    # attach to dt
    NB[, equations[nr_equation]] <- NB_dca$net_benefit
  }
  # convert to data frame
  assign(paste0("NB_", horizon, "y"), NB)
  
  # convert to long format
  NB_long <- NB |>
    tidyr::pivot_longer(cols = -thresholds,
                        names_to = "variable",
                        values_to = "value")
  
  # custom breaks for x-axis and y-axis labels
  max_x <- max(max_tp, max(vlines, na.rm = TRUE))
  default_breaks <- scales::extended_breaks()(c(0, max_x))
  custom_x_breaks <- sort(unique(c(default_breaks, vlines)))
  custom_y_breaks <- seq(min(min_nb), max_nb, ifelse(horizon == 2, 0.02, 0.005))
  
  # Generate plot
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
        "–",
        ifelse(horizon == 2, 29, 59),
        " mL/min/1.73m² (N=",
        eval(parse(
          text = paste0("unique(table_dca_", horizon, "y_1$n)")
        )),
        ")"
      ),
      x = "Threshold Probability",
      y = "Net Benefit"
    ) +
    ggplot2::theme(
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 10),
      legend.margin = ggplot2::margin(t = 0, r = 0, b = 0, l = -20), 
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
    ggplot2::guides(color = ggplot2::guide_legend(
      # breaks = c("Treat_all", "Treat_none", equations),
      ncol = 7,
      byrow = TRUE
    )) +
    ggplot2::geom_vline(xintercept = vlines[1],
                        linetype = "dashed",
                        color = "black") +
    ggplot2::geom_vline(xintercept = vlines[2],
                        linetype = "dashed",
                        color = "black")
  assign(paste0("dca_plot_", horizon, "y"), dca_plot)
}

# save DCA plot
combined_dca_plot <- (dca_plot_2y + dca_plot_5y) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(
    theme = ggplot2::theme(legend.position = "bottom")
  )

# save plot
ggplot2::ggsave(
  filename = "Combined DCA plot.png",
  plot = combined_dca_plot,
  width = 10,
  height = 4,
  dpi = 300
)

# save table
NB_table <- rbind(data.frame(time = 2, NB_2y[NB_2y$thresholds %in% c(0.1, 0.4), ]), data.frame(time = 5, NB_5y[NB_5y$thresholds %in% c(0.03, 0.05), ]))
openxlsx::write.xlsx(NB_table, rowNames = FALSE, file = "Net Benefit Table.xlsx")

################################################################################
### Time-dependent AUC
################################################################################
td_AUC <- openxlsx::read.xlsx("P:/SCREAM2/SCREAM2_Research/Malou Magnani/Final/To Run/Time dependent AUC.xlsx")
colnames(td_AUC)[1] <- "Equation"

# Parse CI strings "[lo; hi]" into numeric columns
parse_ci <- function(x) {
  x <- gsub("\\[|\\]", "", x)
  parts <- strsplit(x, ";\\s*")
  data.frame(lo = as.numeric(sapply(parts, `[`, 1)), hi = as.numeric(sapply(parts, `[`, 2)))
}

ci_2y <- parse_ci(td_AUC$AUC_2y_CI)
ci_5y <- parse_ci(td_AUC$AUC_5y_CI)

td_AUC <- td_AUC |>
  dplyr::mutate(
    AUC_2y = as.numeric(AUC_2y),
    AUC_5y = as.numeric(AUC_5y),
    lo_2y  = ci_2y$lo,
    hi_2y  = ci_2y$hi,
    lo_5y  = ci_5y$lo,
    hi_5y  = ci_5y$hi,
    lab_2y = paste(sprintf("%.3f", AUC_2y), AUC_2y_CI),
    lab_5y = paste(sprintf("%.3f", AUC_5y), AUC_5y_CI),
    row_num  = dplyr::row_number(),
    shade    = row_num %% 2 == 0,
    Equation = factor(Equation, levels = rev(Equation))
  )

# Pivot to long format
td_AUC_long <- td_AUC |>
  tidyr::pivot_longer(
    cols      = c(AUC_2y, AUC_5y),
    names_to  = "horizon",
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

# Plot settings
x_limits <- c(0.91, 1.00)
x_breaks <- seq(0.92, 1.00, by = 0.02)
x_right  <- 1.005   # where the AUC text labels start (right side)
x_left   <- 0.91    # where the equation name labels are placed (left side, inside plot)
n_eq     <- nlevels(td_AUC_long$Equation)

# Shading: one rect per shaded row, spanning full x range
# Built from the long data so it applies to both facets identically
shade_df <- td_AUC_long |>
  dplyr::filter(shade) |>
  dplyr::mutate(y_pos = as.integer(Equation),
                ymin  = y_pos - 0.5,
                ymax  = y_pos + 0.5) |>
  dplyr::distinct(horizon, Equation, .keep_all = TRUE)

# Plot
AUC_plot <- ggplot2::ggplot(td_AUC_long, ggplot2::aes(
  y = Equation,
  x = AUC,
  xmin = lo,
  xmax = hi
)) +
  
  # Shading — xmin/xmax = -Inf/Inf covers the full panel width;
  # with clip = "off" and the label geom placed at x_left inside the scale,
  # the grey band will also sit behind the equation name text.
  ggplot2::geom_rect(
    data = shade_df,
    ggplot2::aes(
      ymin = ymin,
      ymax = ymax,
      xmin = -Inf,
      xmax = Inf
    ),
    inherit.aes = FALSE,
    fill  = "grey",
    color = NA
  ) +
  
  # Dashed reference line at x = 1.00
  ggplot2::geom_vline(
    xintercept = x_limits[2],
    color = "black",
    linewidth = 0.4,
    linetype = "dashed"
  ) +
  
  # Confidence intervals and point estimates
  ggplot2::geom_errorbarh(height = 0.15,
                          linewidth = 0.4,
                          color = "black") +
  ggplot2::geom_point(shape = 15,
                      size = 2,
                      color = "black") +
  
  # AUC labels on the right (outside panel, clip = "off" allows this)
  ggplot2::geom_text(
    ggplot2::aes(x = x_right, label = label),
    hjust = 0,
    size = 2.8,
    color = "black"
  ) +
  
  # Equation name labels placed INSIDE the panel at the far left.
  # Because they live inside the panel, the shading geom_rect (which also
  # spans the full panel) renders behind them automatically.
  ggplot2::geom_text(
    ggplot2::aes(x = x_left, y = Equation, label = Equation),
    hjust = 0,
    size = 2.8,
    color = "black",
    inherit.aes = FALSE
  ) +
  
  ggplot2::scale_x_continuous(
    limits = c(x_limits[1], 1.04),
    labels = c("", x_breaks[-1]),
    breaks = c(x_breaks[2], x_breaks[-1]),
    expand = c(0, 0)
  ) +
  ggplot2::coord_cartesian(
    clip = "off",
    xlim = c(x_limits[1], 1.04),
    ylim = c(0.5, n_eq + 0.5)
  ) +
  ggplot2::facet_wrap( ~ horizon, ncol = 1) +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(
    # Hide the real y-axis entirely — labels are now drawn by geom_text
    axis.line.y        = ggplot2::element_blank(),
    axis.ticks.y       = ggplot2::element_blank(),
    axis.title.y       = ggplot2::element_blank(),
    axis.text.y        = ggplot2::element_blank(),
    axis.title.x       = ggplot2::element_blank(),
    axis.line.x        = ggplot2::element_blank(),
    axis.text.x        = ggplot2::element_text(color = "black"),
    panel.grid         = ggplot2::element_blank(),
    legend.position    = "none",
    plot.margin        = ggplot2::margin(20, 0, 5, 5),
    strip.background   = ggplot2::element_blank(),
    strip.text         = ggplot2::element_text(
      face = "bold",
      size = 10,
      hjust = 1,
      color = "black"
    )
  ) +
  # Add a shorter x-axis line manually
  ggplot2::geom_segment(
    data = data.frame(
      x = x_breaks[2],
      xend = 1.00,
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

# Save
suppressMessages(ggplot2::ggsave(
  "forest_plot_AUC.png", AUC_plot, width = 5, height = 4, dpi = 300
))