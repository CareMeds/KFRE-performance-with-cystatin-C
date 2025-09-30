################################################################################
## This script is intended to validate predictions #############################
## Author: Malou Magnani #######################################################
################################################################################
# remove history
rm(list=ls(all.names=TRUE))

# set seed for reproducibility
set.seed(27)

# set directory to save results
setwd("~/Results/")

# load data sets
load("~/Data/cohort_predictions.RData")

# load functions
source("~/Code/Functions for analyses.R")

# load libraries
library(openxlsx)       # save to excel
library(riskRegression) # evalution metrics
library(Hmisc)          # C-index
library(survival)       # time-to-event analyses
library(cmprsk)         # competing risk
library(boot)           # bootstrapping
library(dplyr)          # data manipulation
library(tidyr)          # data manipulation
library(knitr)          # pretty tables
library(kableExtra)     # pretty tables
library(ggplot2)        # generate plots
library(ggthemes)       # themes of plots
library(ggridges)       # gradient distribution
library(patchwork)      # combine plots
library(cowplot)        # combine plots

################################################################################
### For each equation, calculate performance measures ##########################
################################################################################
horizons <- c(2, 5)
measures <- c("Cstat_5y", "Cstat_2y_CI", "Cstat_5y", "Cstat_5y_CI",
              "Int_2y", "Int_2y_CI", "Int_5y", "Int_5y_CI",
              "Slope_2y", "Slope_2y_CI", "Slope_5y", "Slope_5y_CI",
              "OE_2y", "OE_2y_CI", "OE_5y", "OE_5y_CI",
              "Brier_2y", "Brier_2y_CI", "Brier_5y", "Brier_5y_CI",
              "Scaled_Brier_2y", "Scaled_Brier_2y_CI", "Scaled_Brier_5y", 
              "Scaled_Brier_5y_CI")
equations <- c("ckd_epi_2009_cr", "ckd_epi_2021_cr", "ckd_epi_2012_cys",
               "ckd_epi_2012_cr_cys", "ckd_epi_2021_cr_cys")
model_names <- c("CKD-EPIcr 2009", 
                 "CKD-EPIcr 2021", 
                 "CKD-EPIcys 2012", 
                 "CKD-EPIcrcys 2012", 
                 "CKD-EPIcrcys 2021")
results_df <- data.frame(matrix(nrow=length(equations), ncol=length(measures)))
rownames(results_df) <- equations
colnames(results_df) <- measures
bootstrapping <- TRUE
for (horizon in horizons){
  for (equation in equations){
    # Prognostic index
    PI <- eval(parse(text=paste0("cohort$PI_", equation)))
    
    # Modified Harrell's C-statistic by Wolbers
    discrimination <- c_statistic(PI = PI,
                                  data = cohort,
                                  horizon = horizon,
                                  bootstrap = bootstrapping,
                                  B = 500)
    results_df[equation, paste0("Cstat_", horizon, "y")] <- 
      sprintf("%.3f", discrimination$C_Statistic)
    
    # calibration intercept and slope
    pred_risks <- eval(parse(text=paste0("cohort$risk_", 
                                         horizon, "y_", equation)))
    int_slope <- cal_int_slope(pred_risks = pred_risks, 
                               data = cohort, 
                               horizon = horizon)
    results_df[equation, paste0("Int_", horizon, "y")] <- 
      sprintf("%.3f", int_slope$Intercept)
    results_df[equation, paste0("Int_", horizon, "y_CI")] <-  
      paste0("[", sprintf("%.3f", int_slope$Intercept_CI["2.5 %"]),
             "; ", sprintf("%.3f", int_slope$Intercept_CI["97.5 %"]), "]")
    results_df[equation, paste0("Slope_", horizon, "y")] <- 
      sprintf("%.3f", int_slope$Slope)
    results_df[equation, paste0("Slope_", horizon, "y_CI")] <- 
      paste0("[", sprintf("%.3f", int_slope$Slope_CI["2.5 %"]),
             "; ", sprintf("%.3f", int_slope$Slope_CI["97.5 %"]), "]")
    
    # O/E ratio
    CIF <- eval(parse(text=paste0("cin_", horizon, "y")))
    OE <- oe_ratio(pred_risks = pred_risks, 
                   CIF = CIF, 
                   horizon = horizon, 
                   bootstrap = bootstrapping,
                   B = 500)
    results_df[equation, paste0("OE_", horizon, "y")] <- 
      sprintf("%.3f", OE$O_E_Ratio)
    
    # brier score
    brier <- brier_scores(pred_risks = pred_risks, 
                          data = cohort, 
                          horizon = horizon,
                          bootstrap = bootstrapping,
                          B = 500)
    results_df[equation, paste0("Brier_", horizon, "y")] <- 
      sprintf("%.3f", brier$Brier_Score)
    results_df[equation, paste0("Scaled_Brier_", horizon, "y")] <- 
      sprintf("%.1f", brier$Scaled_Brier_Score*100)
    
    # add CI from bootstrapping if turned on
    if (bootstrapping){
      results_df[equation, paste0("Cstat_", horizon, "y_CI")] <- 
        paste0("[", sprintf("%.3f", discrimination$C_Statistic_CI["Lower"]),
               "; ", sprintf("%.3f", discrimination$C_Statistic_CI["Upper"]), "]")
      results_df[equation, paste0("OE_", horizon, "y_CI")] <- 
        paste0("[", sprintf("%.3f", OE$O_E_CI["Lower"]),
               "; ", sprintf("%.3f", OE$O_E_CI["Upper"]), "]")
      results_df[equation, paste0("Brier_", horizon, "y_CI")] <- 
        paste0("[", sprintf("%.3f", brier$Brier_CI["Lower"]),
               "; ", sprintf("%.3f", brier$Brier_CI["Upper"]), "]")
      results_df[equation, paste0("Scaled_Brier_", horizon, "y_CI")] <- 
        paste0("[", sprintf("%.1f", brier$Scaled_Brier_CI["Lower"]*100),
               "; ", sprintf("%.1f", brier$Scaled_Brier_CI["Upper"]*100), "]")
    }
  }
}
# save to Excel
rownames(results_df) <- model_names
openxlsx::write.xlsx(results_df, rowNames=TRUE, file = "Measures.xlsx")

################################################################################
### Combined calibration plot for all models for both horizons #################
################################################################################
colors <- c("orange", 
            "skyblue", 
            "darkblue", 
            "darkred", 
            "darkgreen")
pred_list_2y <- list(cohort$risk_2y_ckd_epi_2009_cr,
                     cohort$risk_2y_ckd_epi_2021_cr,
                     cohort$risk_2y_ckd_epi_2012_cys,
                     cohort$risk_2y_ckd_epi_2012_cr_cys,
                     cohort$risk_2y_ckd_epi_2021_cr_cys)
pred_list_5y <- list(cohort$risk_5y_ckd_epi_2009_cr,
                     cohort$risk_5y_ckd_epi_2021_cr,
                     cohort$risk_5y_ckd_epi_2012_cys,
                     cohort$risk_5y_ckd_epi_2012_cr_cys,
                     cohort$risk_5y_ckd_epi_2021_cr_cys)
for (horizon in horizons){
  # Initialize calibration data
  calibration_data <- data.frame()
  
  # Loop over models
  pred_list <- eval(parse(text=paste0("pred_list_", horizon, "y")))
  for (i in 1:length(pred_list)) {
    Score <- riskRegression::Score(
      list("model" = pred_list[[i]]),
      formula = eval(parse(text=paste0("Hist(time_to_event_", horizon, 
                                       "y, outcome_", horizon, "y) ~ 1"))),
      cens.method = "pseudo",
      data = cohort,
      times = horizon*365.25,
      outcome = 1,
      conf.int = TRUE,
      plots = "calibration"
    )
    
    # Extract and smooth pseudo-values
    pseudos <- data.frame(Score$Calibration$plotframe) |>
      dplyr::arrange(risk)
    smooth_pseudos <- predict(
      stats::loess(pseudovalue ~ risk, data = pseudos, degree = 1, span = 0.33)
    )
    
    # Store in a data frame
    temp_df <- data.frame(
      risk = pseudos$risk,
      observed = smooth_pseudos,
      model = model_names[i]
    )
    
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
    ggplot2::annotate("segment", x = 0, y = 0, xend = 1, yend = 1, 
                      linetype = "dashed", color = "gray40") +
    ggplot2::scale_color_manual(values = colors, breaks = model_names) +
    ggplot2::scale_y_continuous(breaks = seq(0, 1, by = 0.25), 
                                limits = c(0, 1)) +
    ggplot2::scale_x_continuous(limits = c(0, 1)) +
    ggplot2::theme(
      legend.text = ggplot2::element_text(size = 8), 
      legend.title = ggplot2::element_text(size = 9)) +
    ggplot2::labs(
      title = paste0(horizon, "-year KFRE CKD-EPI equations"),
      x = "Predicted Risks",
      y = "Observed Risks",
      color = "Equation"
    ) +
    ggthemes::theme_clean() +
    ggplot2::theme(
      legend.text = ggplot2::element_text(size = 8), 
      legend.title = ggplot2::element_text(size = 9),
      legend.key.size = ggplot2::unit(0.5, "cm"))
  
  assign(paste0("combined_cal_plot_", horizon), combined_cal_plot)
}
# position two plots side by side
combined_cal_plots <- cowplot::plot_grid(
  combined_cal_plot_2, combined_cal_plot_5,
  ncol = 2, 
  align = "hv",            # align axes
  axis = "tblr",           # match all axes: top, bottom, left, right
  rel_heights = c(1, 1, 1) # adjust row heights if needed
)
# save plot
ggplot2::ggsave(filename = paste("Combined calibration plot.png"),
                plot = combined_cal_plots, width = 12, height = 4, dpi = 300)

################################################################################
### Distribution plots 2-years #################################################
################################################################################
equation_order_2y <- c("risk_2y_ckd_epi_2021_cr_cys", 
                       "risk_2y_ckd_epi_2012_cr_cys", 
                       "risk_2y_ckd_epi_2012_cys", 
                       "risk_2y_ckd_epi_2021_cr", 
                       "risk_2y_ckd_epi_2009_cr")

# long format
cohort_risk_2y_long <- cohort |> 
  dplyr::select(lopnr, 
                risk_2y_ckd_epi_2009_cr, 
                risk_2y_ckd_epi_2021_cr, 
                risk_2y_ckd_epi_2012_cys, 
                risk_2y_ckd_epi_2012_cr_cys, 
                risk_2y_ckd_epi_2021_cr_cys) |> 
  tidyr::pivot_longer(cols = starts_with("risk_"), 
                      names_to = "Equation", 
                      values_to = "Risk") |> 
  dplyr::mutate(Equation = factor(Equation, levels = equation_order_2y),
                Equation = dplyr::recode(Equation,
                                         "risk_2y_ckd_epi_2009_cr" = "CKD-EPIcr 2009",
                                         "risk_2y_ckd_epi_2021_cr" = "CKD-EPIcr 2021",
                                         "risk_2y_ckd_epi_2012_cys" = "CKD-EPIcys 2012",
                                         "risk_2y_ckd_epi_2012_cr_cys" = "CKD-EPIcrcys 2012",
                                         "risk_2y_ckd_epi_2021_cr_cys" = "CKD-EPIcrcys 2021"))

# 0-100% density ridge plot
plot_KFRE_2y <- ggplot2::ggplot(cohort_risk_2y_long, 
                                ggplot2::aes(x = Risk, 
                                             y = Equation, 
                                             fill = stat(x))) +
  ggridges::geom_density_ridges_gradient(scale = 1, 
                                         rel_min_height = 0.000001, 
                                         alpha = 0.8) +  
  ggplot2::scale_fill_viridis_c(name = "Risk", option = "C", 
                                limits = c(0, 0.05), 
                                oob = scales::squish) +  
  ggplot2::scale_x_continuous(limits = c(0, 1), 
                              breaks = seq(0, 1, by = 0.2)) +  
  ggplot2::scale_y_discrete(expand = c(0, 0)) + # uniform spacing
  ggplot2::labs(title = "2-year KFRE risk 0-100%",
                x = "Predicted risk", y = "Equation") +
  ggplot2::theme_minimal(base_size = 12) +  
  ggplot2::theme(panel.spacing = ggplot2::unit(0.5, "lines"), 
                 axis.text.y = ggplot2::element_text(size = 10))

# 0-0.5% Create density ridge plot
plot_KFRE_zoomed_00.5_2y <- ggplot2::ggplot(cohort_risk_2y_long, 
                                            aes(x = Risk,
                                                y = Equation, 
                                                fill = stat(x))) +
  ggridges::geom_density_ridges_gradient(scale = 0.9, 
                                         rel_min_height = 0.000001, 
                                         alpha = 0.8) +  
  ggplot2::scale_fill_viridis_c(name = "Risk", option = "C", 
                                limits = c(0, 0.004), 
                                oob = scales::squish) +  
  ggplot2::scale_x_continuous(limits = c(0, 0.005), 
                              breaks = seq(0, 0.005, by = 0.001)) +  
  ggplot2::scale_y_discrete(expand = c(0, 0)) +  # uniform spacing
  ggplot2::labs(title = "2-year KFRE risk 0-0.5%",
                x = "Risk Score", y = "Equation") +
  ggplot2::theme_minimal(base_size = 12) +  
  ggplot2::theme(
    plot.background = ggplot2::element_rect(fill = alpha("white", 0.8), 
                                            color = NA), 
    panel.background = ggplot2::element_rect(fill = alpha("white", 0.8), 
                                             color = NA),
    legend.title = ggplot2::element_text(size = 8), 
    legend.text = ggplot2::element_text(size = 6),
    legend.key.size = ggplot2::unit(0.5, "lines"),
    panel.spacing = ggplot2::unit(0.5, "lines"), 
    axis.text.y = ggplot2::element_text(size = 10), 
    axis.text.x = ggplot2::element_text(size = 8),  
    axis.title.x = ggplot2::element_text(size = 8), 
    axis.title.y = ggplot2::element_text(size = 8), 
    plot.title = ggplot2::element_text(size = 8))   

# Combine plots
plot_KFRE_combined_2y <- plot_KFRE_2y + 
  ggplot2::annotation_custom(ggplot2::ggplotGrob(plot_KFRE_zoomed_00.5_2y), 
                             xmin = 0.4, xmax = 1.05, # Horizontal placement of the zoomed plot
                             ymin = 3, ymax = 6) + 
  ggplot2::theme_classic() +
  ggplot2::theme(legend.position = "bottom",
                 legend.title = ggplot2::element_text(size = 8),  
                 legend.text = ggplot2::element_text(size = 6),   
                 plot.tag = ggplot2::element_blank(), # Remove tag
                 axis.text.x = ggplot2::element_text(size = 10),
                 axis.text = ggplot2::element_text(size = 6),
                 axis.text.y = ggplot2::element_text(size = 10),
                 plot.background = ggplot2::element_rect(fill = "transparent", 
                                                         color = "transparent"),
                 panel.background = ggplot2::element_rect(fill = "white", 
                                                          color = "transparent"))

################################################################################
### Distribution plots 2-years #################################################
################################################################################
equation_order_5y <- c("risk_5y_ckd_epi_2021_cr_cys", 
                       "risk_5y_ckd_epi_2012_cr_cys", 
                       "risk_5y_ckd_epi_2012_cys", 
                       "risk_5y_ckd_epi_2021_cr", 
                       "risk_5y_ckd_epi_2009_cr")

# long format
cohort_risk_5y_long <- cohort |> 
  dplyr::select(lopnr,
                risk_5y_ckd_epi_2009_cr, 
                risk_5y_ckd_epi_2021_cr, 
                risk_5y_ckd_epi_2012_cys,
                risk_5y_ckd_epi_2012_cr_cys, 
                risk_5y_ckd_epi_2021_cr_cys) |> 
  tidyr::pivot_longer(cols = starts_with("risk_"), 
                      names_to = "Equation", 
                      values_to = "Risk") |> 
  dplyr::mutate(Equation = factor(Equation, levels = equation_order_5y),
                Equation = dplyr::recode(Equation,
                                         "risk_5y_ckd_epi_2009_cr" = "CKD-EPIcr 2009",
                                         "risk_5y_ckd_epi_2021_cr" = "CKD-EPIcr 2021",
                                         "risk_5y_ckd_epi_2012_cys" = "CKD-EPIcys 2012",
                                         "risk_5y_ckd_epi_2012_cr_cys" = "CKD-EPIcrcys 2012",
                                         "risk_5y_ckd_epi_2021_cr_cys" = "CKD-EPIcrcys 2021"))

# 0-100% density ridge plot
plot_KFRE_5y <- ggplot2::ggplot(cohort_risk_5y_long, aes(x = Risk, 
                                                         y = Equation, 
                                                         fill = stat(x))) +
  ggridges::geom_density_ridges_gradient(scale = 1, 
                                         rel_min_height = 0.000001, 
                                         alpha = 0.8) +  
  ggplot2::scale_fill_viridis_c(name = "Risk", option = "C", 
                                limits = c(0, 0.05), 
                                oob = scales::squish) +  
  ggplot2::scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +  
  ggplot2::scale_y_discrete(expand = c(0, 0)) +  # uniform spacing
  ggplot2::labs(title = "5-year KFRE risk 0-100%",
                x = "Risk Score", y = "Equation") +
  ggplot2::theme_minimal(base_size = 12) +  
  ggplot2::theme(panel.spacing = ggplot2::unit(0.5, "lines"), 
                 axis.text.y = ggplot2::element_text(size = 10))

# 0-1% Create density ridge plot
plot_KFRE_zoomed_01_5y <- ggplot2::ggplot(cohort_risk_5y_long, aes(x = Risk, 
                                                                   y = Equation, 
                                                                   fill = stat(x))) +
  ggridges::geom_density_ridges_gradient(scale = 0.9, 
                                         rel_min_height = 0.000001, 
                                         alpha = 0.8) +  
  ggplot2::scale_fill_viridis_c(name = "Risk", option = "C", 
                                limits = c(0, 0.008), 
                                oob = scales::squish) +  
  ggplot2::scale_x_continuous(limits = c(0, 0.01), 
                              breaks = seq(0, 0.01, by = 0.002)) +  
  ggplot2::scale_y_discrete(expand = c(0, 0)) +  # uniform spacing
  ggplot2::labs(title = "5y KFRE Risk Distributions by Equation 0-1%",
                x = "Risk Score", y = "Equation") +
  ggplot2::theme_minimal(base_size = 12) +  
  ggplot2::theme(
    plot.background = ggplot2::element_rect(fill = alpha("white", 0.8), color = NA), 
    panel.background = ggplot2::element_rect(fill = alpha("white", 0.8), color = NA),
    legend.title = ggplot2::element_text(size = 8),  
    legend.text = ggplot2::element_text(size = 6),   
    legend.key.size = ggplot2::unit(0.5, "lines"),   
    panel.spacing = ggplot2::unit(0.5, "lines"), 
    axis.text.y = ggplot2::element_text(size = 10),  
    axis.text.x = ggplot2::element_text(size = 8),   
    axis.title.x = ggplot2::element_text(size = 8),  
    axis.title.y = ggplot2::element_text(size = 8),  
    plot.title = ggplot2::element_text(size = 8))    

# Combine plots
plot_KFRE_combined_5y <- plot_KFRE_5y + 
  ggplot2::annotation_custom(ggplotGrob(plot_KFRE_zoomed_01_5y), 
                             xmin = 0.4, xmax = 1.05, # Horizontal placement of the zoomed plot
                             ymin = 3, ymax = 6) + 
  ggplot2::theme_classic() +
  ggplot2::theme(legend.position = "bottom",
                 legend.title = ggplot2::element_text(size = 8),  
                 legend.text = ggplot2::element_text(size = 6),   
                 plot.tag = ggplot2::element_blank(),  
                 axis.text.x = ggplot2::element_text(size = 10),
                 axis.text = ggplot2::element_text(size = 6),
                 axis.text.y = ggplot2::element_text(size = 10),
                 plot.background = ggplot2::element_rect(fill = "white", 
                                                         color = "white"),
                 panel.background = ggplot2::element_rect(fill = "white", 
                                                          color = "white"))

# combine plots for 2 and 5 years prediction horizon
combined_dist_plot <- plot_KFRE_combined_2y + plot_KFRE_combined_5y + 
  plot_layout(ncol = 1, guides = "collect")
combined_dist_plot & theme(legend.position = "bottom")
ggplot2::ggsave(filename = "Combined distribution plots.png", 
                plot = combined_dist_plot, width = 15, height = 10, dpi = 300)

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
  ggplot2::scale_fill_viridis_c(name = "Risk Difference", 
                                option = "C", 
                                limits = c(-0.0025, 0.0025), 
                                oob = scales::squish) +  
  ggplot2::scale_x_continuous(limits = c(-0.01, 0.01), breaks = seq(-0.01, 0.01, by = 0.001)) +  
  ggplot2::scale_y_discrete(expand = c(0, 0)) +  # uniform spacing 
  ggplot2::geom_vline(xintercept = 0, # reference line at zero difference
                      linetype = "dashed", 
                      color = "black") +  
  ggplot2::labs(title = "2y KFRE risk difference compared to CKD-EPIcr 2009",
                x = "Risk difference (Risk equation - Risk CKD-EPIcr 2009)", 
                y = "Equation") +
  ggplot2::theme_minimal(base_size = 12) +  
  ggplot2::theme(panel.spacing = ggplot2::unit(0.5, "lines"), 
                 axis.text.y = ggplot2::element_text(size = 10))

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
  ggplot2::scale_fill_viridis_c(name = "Risk Difference", 
                                option = "C", 
                                limits = c(-0.0025, 0.0025), 
                                oob = scales::squish) +  
  ggplot2::scale_x_continuous(limits = c(-0.01, 0.01),
                              breaks = seq(-0.01, 0.01, by = 0.001)) +  
  ggplot2::scale_y_discrete(expand = c(0, 0)) +  # uniform spacing
  ggplot2::geom_vline(xintercept = 0, # reference line at zero difference
                      linetype = "dashed", 
                      color = "black") +  
  ggplot2::labs(title = "5y KFRE Risk Differences Compared to CKD-EPIcr 2009",
                x = "Risk Difference (Risk Equation - Risk CKD-EPIcr 2009)", 
                y = "Equation") +
  ggplot2::theme_minimal(base_size = 12) +  
  ggplot2::theme(panel.spacing = ggplot2::unit(0.5, "lines"), 
                 axis.text.y = ggplot2::element_text(size = 10))

# combine plots
combined_plot_diff <- plot_KFRE_diff_2y + plot_KFRE_diff_5y + 
  plot_layout(ncol = 1, guides = "collect")
combined_plot_diff & theme(legend.position = "bottom")
ggplot2::ggsave(filename = "Combined risk difference.png", 
                plot = combined_plot_diff, width = 15, height = 10, dpi = 300)

################################################################################
### eGFR distributions #########################################################
################################################################################
# long format
cohort_egfr_long <- cohort |>
  tidyr::pivot_longer(
    cols = c(
      ckd_epi_2021_cr_cys, ckd_epi_2012_cr_cys,  # crcys equations
      ckd_epi_2012_cys,                          # cys equation
      ckd_epi_2021_cr, ckd_epi_2009_cr),         # cr equations
    names_to = "Equation",
    values_to = "eGFR"
  ) |> 
  dplyr::mutate(Equation = factor(Equation, levels = equations),
                Equation = dplyr::recode(Equation,
                                         "ckd_epi_2009_cr" = "CKD-EPIcr 2009",
                                         "ckd_epi_2021_cr" = "CKD-EPIcr 2021",
                                         "ckd_epi_2012_cys" = "CKD-EPIcys 2012",
                                         "ckd_epi_2012_cr_cys" = "CKD-EPIcrcys 2012",
                                         "ckd_epi_2021_cr_cys" = "CKD-EPIcrcys 2021"))

# create plot
egfr_plot <- ggplot2::ggplot(cohort_egfr_long, 
                             aes(x = eGFR, 
                                 y = Equation, 
                                 fill = stat(x))) +
  ggridges::geom_density_ridges_gradient(scale = 0.9, 
                                         rel_min_height = 0.001, 
                                         alpha = 0.8) +  
  ggplot2::scale_fill_viridis_c(name = "eGFR", option = "C", limits = c(0, 60), 
                                oob = scales::squish) +  
  ggplot2::scale_x_continuous(limits = c(0, 60), breaks = seq(0, 60, by = 10)) + 
  ggplot2::scale_y_discrete(expand = c(0, 0)) +  # uniform spacing 
  ggplot2::labs(title = "eGFR distributions",
                x = "eGFR (mL/min/1.73m²)", y = "Equation") +
  ggplot2::theme_minimal(base_size = 12) +  
  ggplot2::theme(panel.spacing = ggplot2::unit(1.5, "lines"), 
                 axis.text.y = ggplot2::element_text(size = 10)) +
  ggplot2::theme(plot.margin = margin(20, 10, 20, 10)) # Expands top/bottom margins
ggplot2::ggsave(filename = "eGFR distribution.png", 
                plot = egfr_plot, width = 15, height = 10, dpi = 300)

################################################################################
### Reclassification tables ####################################################
################################################################################
# 2y eGFR <30
cohort_risk_2y_reclas <- cohort |>
  dplyr::filter(ckd_epi_2009_cr >= 10 & ckd_epi_2009_cr <= 30)

# 5y eGFR 30-60
cohort_risk_5y_reclas <- cohort |> 
  dplyr::filter(ckd_epi_2009_cr >= 30 & ckd_epi_2009_cr <= 60)

for (equation in equations){
  # categorize 2-year predicted risk below and above 40%
  risk <- eval(parse(text=paste0("cohort_risk_2y_reclas$risk_2y_", equation)))
  cohort_risk_2y_reclas[, paste0("cat_risk_2y_", equation)] <- 
    ifelse(risk > 0.40, ">40%", "≤40%")
  
  # categorize 5-year predicted risk in three categories
  risk <- eval(parse(text=paste0("cohort_risk_5y_reclas$risk_5y_", equation)))
  cohort_risk_5y_reclas[, paste0("cat_risk_5y_", equation)] <- 
    ifelse(risk < 0.03,  "<3%",
           ifelse(risk >= 0.03 & risk <= 0.05, 
                  "3% - 5%", 
                  ">5%"))
}

# Split data based on outcome
cohort_risk_2y_reclas_cases <- subset(cohort_risk_2y_reclas, 
                                      outcome_2y == 1)     
cohort_risk_2y_reclas_control <- subset(cohort_risk_2y_reclas, 
                                        outcome_2y %in% c(0, 2)) 
cohort_risk_5y_reclas_cases <- subset(cohort_risk_5y_reclas, 
                                      outcome_5y == 1)     
cohort_risk_5y_reclas_control <- subset(cohort_risk_5y_reclas, 
                                        outcome_5y %in% c(0, 2)) 

# Make two reclassification tables for 2-year and 5-year predictions
table_2y <- c()
table_5y <- c()
for (equation in equations[-1]){
  # for all patients
  table_2y_all <- generate_reclass_table_2y(data = cohort_risk_2y_reclas, 
                                            risk_categories = eval(parse(text=paste0("cohort_risk_2y_reclas$cat_risk_2y_", equation))), 
                                            title = paste("Reclassification Table all patients after 2 years", equation))
  table_5y_all <- generate_reclass_table_5y(data = cohort_risk_5y_reclas, 
                                            risk_categories = eval(parse(text=paste0("cohort_risk_5y_reclas$cat_risk_5y_", equation))), 
                                            title = paste("Reclassification Table all patients after 5 years", equation))
  
  # for patients with the outcome
  table_2y_cases <- generate_reclass_table_2y(data = cohort_risk_2y_reclas_cases, 
                                              risk_categories = eval(parse(text=paste0("cohort_risk_2y_reclas_cases$cat_risk_2y_", equation))), 
                                              title = paste("Reclassification Table patients with KFRT after 2 years", equation))
  table_5y_cases <- generate_reclass_table_5y(data = cohort_risk_5y_reclas_cases, 
                                              risk_categories = eval(parse(text=paste0("cohort_risk_5y_reclas_cases$cat_risk_5y_", equation))), 
                                              title = paste("Reclassification Table patients with KFRT after 5 years", equation))
  
  # for patients without the outcome
  table_2y_control <- generate_reclass_table_2y(data = cohort_risk_2y_reclas_control, 
                                                risk_categories = eval(parse(text=paste0("cohort_risk_2y_reclas_control$cat_risk_2y_", equation))), 
                                                title = paste("Reclassification Table patients without KFRT after 2 years", equation))
  table_5y_control <- generate_reclass_table_5y(data = cohort_risk_5y_reclas_control, 
                                                risk_categories = eval(parse(text=paste0("cohort_risk_5y_reclas_control$cat_risk_5y_", equation))), 
                                                title = paste("Reclassification Table patients without KFRT after 5 years", equation))
  
  # final table
  table_2y <- cbind(table_2y,
                    rbind(table_2y_all$table, rep("", 2), 
                          table_2y_cases$table, rep("", 2), 
                          table_2y_control$table))
  table_5y <- cbind(table_5y,
                    rbind(table_5y_all$table, rep("", 2), 
                          table_5y_cases$table, rep("", 2), 
                          table_5y_control$table))
}
# write to Excel
openxlsx::write.xlsx(table_2y,
                     file = "Reclassification 2-year predictions.xlsx")
openxlsx::write.xlsx(table_5y,
                     file = "Reclassification 5-year predictions.xlsx")

