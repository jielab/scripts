#Author: Linsey Jackson
#Date: 12/9/2025
#Script name: linear regression models
#investigate heterogeneity of treatment effects (body weight change, waist-to-height ratio change) across baseline predicted risk strata (quartiles) 
#using linear regression models, whereby risk strata were coded numerically to test for linear trends

# ---- Setup ----
library(dplyr)
library(tidyr)
library(emmeans)
library(car)
library(purrr)
library(ggplot2)
library(stringr)
library(grid)   
library(forcats)
library(ggh4x)

# Inputs expected:
# df_joined         : data.frame with endpoints, arms, and strata columns
# strata_columns    : character vector of column names used to define Stratum
# endpoint_columns  : character vector of continuous endpoint column names

# Reset to treatment contrasts
options(contrasts = c("contr.treatment", "contr.poly"))

# choose endpoints
endpoint_columns <- c("PCHG_Weight (kg)","CHG_Weight (kg)", "WHtR_Change") 

# Storage
out_emm    <- list()
out_interaction  <- list()
out_simple <- list()

for (strata_col in strata_columns) {
  
  df_strata <- df_joined %>%
    filter(!is.na(.data[[strata_col]])) %>%
    mutate(
      Stratum = case_when(
        .data[[strata_col]] == "Q1" ~ 1,
        .data[[strata_col]] == "Q2" ~ 2,
        .data[[strata_col]] == "Q3" ~ 3,
        .data[[strata_col]] == "Q4" ~ 4,
        TRUE ~ NA_real_
      ),
      # Set factor with Placebo first (as reference)
      ARM = factor(.data[["ARM"]], 
                   levels = c("Placebo", "TZP 5mg", "TZP 10mg", "TZP 15mg"))
    ) %>%
    droplevels()
  
  contrasts(df_strata$ARM) <- contr.treatment(levels(df_strata$ARM))
  
  for (endpoint_col in endpoint_columns) {
    if (!endpoint_col %in% names(df_strata)) next
    
    # Fit ONE model with all arms
    form <- as.formula(sprintf("`%s` ~ Stratum * ARM", endpoint_col))
    fit  <- lm(form, data = df_strata)
    
    # Extract interaction coefficients (these test linear trend differences vs Placebo)
    coef_summary <- summary(fit)$coefficients
    interaction_rows <- grep("^Stratum:ARM", rownames(coef_summary))
    
    if (length(interaction_rows) > 0) {
      interaction_tests <- tibble(
        Treatment_ARM = gsub("Stratum:ARM", "", rownames(coef_summary)[interaction_rows]),
        Interaction_P_Value = coef_summary[interaction_rows, "Pr(>|t|)"],
        Interaction_Estimate = coef_summary[interaction_rows, "Estimate"],
        Interaction_SE = coef_summary[interaction_rows, "Std. Error"],
        Risk_Score_Strata = strata_col,
        Endpoint = endpoint_col
      )
    } else {
      interaction_tests <- tibble()
    }
    
    # EMMs for plotting
    emm_obj <- emmeans(fit, ~ ARM | Stratum, at = list(Stratum = 1:4))
    emm_cells <- as.data.frame(emm_obj) %>%
      mutate(Risk_Score_Strata = strata_col, Endpoint = endpoint_col)
    
    # Simple effects
    simple_df <- tryCatch({
      pairs(emm_obj) %>% 
        as.data.frame() %>%
        filter(grepl("Placebo", contrast))
    }, error = function(e) NULL)
    
    if (!is.null(simple_df) && nrow(simple_df) > 0) {
      simple_df <- simple_df %>%
        mutate(Risk_Score_Strata = strata_col, Endpoint = endpoint_col)
    } else {
      simple_df <- tibble()
    }
    
    key <- paste(strata_col, endpoint_col, sep = "_")
    out_emm[[key]]    <- emm_cells
    out_interaction[[key]]  <- interaction_tests
    if (nrow(simple_df) > 0) out_simple[[key]] <- simple_df
  }
}

emm_all    <- bind_rows(out_emm)
interaction_all  <- bind_rows(out_interaction)
simple_all <- if (length(out_simple)) bind_rows(out_simple) else tibble()

#plot

# ---- Helpers for plotting ----
safe_name <- function(x) {
  x %>%
    str_replace_all("[^A-Za-z0-9_-]+", "_") %>%
    str_replace_all("_+", "_") %>%
    str_replace("^_|_$", "")
}

library(stringr)

# Create table with interaction p-values (now one per treatment arm)
int_tbl <- interaction_all %>%
  mutate(
    Formatted_Strata = Risk_Score_Strata %>%
      str_replace("_group$", "") %>%              # remove suffix
      str_replace_all("[._]", " ") %>%            # replace both underscores and periods with space
      str_squish()
  ) %>%
  group_by(Risk_Score_Strata, Endpoint, Formatted_Strata) %>%
  summarise(
    Facet_Label = paste0(
      Formatted_Strata, "\n",
      "Interaction p-values (vs Placebo):\n",
      paste(Treatment_ARM, ": p=", signif(Interaction_P_Value, 3), collapse = "\n", sep = "")
    ),
    .groups = "drop"
  )


# -------- plot controls --------
facet_cols        <- 6
base_font_size    <- 9
axis_text_size    <- 7
strip_text_size   <- 6
legend_text_size  <- 8
legend_title_size <- 8
point_size        <- 2.0
line_width        <- 0.6
dodge_w           <- 0.20
err_bar_w         <- 0.15
panel_space_lines <- 0.4
row_height_in     <- 2.0
min_pdf_height    <- 6
pdf_width         <- 12

# ---- Plot: one PDF per Endpoint, facets = Formatted_Strata (single page) ----
for (ep in unique(emm_all$Endpoint)) {
  
  df_ep <- emm_all %>%
    filter(Endpoint == ep) %>%
    mutate(
      Stratum = as.numeric(Stratum),
      ARM     = factor(ARM)
    ) %>%
    left_join(
      int_tbl %>% select(Risk_Score_Strata, Endpoint, Facet_Label, Formatted_Strata),
      by = c("Risk_Score_Strata", "Endpoint")
    )
  
  if (nrow(df_ep) == 0) next
  
  # --- Common y-range across all facets for this endpoint ---
  y_min <- min(df_ep$lower.CL, na.rm = TRUE)
  y_max <- max(df_ep$upper.CL, na.rm = TRUE)
  
  if (!is.finite(y_min) || !is.finite(y_max)) next
  if (y_min == y_max) {
    y_min <- y_min - 0.5
    y_max <- y_max + 0.5
  }
  
  rng   <- y_max - y_min
  pad   <- 0.04 * rng
  ylims <- range(c(0, y_min - pad, y_max + pad), na.rm = TRUE)
  custom_breaks <- if (ep == "PCHG_Weight (kg)") c(-25,-20, -15, -10, -5, 0) else waiver()
  
  # ---- Build labeller ----
  labs_vec <- df_ep %>%
    distinct(Formatted_Strata, Facet_Label) %>%
    filter(!is.na(Facet_Label)) %>%
    { setNames(.$Facet_Label, .$Formatted_Strata) }
  
  n_facets   <- dplyr::n_distinct(df_ep$Formatted_Strata)
  facet_rows <- max(1, ceiling(n_facets / facet_cols))
  pdf_height <- max(min_pdf_height, row_height_in * facet_rows)
  
  pdf(file = paste0("Figures_", safe_name(ep), "_AllRiskScores_Faceted.pdf"),
      width = pdf_width, height = pdf_height)
  
  p <- ggplot(df_ep, aes(x = Stratum, y = emmean, group = ARM, color = ARM)) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                  width = err_bar_w, position = position_dodge(width = dodge_w)) +
    geom_line(linewidth = line_width, position = position_dodge(width = dodge_w)) +
    geom_point(size = point_size, position = position_dodge(width = dodge_w)) +
    scale_color_manual(
      values = c( "Placebo" = "darkgray",  
                  "TZP 5mg"     = "#FBCFC8",  
                  "TZP 10mg"    = "#F58E7D", 
                  "TZP 15mg"    = "#E1251B")
    ) +
    ggh4x::facet_wrap2(
      ~ Formatted_Strata,
      ncol    = facet_cols,
      scales  = "fixed",
      labeller = as_labeller(labs_vec),
      axes    = "all"
    ) +
    scale_y_continuous(
      limits = ylims,
      expand = expansion(mult = 0),
      breaks = custom_breaks
    ) +
    labs(
      title    = ep,
      subtitle = "Each panel shows a risk-score group; lines/points colored by Treatment Arm",
      x = "Risk Score Strata",
      y = "Adjusted Mean (EMM)",
      color = "Treatment Arm"
    ) +
    theme_minimal(base_size = base_font_size) +
    theme(
      plot.title    = element_text(face = "bold"),
      axis.text.x   = element_text(angle = 45, hjust = 1, size = axis_text_size),
      axis.text.y   = element_text(size = axis_text_size),
      strip.text    = element_text(size = strip_text_size, face = "bold"),
      panel.spacing = unit(panel_space_lines, "lines"),
      legend.position = "bottom",
      legend.text   = element_text(size = legend_text_size),
      legend.title  = element_text(size = legend_title_size),
      legend.key.size = unit(0.5, "lines"),
      plot.margin   = margin(10, 10, 20, 10)
    )
  
  print(p)
  dev.off()
}