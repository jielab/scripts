#############################################################
####                 Risk stratification                 ####
#############################################################

# This script generates:
## Fig 5A  : Absolute observed 10-y event rates across deciles of estimated risk
## Fig 5B  : Estimated risk distribution across BMI categories
## Fig 5C  : Cumulative incidences across quintiles of CV death risk
## SFig 10 : FDR-DR plots
## SFig 14 : Rate ratios 50% top / 50% bottom risk vs. obesity / overweight
## Data for: Fig 4C,D,E
## Data for: Supp Fig 8,9

rm(list = ls())
options(stringsAsFactors = F)
options(scipen = 1)
setwd("<path>")

## --> packages needed <-- ##
require(data.table)
require(tidyverse)
require(readxl)
require(dplyr)
require(magrittr)
require(ggsci)
require(patchwork)
require(ggridges)
require(survminer)
require(survival)

## import data
ukb.phe <-  fread("<path>")
ukb.phe.original <- ukb.phe

## import labels
lab.phe <- fread("<path>")

## define outcomes
outc_list <- c("<outcomes>")

## define follow-up
cdat_list <- c("<followup>")

## rename outcomes
new_labels <- c("<nicenames>")

## reorder outcomes
custom_order <- c("<nicenames>")

########################################################
####            Absolute risks  (Fig 5A)            ####
########################################################

## data frames to store
merged_results <- data.frame()
results_bmi_cat <- data.frame()
results_lin.pred <- data.frame()

# loop across outcomes
for (i in seq_along(outc_list)) {
  outc <- outc_list[i]
  cdat <- cdat_list[i]
  
  set.seed(123)
  
  ## remove participants with event during first 6 months
  ukb.phe <-
    ukb.phe.original %>% filter(!(((
      eval(as.name(cdat)) < 0.5
    )) & (eval(as.name(
      outc
    )) == 1)))
  
  ukb.phe <- as.data.frame(ukb.phe)
  ##--  create data splits  --##
  # validation
  v.idx <- ukb.phe$f.eid[ukb.phe$split == "v.idx"]
  
  ## load models
  loaded_object_name <-
    load(file = paste0("<path>", outc, ".RData"))
  feature.com <- get(loaded_object_name)
  
  ukb.test <- ukb.phe %>% filter(f.eid %in% v.idx)
  ukb.test <- as.data.table(ukb.test)
  
  ## join linear predictor to dataset
  ukb.test[, lin.pred := feature.com$linear.predictor]
  
  ## absolute risks across quintiles
  res <- ukb.test %>%
    mutate(decile = ntile(lin.pred, 5)) %>%
    mutate(bmi_category = case_when(
      bmi >= 27 & bmi < 30 ~ "Overweight (BMI 27-30)",
      bmi >= 30 ~ "Obesity (BMI >30)"
    )) %>%
    group_by(decile) %>%
    summarise(count = sum(eval(as.name(outc)) == 1, na.rm = TRUE),
              total_participants = n()) %>%
    mutate(
      percentage = (count / total_participants) * 100,
      decile = as.factor(decile),
      outcome = outc
    )
  
  ## absolute risks across 50% top vs bottom
  res_half <- ukb.test %>%
    mutate(decile_50 = ntile(lin.pred, 2)) %>%
    group_by(decile_50) %>%
    summarise(count = sum(eval(as.name(outc)) == 1, na.rm = TRUE),
              total_participants = n()) %>%
    mutate(
      percentage = (count / total_participants) * 100,
      decile_50 = as.factor(decile_50),
      outcome = outc
    )
  
  ## absolute risks obesity and overweight
  res_BMI_based <- ukb.test %>%
    mutate(bmi_category = case_when(
      bmi >= 27 & bmi < 30 ~ "Overweight (BMI 27-30)",
      bmi >= 30 ~ "Obesity (BMI >30)"
    )) %>%
    group_by(bmi_category) %>%
    summarise(count = sum(eval(as.name(outc)) == 1, na.rm = TRUE),
              total_participants = n()) %>%
    mutate(
      percentage = (count / total_participants) * 100,
      bmi_category = as.factor(bmi_category),
      outcome = outc
    )
  
  ## absolute risks across quintiles stratified by overweight/obesity
  res_overweight_obesity <- ukb.test %>%
    mutate(decile = ntile(lin.pred, 5)) %>%
    mutate(category = case_when(
      bmi >= 27 & bmi < 30 ~ "Overweight (BMI 27-30)",
      bmi >= 30 ~ "Obesity (BMI >30)"
    )) %>%
    group_by(decile, category) %>%
    summarise(count = sum(eval(as.name(outc)) == 1, na.rm = TRUE),
              total_participants = n()) %>%
    mutate(
      percentage = (count / total_participants) * 100,
      decile = as.factor(decile),
      outcome = outc
    )
  
  ## rate ratios top vs bottom 20%
  ratios <- res %>%
    group_by(outcome) %>%
    summarise(ratio_decile5_to_decile1 = percentage[decile == 5] / percentage[decile == 1])
  
  ## rate ratios top vs bottom 50%
  ratios_half <- res_half %>%
    group_by(outcome) %>%
    summarise(ratio_50_50 = percentage[decile_50 == 2] / percentage[decile_50 == 1])
  
  ## rate ratios obesity vs overweight
  ratios_BMI_based <- res_BMI_based %>%
    group_by(outcome) %>%
    summarise(ratio_ob_ow = percentage[bmi_category == "Obesity (BMI >30)"] / percentage[bmi_category == "Overweight (BMI 27-30)"])
  
  ## estimated risks
  res_lin.pred <- data.frame(outcome = outc,
                             lin.pred = ukb.test$lin.pred,
                             bmi = ukb.test$bmi)
  
  ## merge the ratio back into the original dataset
  final_result <- res %>%
    left_join(ratios, by = c("outcome" = "outcome")) %>%
    left_join(ratios_half, by = "outcome") %>%
    left_join(ratios_BMI_based, by = "outcome")
  
  results_lin.pred <- rbind(results_lin.pred, res_lin.pred)
  results_bmi_cat <- rbind(results_bmi_cat, res_overweight_obesity)
  merged_results <- rbind(merged_results, final_result)
}

## rename outcomes
merged_results$outcome <- new_labels[merged_results$outcome]
## reorder outcomes
merged_results$outcome <-
  factor(merged_results$outcome, levels = unique(custom_order))

## rename outcomes
results_bmi_cat$outcome <- new_labels[results_bmi_cat$outcome]
## reorder outcomes
results_bmi_cat$outcome <-
  factor(results_bmi_cat$outcome, levels = unique(custom_order))

## rename outcomes
results_lin.pred$outcome <- new_labels[results_lin.pred$outcome]
## reorder outcomes
results_lin.pred$outcome <-
  factor(results_lin.pred$outcome, levels = unique(custom_order))

results_bmi_cat$category <-
  factor(results_bmi_cat$category,
         levels = c("Overweight (BMI 27-30)", "Obesity (BMI >30)"))

## plot absolute risks
a <- merged_results %>%
  ggplot(aes(x = decile, y = percentage)) +
  facet_wrap(~ outcome, scales = "free_y", ncol = 3) +
  geom_bar(
    data = merged_results %>% filter(decile != 5),
    aes(x = decile, y = percentage),
    stat = "identity",
    position = position_dodge2(),
    fill = "#80796BFF",
    color = "black",
    width = .8,
    alpha = .9
  ) +
  geom_bar(
    data = results_bmi_cat %>% filter(decile == 5),
    aes(x = decile, y = percentage, fill = category),
    stat = "identity",
    position = position_dodge2(),
    color = "black",
    width = .9,
    linewidth = .25,
    alpha = .9
  ) +
  scale_fill_manual(
    values = c("#DF8F44FF", "#374E55FF"),
    guide = guide_legend(
      ncol = 1,
      byrow = TRUE,
      keyheight = unit(0.2, "cm"),
      keywidth = unit(0.4, "cm")
    )
  ) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(x = "Estimated risk quintiles", y = "Actual 10-year event risks") +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    legend.position = "top",
    legend.key.size = unit(0.4, 'cm'),
    legend.spacing = unit(0, 'cm'),
    strip.text = element_text(color = "white"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.background = element_rect(
      color = "black",
      fill = "black",
      linetype = "solid"
    )
  )


#############################################
####    5050 RR vs. OB-OW RR (SFig 14)   ####
#############################################
RR_quintiles <- merged_results %>%
  group_by(outcome) %>%
  mutate(
    rate_decile5 = percentage[decile == "5"],
    rr_5_vs_4 = rate_decile5 / percentage[decile == "4"],
    rr_5_vs_3 = rate_decile5 / percentage[decile == "3"],
    rr_5_vs_2 = rate_decile5 / percentage[decile == "2"],
    rr_5_vs_1 = rate_decile5 / percentage[decile == "1"]
  ) %>%
  ungroup() %>%
  filter(decile == 5) %>%
  dplyr::select(outcome, rr_5_vs_1, rr_5_vs_2, rr_5_vs_3, rr_5_vs_4)

RR_50 <-
  merged_results %>% filter(decile == 5) %>% dplyr::select(outcome, ratio_50_50)

rate_ratios <- RR_50 %>% left_join(RR_quintiles, by = "outcome") %>%
  arrange(desc(ratio_50_50))

fwrite(rate_ratios, "<path>/rateratios.csv")

df_ratios <- merged_results %>% filter(decile == 5) %>%
  select(outcome, ratio_50_50, ratio_ob_ow)

df_ratios <- df_ratios %>%
  mutate(outcome = reorder(outcome, -ratio_50_50))

df_long <- df_ratios %>%
  pivot_longer(c(ratio_50_50, ratio_ob_ow),
               names_to  = "ratio_type",
               values_to = "ratio") %>%
  mutate(ratio_type = recode(ratio_type,
                             ratio_50_50 = "50/50",
                             ratio_ob_ow = "OB/OW"))


pdf("<path>/5050_OBOW_RateRatios.pdf",
    width = 7,
    height = 6)
ggplot() +
  geom_hline(yintercept = 1,
             linetype = "dashed",
             color = "grey60") +
  geom_segment(
    data = df_ratios,
    aes(
      x = outcome,
      xend = outcome,
      y = ratio_ob_ow + 0.3,
      yend = ratio_50_50 - 0.3
    ),
    arrow = arrow(length = unit(1.5, "mm")),
    linewidth = 0.6,
    color = "grey40"
  ) +
  geom_point(
    data = df_long,
    aes(
      x = outcome,
      y = ratio,
      fill = ratio_type,
      shape = ratio_type
    ),
    size = 3
  ) +
  geom_text(
    data = df_long,
    aes(
      x = outcome,
      y = ratio,
      label = sprintf("%.2f", ratio),
      vjust = ifelse(ratio_type == "OB/OW", 1.8, -1),
      color = ratio_type
    ),
    size = 3.2,
    show.legend = FALSE
  ) +
  scale_shape_manual(values = c(`OB/OW` = 21, `50/50` = 23)) +
  scale_color_manual(values = c(`OB/OW` = "grey20", `50/50` = "black")) +
  scale_fill_manual(values = c("#0072B5FF", "#FFDC91FF")) +
  labs(x = "Outcome",
       y = "Rate ratio",
       shape = "Ratio",
       fill  = "Ratio") +
  ylim(0, 15) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_blank()
  )
dev.off()

########################################################
####            Risk Distribution (Fig 5B)          ####
########################################################
## reorder outcomes
results_lin.pred$outcome <-
  factor(results_lin.pred$outcome, levels = unique(rev(custom_order)))

## create detailed BMI categories
results_lin.pred <- results_lin.pred %>%
  mutate(
    bmi_category = case_when(
      bmi >= 27 & bmi < 30 ~ "Overweight (BMI 27-30)",
      bmi >= 30 & bmi < 35 ~ "Class 1 Obesity (BMI 30-30)",
      bmi >= 35 & bmi < 40 ~ "Class 2 Obesity (BMI 35-35)",
      bmi >= 40 ~ "Class 3 Obesity (BMI >40)",
      TRUE ~ NA_character_
    )
  )

## colors
jama_colors <- pal_jama("default")(7)
reordered_colors <- jama_colors[c(3, 4, 1, 2)]

b <- results_lin.pred %>%
  ggplot(aes(x = lin.pred, y = outcome, fill = bmi_category)) +
  stat_density_ridges(quantile_lines = TRUE,
                      alpha = .7,
                      color = NA) +
  labs(x = "Estimated risk", y = "") +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_blank(),
    legend.position = "top",
    legend.key.size = unit(0.25, 'cm'),
    legend.spacing = unit(0, 'cm'),
    panel.spacing = unit(0.1, "lines"),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(
      color = "black",
      fill = "black",
      linetype = "solid"
    )
  ) +
  scale_fill_manual(
    values = reordered_colors,
    guide = guide_legend(
      ncol = 1,
      byrow = TRUE,
      keyheight = unit(0.2, "cm"),
      keywidth = unit(0.4, "cm")
    )
  )

## plot with overweight vs BMI >45 (Revision Comment)
pdf("<path>/reviewerfig_BMI45.pdf",
    height = 6,
    width = 7)
results_lin.pred %>%
  mutate(
    bmi_category = case_when(
      bmi >= 27 & bmi < 30 ~ "(BMI 27-30)",
      bmi >= 45 ~ "(BMI >45)",
      TRUE ~ NA_character_
    )
  ) %>% filter(!is.na(bmi_category)) %>%
  ggplot(aes(x = lin.pred, y = outcome, fill = bmi_category)) +
  geom_boxplot(outlier.size = 1) +
  labs(x = "Estimated risk", y = "") +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_blank(),
    legend.position = "right",
    legend.key.size = unit(0.25, 'cm'),
    legend.spacing = unit(0, 'cm'),
    panel.spacing = unit(0.1, "lines"),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(
      color = "black",
      fill = "black",
      linetype = "solid"
    )
  ) +
  scale_fill_manual(
    values = c("darkred", "#DF8F44FF"),
    guide = guide_legend(
      ncol = 1,
      byrow = TRUE,
      keyheight = unit(0.2, "cm"),
      keywidth = unit(0.4, "cm")
    )
  )
dev.off()

#################################################################
####    Pie charts proportions of T2D/overweight (Fig 5A)    ####
#################################################################

## proportion of obesity/overweight individuals in decile 5
bmi_pie_data <- results_bmi_cat %>%
  filter(decile == 5) %>%
  group_by(outcome) %>%
  summarise(
    total_obesity = sum(total_participants[category == "Obesity (BMI >30)"]),
    total_overweight = sum(total_participants[category == "Overweight (BMI 27-30)"])
  ) %>%
  mutate(
    obesity_prop = total_obesity / (total_obesity + total_overweight),
    overweight_prop = total_overweight / (total_obesity + total_overweight)
  ) %>%
  pivot_longer(
    cols = c(obesity_prop, overweight_prop),
    names_to = "category",
    values_to = "percentage"
  ) %>%
  mutate(category = recode(
    category,
    "obesity_prop" = "Obesity",
    "overweight_prop" = "Overweight"
  ))

a <- bmi_pie_data %>%
  ggplot(aes(x = "", y = percentage, fill = category)) +
  geom_bar(
    stat = "identity",
    width = 1,
    color = "black",
    alpha = 0.9
  ) +
  scale_fill_manual(values = c("#374E55FF", "#DF8F44FF")) +
  coord_polar("y", start = 0) +
  facet_wrap(~ outcome, ncol = 3) +
  labs(title = "Proportion of overweight and obesity in decile 5",
       fill = "BMI Category") +
  theme_void() +
  theme(legend.position = "none")

## plot and save
pdf("graphics/piechartsbmiandt2d.pdf",
    width = 7,
    height = 7)
a
dev.off()

### ---- Metrics for the text ---- ###
## median (range) proportion of overweight
round(median(as.data.table(bmi_pie_data)[category == "Overweight"]$percentage), 4)
round(range(as.data.table(bmi_pie_data)[category == "Overweight"]$percentage), 3)
# 0.399 (0.198 - 0.455)

########################################################
####     Cumulative incidence CV Death (Fig 5C)     ####
########################################################
## set seed
set.seed(123)

## remove participants with event during first 6 months
ukb.phe <-
  ukb.phe.original %>% filter(!(fol.death < 0.5 & CVdeath == 1))

ukb.phe <- as.data.frame(ukb.phe)
##--  create data splits  --##
# validation
v.idx <- ukb.phe$f.eid[ukb.phe$split == "v.idx"]
ukb.phe  <- as.data.frame(ukb.phe)

## save results from model
loaded_object_name <-
  load(file = paste0("<path>"))
feature.com <- get(loaded_object_name)

ukb.test <- ukb.phe %>% filter(f.eid %in% v.idx)
ukb.test <- as.data.table(ukb.test)
ukb.test[, lin.pred := feature.com$linear.predictor]

## quintiles
quantiles <-
  quantile(ukb.test$lin.pred, probs = c(0.2, 0.4, 0.6, 0.8))
ukb.test$group <- with(ukb.test,
                       ifelse(
                         lin.pred <= quantiles[1],
                         "Bottom 20%",
                         ifelse(
                           lin.pred <= quantiles[2],
                           "20-40%",
                           ifelse(
                             lin.pred <= quantiles[3],
                             "40-60%",
                             ifelse(lin.pred <= quantiles[4],
                                    "60-80%",
                                    "Top 20%")
                           )
                         )
                       ))

## remove rest of participants
ukb.filtered <- ukb.test %>% filter(!is.na(group))

## fit survival model for cumulative incidence
surv_fit <-
  survfit(Surv(fol.death, CVdeath) ~ group, data = ukb.filtered)

## tidy survival data for ggplot
surv_data <- surv_summary(surv_fit) %>%
  mutate(group = strata) %>%
  separate(group, into = c("Group", "Label"), sep = "=")

## plot cumulative incidence
c <-
  ggplot(surv_data, aes(
    x = time,
    y = 1 - surv,
    color = Label,
    fill = Label
  )) +
  geom_step(size = .5) +
  geom_ribbon(
    data = surv_data,
    aes(
      ymin = 1 - lower,
      ymax = 1 - upper,
      fill = Label
    ),
    alpha = 0.3,
    color = NA
  ) +
  labs(x = "Follow-up (years)",
       y = "Cumulative incidence for\ncardiovascular death (%)") +
  theme_light() +
  scale_fill_manual(
    values = c(
      "#374E55FF",
      "#3C5488FF",
      "#4DBBD5FF",
      "#00A087FF",
      "#E64B35FF"
    ),
    guide = guide_legend(
      ncol = 2,
      byrow = TRUE,
      keyheight = unit(0.2, "cm"),
      keywidth = unit(0.4, "cm")
    )
  ) +
  scale_color_manual(values = c(
    "#374E55FF",
    "#3C5488FF",
    "#4DBBD5FF",
    "#00A087FF",
    "#E64B35FF"
  )) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_x_continuous(breaks = seq(0, 10, by = 2),
                     expand = expansion(mult = c(0, 0.04))) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )


layout <- "
AAAAAAABBBB
AAAAAAABBBB
AAAAAAABBBB
AAAAAAABBBB
AAAAAAABBBB
AAAAAAACCCC
AAAAAAACCCC
"

## plot and save
pdf("<path>/Fig_5.pdf", height = 7.25, width = 7.5)
a + b + c +
  plot_layout(design = layout) &
  theme(text = element_text(size = 8))
dev.off()

########################################################
####                DRs at FDR 5-95                 ####
########################################################

## data frame to store results
all_outcomes_results <- data.frame()

## compare three models
model_names <- c("LOFO.mean", "top.+PGS", "ASCVD")

## do for each outcome
for (i in seq_along(outc_list)) {
  outc <- outc_list[i]
  cdat <- cdat_list[i]
  
  set.seed(123)
  ## filter out participants with event during the first 6 months
  ukb.phe <-
    ukb.phe.original %>% filter(!(((
      eval(as.name(cdat)) < 0.5
    )) & (eval(as.name(
      outc
    )) == 1)))
  
  ukb.phe <- as.data.frame(ukb.phe)
  ##--  create data splits  --##
  # validation
  v.idx <- ukb.phe$f.eid[ukb.phe$split == "v.idx"]
  # convert to data frame
  ukb.phe  <- as.data.frame(ukb.phe)
  
  ## do for each model
  for (model_name in model_names) {
    ## load the model
    model_path <-
      paste0("<path>",
             model_name,
             ".",
             outc,
             ".RData")
    loaded_object_name <- load(file = model_path)
    feature.com <- get(loaded_object_name)
    
    ## get and match linear predictor of test set
    ukb.test <- ukb.phe %>% filter(f.eid %in% v.idx)
    ukb.test <- as.data.table(ukb.test)
    ukb.test[, lin.pred := feature.com$linear.predictor]
    
    ## normalize the linear predictor
    lp.test <- ukb.test$lin.pred
    if (min(lp.test) < 0) {
      lp.test <-
        (lp.test + abs(min(lp.test))) / max(lp.test + abs(min(lp.test)))
    } else {
      lp.test <- lp.test / max(lp.test)
    }
    
    ## calculate FPR and DR for each threshold
    thresholds <- seq(0.0001, 1, by = 0.0001)
    thresh.grid <- lapply(thresholds, function(threshold) {
      pred.case <- ifelse(lp.test > threshold, 1, 0)
      
      tp <- sum(ukb.test[[outc]] == 1 & pred.case == 1)
      fp <- sum(ukb.test[[outc]] == 0 & pred.case == 1)
      fn <- sum(ukb.test[[outc]] == 1 & pred.case == 0)
      tn <- sum(ukb.test[[outc]] == 0 & pred.case == 0)
      
      fpr <- if ((tn + fp) > 0)
        fp / (tn + fp)
      else
        NA
      dr <- if ((tp + fn) > 0)
        tp / (tp + fn)
      else
        NA
      
      data.frame(threshold = threshold,
                 fpr = fpr * 100,
                 dr = dr * 100)
    })
    thresh.grid <- do.call(rbind, thresh.grid)
    
    ## select closest thresholds to target FPR levels
    fpr_targets <-
      c(5,
        10,
        15,
        20,
        25,
        30,
        35,
        40,
        45,
        50,
        55,
        60,
        65,
        70,
        75,
        80,
        85,
        90,
        95)
    
    thresh.selected <- lapply(fpr_targets, function(target) {
      thresh.grid[which.min(abs(thresh.grid$fpr - target)), ]
    })
    thresh.selected <- do.call(rbind, thresh.selected)
    thresh.selected$FPR.true <- fpr_targets
    thresh.selected$outcome <- outc
    thresh.selected$model <- model_name
    
    ## add to cumulative results
    all_outcomes_results <-
      rbind(all_outcomes_results, thresh.selected)
  }
}

## replace outcome names and set custom order
all_outcomes_results$outcome <-
  new_labels[all_outcomes_results$outcome]
all_outcomes_results$outcome <-
  factor(all_outcomes_results$outcome, levels = custom_order)

all_outcomes_results$model <- factor(all_outcomes_results$model,
                                     levels = c("top.+PGS", "LOFO.mean", "ASCVD"))

all_outcomes_results[all_outcomes_results$model == "top.+PGS"]$model <-
  "Extended outcome-specific model"
all_outcomes_results[all_outcomes_results$model == "LOFO.mean"]$model <-
  "Unified core clinical model"
all_outcomes_results[all_outcomes_results$model == "ASCVD"]$model <-
  "ASCVD score model"

## reorder outcomes
all_outcomes_results$outcome <-
  factor(all_outcomes_results$outcome, levels = unique(custom_order))

## save results
fwrite(
  all_outcomes_results,
  "<path>/FDRresults.txt",
  row.names = F,
  na = NA
)

## plot DR/FDR
pdf("<path>/FDRPlots.pdf", width = 7, height = 8)
ggplot(data = all_outcomes_results, aes(
  x = FPR.true,
  y = dr,
  fill = model,
  color = model
)) +
  facet_wrap(~ outcome, ncol = 5) +
  geom_line() +
  geom_point(size = .5,
             shape = 21,
             color = "black") +
  labs(x = "False positive rate (FPR %)", y = "Detection rate (DR %)") +
  ylim(0, 100) +
  scale_color_jama() +
  scale_fill_jama() +
  theme_light() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 10),
    strip.text = element_text(color = "white"),
    strip.background = element_rect(
      color = "black",
      fill = "black",
      linetype = "solid"
    )
  )
dev.off()

########################################################
####       LR, pre-/post-test prob, at FDR 10       ####
########################################################

## data frame to store results
all_outcomes_results <- data.frame()

model_names <-
  c("LOFO.mean", "ASCVD", "SCORE2", "MHO_Zembic", "AgeSexBMI")

## across each outcome
for (i in seq_along(outc_list)) {
  outc <- outc_list[i]
  cdat <- cdat_list[i]
  
  set.seed(123)
  ## filter participants with <6 months follow-up
  ukb.phe <-
    ukb.phe.original %>% filter(!(((
      eval(as.name(cdat)) < 0.5
    )) & (eval(as.name(
      outc
    )) == 1)))
  
  ##--  create data splits  --##
  # validation
  v.idx <- ukb.phe$f.eid[ukb.phe$split == "v.idx"]
  # convert to data frame
  ukb.phe  <- as.data.frame(ukb.phe)
  
  # do for each model
  for (model_name in model_names) {
    # load the model
    model_path <-
      paste0("<path>",
             model_name,
             ".",
             outc,
             ".RData")
    loaded_object_name <- load(file = model_path)
    feature.com <- get(loaded_object_name)
    
    # get linear predictor from each model on test set
    ukb.test <- ukb.phe %>% filter(f.eid %in% v.idx)
    ukb.test <- as.data.table(ukb.test)
    ukb.test[, lin.pred := feature.com$linear.predictor]
    
    # compute pre-test probability
    pre_test_prob <- sum(ukb.test[[outc]] == 1) / nrow(ukb.test)
    
    # normalize the linear predictor
    lp.test <- ukb.test$lin.pred
    if (min(lp.test) < 0) {
      lp.test <-
        (lp.test + abs(min(lp.test))) / max(lp.test + abs(min(lp.test)))
    } else {
      lp.test <- lp.test / max(lp.test)
    }
    
    # calculate FPR and DR for each threshold
    thresholds <- seq(0.0001, 1, by = 0.0001)
    thresh.grid <- lapply(thresholds, function(threshold) {
      pred.case <- ifelse(lp.test > threshold, 1, 0)
      
      tp <- sum(ukb.test[[outc]] == 1 & pred.case == 1)
      fp <- sum(ukb.test[[outc]] == 0 & pred.case == 1)
      fn <- sum(ukb.test[[outc]] == 1 & pred.case == 0)
      tn <- sum(ukb.test[[outc]] == 0 & pred.case == 0)
      
      fpr <- if ((tn + fp) > 0)
        fp / (tn + fp)
      else
        NA
      dr <- if ((tp + fn) > 0)
        tp / (tp + fn)
      else
        NA
      
      data.frame(threshold = threshold,
                 fpr = fpr * 100,
                 dr = dr * 100)
    })
    thresh.grid <- do.call(rbind, thresh.grid)
    
    # select closest thresholds to target FPR levels
    fpr_targets <- c(10)
    thresh.selected <- lapply(fpr_targets, function(target) {
      thresh.grid[which.min(abs(thresh.grid$fpr - target)), ]
    })
    thresh.selected <- do.call(rbind, thresh.selected)
    thresh.selected$FPR.true <- fpr_targets
    thresh.selected$outcome <- outc
    thresh.selected$model <- model_name
    
    # calculate likelihood ratios (LR+ and LR-)
    thresh.selected$LR_pos <- ifelse(thresh.selected$fpr > 0,
                                     thresh.selected$dr / thresh.selected$fpr,
                                     NA)
    thresh.selected$LR_neg <-
      ifelse((100 - thresh.selected$fpr) > 0,
             (100 - thresh.selected$dr) / (100 - thresh.selected$fpr),
             NA
      )
    
    # calculate sensitivity and specificity
    sensitivity <- thresh.selected$dr / 100
    specificity <- 1 - (thresh.selected$fpr / 100)
    
    # calculate post-test probability
    post_test_prob <-
      sensitivity * pre_test_prob / ((sensitivity * pre_test_prob) + ((1 - specificity) * (1 - pre_test_prob)))
    
    # add positive pre- and post-test probabilities to results
    thresh.selected$pre_test_prob <- pre_test_prob
    thresh.selected$post_test_prob <- post_test_prob
    
    # add negative pre- and post-test probabilities to results
    thresh.selected$pre_test_prob_neg <- 1 - pre_test_prob
    thresh.selected$post_test_prob_neg <-
      (specificity * thresh.selected$pre_test_prob_neg) /
      ((specificity * thresh.selected$pre_test_prob_neg) + ((1 - sensitivity) * pre_test_prob))
    
    # add to cumulative results
    all_outcomes_results <-
      rbind(all_outcomes_results, thresh.selected)
  }
}


## replace outcome names and set custom order
all_outcomes_results$outcome <-
  new_labels[all_outcomes_results$outcome]
all_outcomes_results$outcome <-
  factor(all_outcomes_results$outcome, levels = custom_order)

all_outcomes_results$model <- factor(all_outcomes_results$model,
                                     levels = rev(c(
                                       "LOFO.mean", "ASCVD", "SCORE2", "MHO_Zembic", "AgeSexBMI"
                                     )))

## write to txt
fwrite(
  all_outcomes_results,
  "<path>/NNT_results.txt",
  row.names = F,
  na = NA
)