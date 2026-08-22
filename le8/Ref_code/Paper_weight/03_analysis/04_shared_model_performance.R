#############################################################
####                      Figure 4                       ####
#############################################################

# This script generates:
## Fig 4A : Comparison of shared model with outcome-specific models and established predictors
## Fig 4B : Features coefs heatmap
## Fig 4C : External validation EPIC-N
## Fig 4D : Fagan's nomograms for positive test (Type 2 Diabetes)
## Fig 4E : Fagan's nomograms for positive test (Chronic Renal Disease)
## SFig 11: Fagan's nomograms
## SFig 12: Fagan's nomograms
## ST 7   : C-Index comparison across models
## ST 3   : Case numbers

rm(list = ls())
setwd("<path>")

## --> packages needed <-- ##
library(data.table)
library(tidyverse)
library(dplyr)
library(magrittr)
library(ggrepel)
library(patchwork)
library(tidytext)
library(ggsci)
library(ggpubr)
library(openxlsx)

## read labs
lab.set <- fread("<path>")

## define outcomes
outc_list <- c("<outcomes>")

## define outcome labels
new_labels <- c("<nicenames>")

######################################################################
####     Comparison shared model with extended models (SFig 7)    ####
######################################################################
## list to store performances
cindex_data_list <- list()

## list to store top 20 variables
var_selection_simple <- list()
var_selection_extended <- list()

## loop through each outcome
for (outc in outc_list) {
  ## results from clinical features
  feat.sel.simple <-
    fread(paste0("<path>",
                 outc,
                 ".txt"))
  
  ## features to take forward
  feat.sel.simple <- feat.sel.simple[order(-sel)]
  jj <- grep("rand", feat.sel.simple$id)
  var.sel.simple <- feat.sel.simple[1:(jj[1] - 1)]
  
  ## results from extended features
  feat.sel.extended <-
    fread(paste0("<path>", outc, ".txt"))
  
  ## features to take forward
  feat.sel.extended <- feat.sel.extended[order(-sel)]
  jj <- grep("rand", feat.sel.extended$id)
  var.sel.extended <- feat.sel.extended[1:(jj[1] - 1)]
  
  ## load models
  feature.top.basic <-
    get(load(file = paste0("<path>",
                           outc,
                           ".RData")))
  feature.top.disdrug <-
    get(load(file = paste0("<path>",
                           outc,
                           ".RData")))
  feature.top.biomarker <-
    get(load(file = paste0("<path>",
                           outc,
                           ".RData")))
  feature.top.cardiopul <-
    get(load(file = paste0("<path>",
                           outc,
                           ".RData")))
  feature.top.bodycomp <-
    get(load(file = paste0("<path>",
                           outc,
                           ".RData")))
  feature.top.nmr <-
    get(load(file = paste0("<path>",
                           outc,
                           ".RData")))
  feature.top.pgs <-
    get(load(file = paste0("<path>",
                           outc,
                           ".RData")))
  
  ## prepare dataframe for performances
  cindex_data <- data.frame(
    Outcome = outc,
    FeatureSet = factor(
      c(
        "General health & lifestyle",
        "+Diseases & drug intake",
        "+Clinical blood biomarkers",
        "+Cardiopulmonary parameters",
        "+Body composition",
        "+Plasma metabolites",
        "+Polygenic scores"
      ),
      levels = c(
        "+Polygenic scores",
        "+Plasma metabolites",
        "+Body composition",
        "+Cardiopulmonary parameters",
        "+Clinical blood biomarkers",
        "+Diseases & drug intake",
        "General health & lifestyle"
      )
    ),
    Cindex.mean = c(
      feature.top.basic$Cindex.mean,
      feature.top.disdrug$Cindex.mean,
      feature.top.biomarker$Cindex.mean,
      feature.top.cardiopul$Cindex.mean,
      feature.top.bodycomp$Cindex.mean,
      feature.top.nmr$Cindex.mean,
      feature.top.pgs$Cindex.mean
    ),
    ci.low = c(
      feature.top.basic$ci.low,
      feature.top.disdrug$ci.low,
      feature.top.biomarker$ci.low,
      feature.top.cardiopul$ci.low,
      feature.top.bodycomp$ci.low,
      feature.top.nmr$ci.low,
      feature.top.pgs$ci.low
    ),
    ci.upp = c(
      feature.top.basic$ci.upp,
      feature.top.disdrug$ci.upp,
      feature.top.biomarker$ci.upp,
      feature.top.cardiopul$ci.upp,
      feature.top.bodycomp$ci.upp,
      feature.top.nmr$ci.upp,
      feature.top.pgs$ci.upp
    )
  )
  
  ## append to the list
  cindex_data_list[[outc]] <- cindex_data
  
  ## Prepare variable selection data
  lab.set <- fread("<path>")
  
  ## create color gradient for labels
  cat.col <- data.table(category = unique(c(lab.set$category)),
                        cl = colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(16))
  cat.col.vector <- setNames(cat.col$cl, cat.col$category)
  
  ## variable selection simple
  var.sel.simple <-
    left_join(var.sel.simple[1:20, c(1, 253)], lab.set[, c(4, 6, 9)], by =
                c("id" = "short_name_new"))
  
  var.sel.simple <- var.sel.simple[, Outcome := outc]
  var.sel.simple <- var.sel.simple[1:20]
  
  ## variable selection extended
  var.sel.extended <-
    left_join(var.sel.extended[1:20, c(1, 253)], lab.set[, c(4, 6, 9)], by =
                c("id" = "short_name_new"))
  
  var.sel.extended <- var.sel.extended[, Outcome := outc]
  var.sel.extended <- var.sel.extended[1:20]
  
  ## append to list
  var_selection_simple[[outc]] <- var.sel.simple
  var_selection_extended[[outc]] <- var.sel.extended
  
}

## combine all C-index data
cindex_data_all <- as.data.table(bind_rows(cindex_data_list))

## replace outcome names
cindex_data_all$Outcome <- new_labels[cindex_data_all$Outcome]

## order outcomes based on mean C-index
ordered_outcomes <-
  cindex_data_all[FeatureSet == "+Polygenic scores"][order(-Cindex.mean), Outcome]

## reorder facets based on performance
cindex_data_all$Outcome <-
  factor(cindex_data_all$Outcome, levels = unique(ordered_outcomes))

## prepare dataset for ribbon
cindex_data_ribbon <- cindex_data_all %>%
  group_by(Outcome) %>%
  mutate(ymin = min(Cindex.mean), ymax = max(Cindex.mean))

## change order of feature set
cindex_data_all$FeatureSet <- factor(cindex_data_all$FeatureSet,
                                     levels = rev(levels(cindex_data_all$FeatureSet)))

## load results from shared model
all_res <- fread("<path>")

## restrict outcome-specific models to the extended best possible models
cindex_data_full <-
  cindex_data_all %>% filter(FeatureSet %in% "+Polygenic scores") %>%
  dplyr::select(Outcome, Cindex.mean, ci.low, ci.upp) %>%
  mutate(Model = "Extended outcome-wise models") %>%
  mutate(Nfeatures = 20)

## join shared model performances with outcome-specific models
al_res <- rbind(all_res, cindex_data_full, fill = T)

## new order
ordered_levels <-
  al_res$Outcome[al_res$Model == "LOFO20_Mean"][order(-al_res$Cindex.mean[al_res$Model == "LOFO20_Mean"])]

## reorder
al_res$Outcome <- factor(al_res$Outcome, levels = ordered_levels)
al_res$Models <- "OBSCORE vs.\noutcome-specific models"

## prepare ribbon dataset to compare
al_res_ribbon <- al_res %>%
  group_by(Outcome) %>%
  mutate(xmin = min(Cindex.mean), xmax = max(Cindex.mean))

## plot in supplements different version instead during revision
a1 <-
  ggplot(al_res, aes(x = Cindex.mean, y = Outcome, fill = as.factor(Model))) +
  facet_grid(Outcome ~ Models, scales = "free") +
  geom_rect(
    data = al_res_ribbon,
    aes(
      ymin = -Inf,
      ymax = Inf,
      xmin = xmin,
      xmax = xmax,
      group = Outcome
    ),
    fill = "#C52C4B80",
    alpha = 0.25
  ) +
  geom_pointrange(
    aes(
      shape = as.factor(Model),
      xmin = ci.low,
      xmax = ci.upp
    ),
    position = position_dodge(width = 1),
    linewidth = 0.25,
    size = .3,
    stroke = .5
  ) +
  geom_text(
    aes(label = round(Cindex.mean, 2), color = as.factor(Model)),
    position = position_dodge(width = 1),
    hjust = 1.5,
    size = 2,
    color = "black"
  ) +
  labs(y = "", x = "C - Index (95% CI)") +
  scale_fill_manual(values = c("#EE4C97FF", "#0072B5FF")) +
  scale_shape_manual(values = c(21, 23)) +
  theme_light() +
  xlim(0.4, 1) +
  geom_vline(
    xintercept = 0.5,
    linetype = "dotted",
    color = "black",
    size = 0.5
  ) +
  theme(
    legend.title = element_blank(),
    legend.position = c(.805, 0.15),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    strip.text.y = element_blank(),
    strip.background.y = element_blank(),
    strip.text.x = element_text(color = "white", angle = 0),
    strip.background = element_rect(
      color = "black",
      fill = "black",
      linetype = "solid"
    ),
    panel.spacing = unit(0.05, "lines")
  )

### ---- metrics for the text ---- ###
## delta OBSCORE and extended outcome-specific models
round(median((
  al_res[Model == "LOFO20_Mean"] %>%
    left_join(al_res[Model == "Extended outcome-wise models"], by =
                "Outcome") %>%
    mutate(delta = Cindex.mean.y - Cindex.mean.x)
)$delta), 3)
round(range((
  al_res[Model == "LOFO20_Mean"] %>%
    left_join(al_res[Model == "Extended outcome-wise models"], by = "Outcome") %>%
    mutate(delta = Cindex.mean.y - Cindex.mean.x)
)$delta), 3)
# 0.012 (-0.006 - 0.024)

#########################################################################
####    Comparison shared model with established models (Fig 4A)    ####
#########################################################################

AgeSexBMI_res <- fread("<path>")
MHO_Zembic_res <- fread("<path>")
ASCVD_res <- fread("<path>")
SCORE2_res <- fread("<path>")

comparison_res <-
  rbind(ASCVD_res, MHO_Zembic_res, SCORE2_res, AgeSexBMI_res, all_res)
comparison_res$Models <- "OBSCORE vs.\nestablished predictors"

## reorder the outcomes based on the maximum Cindex.mean.y
comparison_res <- comparison_res %>%
  mutate(Outcome = fct_reorder(Outcome, Cindex.mean, .desc = TRUE))

## reorder
outcome_order <-
  c("AgeSexBMI",
    "MetS_Grundy",
    "MHO_Zembic",
    "SCORE2",
    "ASCVD",
    "LOFO20_Mean")

## set factor levels
comparison_res <- comparison_res %>%
  mutate(Model = factor(Model, levels = outcome_order))

comparison_ribbon <- comparison_res %>%
  group_by(Outcome) %>%
  mutate(xmin = min(Cindex.mean), xmax = max(Cindex.mean))

ordered_levels <-
  comparison_res$Outcome[comparison_res$Model == "LOFO20_Mean"][order(-comparison_res$Cindex.mean[comparison_res$Model == "LOFO20_Mean"])]

## apply to all rows
comparison_res$Outcome <-
  factor(comparison_res$Outcome, levels = ordered_levels)

## plot with reordered facets
a2 <- comparison_res %>%
  ggplot() +
  facet_grid(Outcome ~ Models, scales = "free") +
  geom_rect(
    data = comparison_ribbon,
    aes(
      ymin = -Inf,
      ymax = Inf,
      xmin = xmin,
      xmax = xmax,
      group = Outcome
    ),
    fill = "#0099B4FF",
    alpha = 0.025
  ) +
  geom_pointrange(
    aes(
      fill = as.factor(Model),
      y = Model,
      x = Cindex.mean,
      xmin = ci.low,
      xmax = ci.upp,
      shape = as.factor(Model)
    ),
    linewidth = 0.25,
    size = .3,
    stroke = .5
  ) +
  labs(y = "", x = "C - Index (95% CI)") +
  scale_fill_manual(values = c(
    "#FFDC91FF",
    "#20854EFF",
    "#7876B1FF",
    "#E18727FF",
    "#0072B5FF"
  )) +
  scale_shape_manual(values = c(21, 21, 21, 21, 23), guide = "none") +
  theme_light() +
  xlim(0.4, 1) +
  geom_vline(
    xintercept = 0.5,
    linetype = "dotted",
    color = "black",
    size = 0.5
  ) +
  scale_y_discrete(expand = expansion(mult = .35)) +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 5),
    legend.position = c(.805, 0.15),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    strip.text.y = element_text(
      color = "white",
      angle = 0,
      hjust = 0
    ),
    strip.text.x = element_text(color = "white", angle = 0),
    strip.background = element_rect(
      color = "black",
      fill = "black",
      linetype = "solid"
    ),
    panel.spacing = unit(0.05, "lines")
  ) +
  guides(fill = guide_legend(reverse = TRUE),
         shape = guide_legend(reverse = TRUE))

### ---- metrics for the text ---- ###
## median (range) delta OBSCORE ASCVD
round(median((
  comparison_res[Model == "ASCVD"] %>%
    left_join(comparison_res[Model == "LOFO20_Mean"], by = "Outcome") %>%
    mutate(delta = Cindex.mean.y - Cindex.mean.x) %>% arrange(desc(delta))
)$delta), 3)
round(range((
  comparison_res[Model == "ASCVD"] %>%
    left_join(comparison_res[Model == "LOFO20_Mean"], by = "Outcome") %>%
    mutate(delta = Cindex.mean.y - Cindex.mean.x) %>% arrange(desc(delta))
)$delta), 3)
# 0.027 (-0.004 - 0.146)

##############################################################
####    Shared model features and importances (Fig 4B)    ####
##############################################################
## load LOFO feature importances
LOFO_res <- fread("<path>")

## attach proper feature names
LOFO_res %<>% left_join(., lab.set[, c(6, 9:10)],
                        by = c("ExcludedFeature" = "short_name_new"))

## some renaming
LOFO_res %<>%
  mutate(
    label_new = ifelse(
      label_new == "Long-standing illness, disability or infirmity Yes",
      "Long-standing illness",
      label_new
    )
  ) %>%
  mutate(label_new = ifelse(label_new == "Glycated haemoglobin (HbA1c)",
                            "HbA1c",
                            label_new)) %>%
  mutate(label_new = ifelse(label_new == "Sex Female",
                            "Sex",
                            label_new))

## get mean feature importance across outcomes, pick top 20 features with highest mean
top <- LOFO_res %>% filter(ExcludedFeature != "None") %>%
  group_by(ExcludedFeature) %>%
  summarise(mean_Cindex_diff = mean(Cindex.diff, na.rm = TRUE)) %>%
  arrange(desc(-mean_Cindex_diff)) %>% slice(1:20) %>% pull(ExcludedFeature)

## get relative importance based on highest
featureimportance <-
  LOFO_res %>% filter(ExcludedFeature %in% top) %>%
  group_by(label_new) %>%
  summarise(mean_Cindex_diff = mean(Cindex.diff, na.rm = TRUE)) %>%
  mutate(facet = "Mean importance") %>%
  mutate(Importance = mean_Cindex_diff / -0.00861) %>% dplyr::select(label_new, facet, Importance)

## order based on importance
ordering <- featureimportance %>%
  filter(facet == "Mean importance") %>%
  arrange(desc(-Importance)) %>%
  pull(label_new)

## create this for facet
featureimportance$facet_ <- "Feature importance\nacross outcomes"

## function for broken axis
trans <- function(x) {
  pmin(x, 0.575) + 0.05 * pmax(x - 0.6, 0)
}

## define y-axis breaks and labels
xticks <- c(0, 0.2, 0.4, 0.6, 1)

## transform data onto the display scale
featureimportance$Importance_t <-
  trans(featureimportance$Importance)

## plot deprecated during revision
b <- featureimportance %>% filter(facet == "Mean importance") %>%
  ggplot(aes(x = facet,
             y = factor(label_new, levels = ordering))) +
  facet_grid( ~ facet_) +
  geom_point(aes(size = Importance, ),
             fill = "#0072B5FF",
             shape = 22,
  ) +
  labs(y = "", x = "") +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 5),
    legend.key.size = unit(0.1, 'cm'),
    legend.key.spacing = unit(0, 'cm'),
    legend.margin = margin(0, 10, 0, 0),
    legend.box.margin = margin(-10,-10,-10,-10),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.y = element_text(angle = 0, hjust = 1),
    axis.text.x = element_blank(),
    strip.text.y = element_text(color = "white", angle = 0),
    strip.text.x = element_text(color = "white", angle = 0),
    strip.background = element_rect(
      color = "black",
      fill = "black",
      linetype = "solid"
    ),
    panel.spacing = unit(0.2, "lines")
  )

## list to store coefficients
coef_list <- list()
## read coefficients
for (outc in outc_list) {
  model_name <- paste0("<path>", outc)
  model <- get(load(paste0("<path>/", model_name, ".RData")))
  
  coefs <- model$opt.coefficients
  
  coef_dt <- data.table(variable = rownames(coefs),
                        value = as.numeric(coefs))
  setnames(coef_dt, "value", outc)
  coef_list[[outc]] <- coef_dt
}

## merge based on features
coef_merged <-
  Reduce(function(x, y)
    merge(x, y, by = "variable", sort = FALSE), coef_list)

## write coefs
write.xlsx(coef_merged, "<path>/OBSCORE_Coefficients.xlsx")

## revised Figure 4B
coef_long <-
  melt(
    coef_merged,
    id.vars = "variable",
    variable.name = "Outcome",
    value.name = "Coefficient"
  )

coef_long %<>% left_join(., lab.set[, c(6, 9:10)],
                         by = c("variable" = "short_name_new"))

coef_long$Outcome <- as.character(coef_long$Outcome)

coef_long$Outcome <- new_labels[coef_long$Outcome]


coef_long$label_new <- factor(coef_long$label_new,
                              levels = coef_long$label_new[order(rowMeans(as.matrix(coef_merged[,-1, with = FALSE]), na.rm = TRUE))])

## new plot after revision
b2 <-
  ggplot(coef_long, aes(
    x = factor(Outcome, levels = ordered_levels),
    y = factor(label_new, levels = ordering),
    fill = Coefficient
  )) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "#4dbbd0",
    mid = "white",
    high = "#e64b00",
    midpoint = 0
  ) +
  theme_light() +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = .5
    ),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_rect(
      color = "black",
      fill = "black",
      linetype = "solid"
    )
  ) +
  labs(fill = "Coefficients")

##############################################################
####          External validation EPIC-N (Fig 4C)         ####
##############################################################
## OBSCORE results (18 features)
obscore_light <- fread("<path>")
## EPIC-N results (18 features)
epic_n <- fread("<path>")
## merge internal & external validation
validation <-
  obscore_light %>% left_join(., epic_n, by = c("Outcome" = "disease.label"))

## plot validation
validation_plot <-
  validation %>%
  filter(n.cases > 20) %>%
  filter(!is.na(obscore.cindex.mean)) %>%
  ggplot(aes(Cindex.mean, obscore.cindex.mean)) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "black"
  ) +
  geom_smooth(
    method = "lm",
    alpha = .2,
    color = "grey80",
    fill = "grey80"
  ) +
  geom_errorbarh(
    aes(xmin = ci.low, xmax = ci.upp, fill = Outcome),
    linewidth = 0.25,
    shape = 21,
    size = 0.35
  ) +
  geom_errorbar(
    aes(ymin = obscore.cindex.ci.lo, ymax = obscore.cindex.ci.upp, fill = Outcome),
    linewidth = 0.25,
    shape = 21,
    size = 0.35
  ) +
  geom_point(
    shape = 21,
    fill = "skyblue",
    size = 2,
    alpha = .9
  ) +
  geom_text_repel(
    aes(label = Outcome),
    size = 1.5,
    hjust = 1,
    vjust = 1,
    max.overlaps = Inf,
    segment.size = 0.3,
    segment.color = "gray70"
  ) +
  stat_cor(method = "pearson", size = 2) +
  labs(x = "C-index (95% CI) UK Biobank", y = "C-index (95% CI) EPIC-Norfolk") +
  theme_light() +
  theme(
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = "none"
  )

######################################################################
####     Fagan's nomograms for T2D and Kidney disease (Fig 4)     ####
######################################################################
all_outcomes_results <- fread("<path>")

## calculate incidence (prevalence), odds and log odds
dd_line <-
  all_outcomes_results %>%
  mutate(
    prevalence = pre_test_prob * 100,
    pre_test_odds = pre_test_prob / (1 - pre_test_prob),
    post_test_odds = post_test_prob / (1 - post_test_prob),
    log_pre_test_odds = log10(pre_test_odds),
    log_post_test_odds = log10(post_test_odds)
  ) %>%
  gather(key = "category",
         value = "log_odds",
         log_pre_test_odds,
         log_post_test_odds) %>%
  mutate(x = ifelse(category == "log_pre_test_odds", 0, 1))  ## x-axis positions

## transform dataset for plotting
dd_line_transformed <- dd_line %>%
  rename(disease = "outcome",
         model = "model",
         FPR = "FPR.true") %>%
  dplyr::select(
    disease,
    model,
    prevalence,
    FPR,
    dr,
    LR_pos,
    pre_test_prob,
    pre_test_odds,
    post_test_prob
  )

## rename
dat <- dd_line_transformed

## defined functions
odds         <- function(p) {
  # function converts probability into odds
  o <- p / (1 - p)
  return(o)
}

logodds      <- function(p) {
  # function returns logodds for a probability
  lo <- log10(p / (1 - p))
  return(lo)
}

logodds_to_p <- function(lo) {
  # function goes from logodds back to a probability
  o <- 10 ^ lo
  p <- o / (1 + o)
  return(p)
}

p2percent <- function(p) {
  # function turns numeric probability into string percentage
  scales::percent(signif(p, digits = 3))
}

## define axis ticks
ticks_prob    <- c(1, 2, 5, 10, 20, 30,
                   40, 50, 60, 70, 80, 90, 95, 99)

## convert % to odds
ticks_odds    <- odds(ticks_prob / 100)

## convert % to logodds
ticks_logodds <- logodds(ticks_prob / 100)

## select the likelihood ratios of interest (for the middle y-axis)
ticks_lrs     <- c(0.01, 0.05, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100)

## log10 them since plot is in logodds space
ticks_log_lrs <- log10(ticks_lrs)

## fixing particular x-coordinates
left     <- 0
right    <- 1
middle   <- 0.5
midright <- 0.75

## initially these are expressed as probabilities
df <- dat %>%
  gather(key = "category", value = "prob", c(7, 9)) %>%
  mutate(x = ifelse(category == "pre_test_prob" , left, right)) #pre_test_odds

adj_min      <- range(ticks_logodds)[1]
adj_max      <- range(ticks_logodds)[2]
adj_diff     <- adj_max - adj_min
scale_factor <- abs(adj_min) - adj_diff / 2

## convert probabilities to logodds for plotting
df$lo_y  <-
  ifelse(df$x == left,
         logodds(1 - df$prob) - scale_factor,
         logodds(df$prob))

rescale   <- range(ticks_logodds) + abs(adj_min) - adj_diff / 2
rescale_x_breaks  <- ticks_logodds + abs(adj_min) - adj_diff / 2

## refactor levels
df$model <- factor(df$model,
                   levels = rev(c(
                     "LOFO.mean", "ASCVD", "SCORE2", "MHO_Zembic", "AgeSexBMI"
                   )))

## plot for T2D
d <- df %>% filter(disease %in% c("Type 2 Diabetes")) %>%
  ggplot() +
  geom_line(aes(x = x, y = lo_y, color = model), size = .3) +
  geom_vline(xintercept = middle,
             linewidth = .3,
             color = "grey") +
  annotate(
    geom = "text",
    x = rep(middle + .075, length(ticks_log_lrs)),
    y = (ticks_log_lrs - scale_factor) / 2,
    label = ticks_lrs,
    size = rel(10 / 5)
  ) +
  annotate(
    geom = "point",
    x = rep(middle, length(ticks_log_lrs)),
    y = (ticks_log_lrs - scale_factor) / 2,
    size = 0.5,
    color = "grey"
  ) +
  scale_color_manual(values = c(
    "#FFDC91FF",
    "#20854EFF",
    "#7876B1FF",
    "#E18727FF",
    "#0072B5FF"
  )) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = rescale,
    breaks = -rescale_x_breaks,
    labels = ticks_prob,
    name = "Pre-test probability",
    sec.axis = sec_axis(
      transform = ~ .,
      name = "Post-test probability",
      labels = ticks_prob,
      breaks = ticks_logodds
    )
  ) +
  theme_light() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(angle = 90),
    axis.title.y.right = element_text(angle = 90),
    axis.line = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "none"
  )

## plot for chronic renal disease
e <- df %>% filter(disease %in% c("Chronic Renal Failure")) %>%
  ggplot() +
  geom_line(aes(x = x, y = lo_y, color = model), size = .3) +
  geom_vline(xintercept = middle,
             linewidth = .3,
             color = "grey") +
  annotate(
    geom = "text",
    x = rep(middle + .075, length(ticks_log_lrs)),
    y = (ticks_log_lrs - scale_factor) / 2,
    label = ticks_lrs,
    size = rel(10 / 5)
  ) +
  annotate(
    geom = "point",
    x = rep(middle, length(ticks_log_lrs)),
    y = (ticks_log_lrs - scale_factor) / 2,
    size = 0.5,
    color = "grey"
  ) +
  scale_color_manual(values = c(
    "#FFDC91FF",
    "#20854EFF",
    "#7876B1FF",
    "#E18727FF",
    "#0072B5FF"
  )) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = rescale,
    breaks = -rescale_x_breaks,
    labels = ticks_prob,
    name = "Pre-test probability",
    sec.axis = sec_axis(
      transform = ~ .,
      name = "Post-test probability",
      labels = ticks_prob,
      breaks = ticks_logodds
    )
  ) +
  theme_light() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(angle = 90),
    axis.title.y.right = element_text(angle = 90),
    axis.line = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "none"
  )

layout = "
AAAAABBB
"

## plot and save
pdf("<path>/Fig_4_1.pdf", height = 5, width = 8.75)
(a2 + b2 + plot_layout(design = layout) &
    theme(text = element_text(size = 8)))
dev.off()

layout = "
AABC
"

## plot and save
pdf("<path>/Fig_4_2.pdf", height = 3, width = 7.5)
validation_plot + d + e + plot_layout(design = layout) &
  theme(text = element_text(size = 7))
dev.off()

## function to plot nomograms for all diseases
plot_disease <- function(df, diseases) {
  plots <- lapply(diseases, function(d) {
    df %>%
      filter(disease == d) %>%
      ggplot() +
      geom_line(aes(x = x, y = lo_y, color = model), size = .3) +
      geom_vline(xintercept = middle,
                 linewidth = .3,
                 color = "grey") +
      annotate(
        geom = "text",
        x = rep(middle + .075, length(ticks_log_lrs)),
        y = (ticks_log_lrs - scale_factor) / 2,
        label = ticks_lrs,
        size = rel(10 / 5)
      ) +
      annotate(
        geom = "point",
        x = rep(middle, length(ticks_log_lrs)),
        y = (ticks_log_lrs - scale_factor) / 2,
        size = 0.5,
        color = "grey"
      ) +
      scale_color_manual(values = c(
        "#FFDC91FF",
        "#20854EFF",
        "#7876B1FF",
        "#E18727FF",
        "#0072B5FF"
      )) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(
        expand = c(0, 0),
        limits = rescale,
        breaks = -rescale_x_breaks,
        labels = ticks_prob,
        name = "Pre-test probability",
        sec.axis = sec_axis(
          transform = ~ .,
          name = "Post-test probability",
          labels = ticks_prob,
          breaks = ticks_logodds
        )
      ) +
      theme_light() +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(angle = 90),
        axis.title.y.right = element_text(angle = 90),
        axis.line = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "right"
      ) +
      ggtitle(d)
  })
  
  ## combine plots
  final_plot <-
    wrap_plots(plots) & theme(text = element_text(size = 7))
  return(final_plot)
}


diseases_to_plot <-
  c(
    "Chronic Renal Failure",
    "Gout",
    "Type 2 Diabetes",
    "CV Death",
    "Sleep Apnea",
    "NAFLD",
    "Hypercholesterolemia",
    "Ischemic Heart",
    "MI"
  )

## plot and save
pdf("<path>/FaganNomograms1.pdf",
    width = 7,
    height = 7)
plot_disease(df, diseases_to_plot)
dev.off()

diseases_to_plot <-
  c(
    "MACE",
    "Stroke",
    "Coronary Athero",
    "Angina Pectoris",
    "Hypertension",
    "Cholelithiasis",
    "Arthropathy",
    "Diaph Hernia",
    "GERD"
  )

## plot and save
pdf("<path>/FaganNomograms2.pdf",
    width = 7,
    height = 7)
plot_disease(df, diseases_to_plot)
dev.off()

## combine datasets into one
combined_res <-
  rbind(all_res, ASCVD_res, SCORE2_res, MHO_Zembic_res, AgeSexBMI_res)[, c(1:4, 6)]

## create a column with mean Cindex and confidence intervals
combined_res[, Cindex_combined := sprintf("%.3f (%.3f - %.3f)", Cindex.mean, ci.low, ci.upp)]

## order models in table
model_order <-
  c("LOFO20_Mean", "ASCVD", "SCORE2", "MHO_Zembic", "AgeSexBMI")
combined_res[, Model := factor(Model, levels = model_order)]

## reshape to wide for displaying
wide_table <-
  dcast(combined_res, Outcome ~ Model, value.var = "Cindex_combined")

## save table
fwrite(wide_table, "<path>/Comparison_CIndex.csv")

######################################################################
####    Scatter plot comparing shared and extended (Supp Fig 6)   ####
######################################################################
## plot and save
pdf("<path>/SuppFig_6.pdf",
    height = 5,
    width = 6)
all_res %>%
  left_join(cindex_data_full, by = "Outcome") %>%
  ggplot(aes(Cindex.mean.x, Cindex.mean.y)) +
  geom_abline() +
  geom_abline(intercept = 0.025,
              color = "grey",
              linetype = "dashed") +
  geom_abline(intercept = -0.025,
              color = "grey",
              linetype = "dashed") +
  geom_ribbon(
    data = data.frame(Cindex.mean.x = seq(x_limits[1], x_limits[2], length.out = 100)),
    aes(
      x = Cindex.mean.x,
      ymin = Cindex.mean.x - 0.025,
      ymax = Cindex.mean.x + 0.025
    ),
    inherit.aes = FALSE,
    fill = "grey80",
    alpha = 0.5
  ) +
  geom_pointrange(
    aes(
      fill = as.factor(Outcome),
      xmin = ci.low.x,
      xmax = ci.upp.x
    ),
    linewidth = 0.25,
    shape = 21,
    size = .4,
    alpha = .7
  ) +
  geom_pointrange(
    aes(
      fill = as.factor(Outcome),
      ymin = ci.low.y,
      ymax = ci.upp.y
    ),
    linewidth = 0.25,
    shape = 21,
    size = .4,
    alpha = .7
  ) +
  labs(x = "C - Index (95% CI) shared unified clinical model", y = "C - Index (95% CI) outcome-specific extended models") +
  theme_light() +
  xlim(c(0.4, 1)) +
  ylim(c(0.4, 1)) +
  geom_vline(
    xintercept = 0.5,
    linetype = "dotted",
    color = "black",
    size = 0.5
  ) +
  geom_hline(
    yintercept = 0.5,
    linetype = "dotted",
    color = "black",
    size = 0.5
  ) +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 5),
    legend.key.size = unit(0.5, 'cm'),
    legend.key.spacing = unit(0, 'cm'),
    legend.margin = margin(0, 10, 0, 0),
    legend.box.margin = margin(-10,-10,-10,-10),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    strip.text.y = element_text(
      color = "white",
      angle = 0,
      hjust = 0
    ),
    strip.text.x = element_text(color = "white", angle = 0),
    strip.background = element_rect(
      color = "black",
      fill = "black",
      linetype = "solid"
    ),
    panel.spacing = unit(0.2, "lines")
  ) & theme(text = element_text(size = 9))
dev.off()

######################################################################
####     Top 20 clinical and extended features (Supp Fig 4,5)     ####
######################################################################

## combine all simple variables
var_selection_simple <- bind_rows(var_selection_simple)
## replace the outcomes
var_selection_simple$Outcome <-
  new_labels[var_selection_simple$Outcome]

## combine all extended variables
var_selection_extended <- bind_rows(var_selection_extended)
## replace the outcomes
var_selection_extended$Outcome <-
  new_labels[var_selection_extended$Outcome]

var_selection_simple %<>% left_join(., lab.set[, 9:10], by = c("id" = "short_name_new"))
var_selection_extended %<>% left_join(., lab.set[, 9:10], by = c("id" =
                                                                   "short_name_new"))

## some renaming
var_selection_simple <- var_selection_simple %>%
  mutate(
    label_new = ifelse(
      label_new == "Caffeine |Acetaminophen |Codeine|Acetaminophen |Codeine |Ketorolac|Acetaminophen |Codeine",
      "Pain killer",
      label_new
    )
  )

## some renaming
var_selection_extended <- var_selection_extended %>%
  mutate(
    label_new = ifelse(
      label_new == "Caffeine |Acetaminophen |Codeine|Acetaminophen |Codeine |Ketorolac|Acetaminophen |Codeine",
      "Pain killer",
      label_new
    )
  )

## calculate the number of unique outcomes each ID appears in
id_order_simple <- var_selection_simple %>%
  group_by(label_new) %>%
  summarise(outcome_count = n_distinct(Outcome)) %>%
  arrange(desc(-outcome_count)) %>%
  pull(label_new)

## calculate the number of unique outcomes each ID appears in
id_order_extended <- var_selection_extended %>%
  group_by(label_new) %>%
  summarise(outcome_count = n_distinct(Outcome)) %>%
  arrange(desc(-outcome_count)) %>%
  pull(label_new)

## clinical features
a <-
  var_selection_simple %>% filter(label_new %in% id_order_simple) %>%
  mutate(label_new = factor(label_new, levels = id_order_simple)) %>%
  ggplot(aes(x = Outcome, y = label_new, fill = category)) +
  geom_tile(shape = 21,
            size = .4,
            color = "black") +
  labs(x = "", y = "") +
  theme_light() +
  coord_fixed() +
  scale_x_discrete(position = "top") +
  scale_fill_manual(
    values = cat.col.vector,
    guide = guide_legend(
      ncol = 1,
      byrow = TRUE,
      keyheight = unit(0.2, "cm"),
      keywidth = unit(0.2, "cm")
    )
  ) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 6),
    legend.title = element_blank(),
    axis.text.x = element_text(
      angle = 90,
      hjust = 0,
      vjust = 0.5
    ),
    axis.text = element_text(size = 5),
    strip.text = element_text(
      color = "white",
      angle = 90,
      hjust = 0
    ),
    strip.background = element_rect(
      color = "black",
      fill = "black",
      linetype = "solid"
    ),
    panel.spacing = unit(0.1, "lines"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    text = element_text(size = 5)
  )

## extended features
b <-
  var_selection_extended %>% filter(label_new %in% id_order_extended) %>%
  mutate(label_new = factor(label_new, levels = id_order_extended)) %>%
  ggplot(aes(x = Outcome, y = label_new, fill = category)) +
  geom_tile(shape = 21,
            size = .4,
            color = "black") +
  labs(x = "", y = "") +
  theme_light() +
  coord_fixed() +
  scale_x_discrete(position = "top") +
  scale_fill_manual(
    values = cat.col.vector,
    guide = guide_legend(
      ncol = 1,
      byrow = TRUE,
      keyheight = unit(0.2, "cm"),
      keywidth = unit(0.2, "cm")
    )
  ) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 6),
    legend.title = element_blank(),
    axis.text.x = element_text(
      angle = 90,
      hjust = 0,
      vjust = .5
    ),
    axis.text = element_text(size = 5),
    strip.text = element_text(
      color = "white",
      angle = 90,
      hjust = 0
    ),
    strip.background = element_rect(
      color = "black",
      fill = "black",
      linetype = "solid"
    ),
    panel.spacing = unit(0.1, "lines"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    text = element_text(size = 5)
  )

## plot and save clinical variable selection
pdf("<path>/SuppFig_4.pdf",
    height = 10 ,
    width = 6)
a & theme(text = element_text(size = 5))
dev.off()

## plot and save extended variable selection
pdf("<path>/SuppFig_5.pdf",
    height = 10 ,
    width = 6)
b & theme(text = element_text(size = 5))
dev.off()

###############################################################
####  Supplementary Table for number of individuals (ST3)  ####
###############################################################
## list to store results
case_numbers_list <- list()

## load case numbers for each outcomes
for (outcome in outc_list) {
  file_path <-
    paste0("<path>", outcome, ".txt")
  if (file.exists(file_path)) {
    case_data <- fread(file_path)
    case_data[, Outcome := outcome]
    case_numbers_list[[outcome]] <- case_data
  } else {
    warning(paste("File not found for outcome:", outcome))
  }
}

## combine all
combined_case_numbers <-
  rbindlist(case_numbers_list, use.names = TRUE, fill = TRUE)

## replace names
combined_case_numbers$Outcome <-
  new_labels[combined_case_numbers$Outcome]

## save
output_file <- "<path>/combined_case_numbers.csv"
fwrite(combined_case_numbers, output_file)