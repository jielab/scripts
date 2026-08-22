#############################################################
####                 Unimodal performances               ####
#############################################################

# This script generates:
## Fig 2 : Performances of unimodal models
## SFig 2: A detailed plot with performances of unimodal models
## SFig 3: Comparison of unimodal models with 20 features vs. all non-random features
## ST 4  : A table with unimodal performances

rm(list = ls())
setwd("<path>")

## --> packages needed <-- ##
library(data.table)
library(tidyverse)
library(readxl)
library(dplyr)
library(magrittr)
library(ggrepel)
library(patchwork)
library(ggsci)

## define outcomes
outc_list <- c("<outcomes>")

## define outcome labels
new_labels <- c("<nicenames>")

#############################################
####    Unimodal performances (Fig 2)    ####
#############################################
## list to store results
cindex_individual_list <- list()

## loop through each outcome
for (outc in outc_list) {
  ## load models
  feature.top.basic <-
    get(load(file = paste0("<path>", outc, ".RData")))
  feature.top.disdrug <-
    get(load(file = paste0("<path>", outc, ".RData")))
  feature.top.biomarker <-
    get(load(file = paste0("<path>", outc, ".RData")))
  feature.top.cardiopul <-
    get(load(file = paste0("<path>", outc, ".RData")))
  feature.top.bodycomp <-
    get(load(file = paste0("<path>", outc, ".RData")))
  feature.top.nmr <-
    get(load(file = paste0("<path>", outc, ".RData")))
  feature.top.pgs <-
    get(load(file = paste0("<path>", outc, ".RData")))
  
  ## create dataframe of performances
  cindex_data <- data.frame(
    Outcome = outc,
    FeatureSet = factor(
      c(
        "General health & lifestyle",
        "Diseases & drug intake",
        "Clinical blood biomarkers",
        "Cardiopulmonary parameters",
        "Body composition",
        "Plasma metabolites",
        "Polygenic scores"
      ),
      levels = c(
        "Polygenic scores",
        "Plasma metabolites",
        "Body composition",
        "Cardiopulmonary parameters",
        "Clinical blood biomarkers",
        "Diseases & drug intake",
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
  cindex_individual_list[[outc]] <- cindex_data
  
}
## combine all C-index data
cindex_individual_all <-
  as.data.table(bind_rows(cindex_individual_list))

## replace outcome names
cindex_individual_all$Outcome <-
  new_labels[cindex_individual_all$Outcome]

## define plot order based on decreasing mean C-index
ordered_outcomes <-
  cindex_individual_all[FeatureSet == "Polygenic scores"][order(-Cindex.mean), Outcome]

## reorder outcome order
cindex_individual_all$Outcome <-
  factor(cindex_individual_all$Outcome, levels = unique(ordered_outcomes))

## prepare ribbon shading dataset
cindex_individual_ribbon <- cindex_individual_all %>%
  group_by(Outcome) %>%
  mutate(ymin = min(Cindex.mean), ymax = max(Cindex.mean))

## create a dummy column for facet
cindex_individual_all$Features <- "Performance per outcome"

## get colors from the nejm palette and reorder them
nejm_colors <- pal_nejm("default")(8)
reordered_colors <-
  nejm_colors[c(7, 6, 5, 4, 3, 2, 1)] # Change the indices as needed

cindex_individual_all$FeatureSet <-
  factor(cindex_individual_all$FeatureSet,
         levels = rev(levels(cindex_individual_all$FeatureSet)))

## plot
a <-
  ggplot(cindex_individual_all, aes(x = Cindex.mean, y = FeatureSet)) +
  facet_grid(FeatureSet ~ Features, scales = "free") +
  geom_boxplot(aes(fill = as.factor(FeatureSet)),
               alpha = 0.7,
               outlier.color = NA) +
  geom_pointrange(
    aes(
      fill = as.factor(FeatureSet),
      xmin = ci.low,
      xmax = ci.upp
    ),
    position = position_dodge(width = 1),
    linewidth = 0.25,
    shape = 21,
    size = .35,
    alpha = .7
  ) +
  labs(y = "", x = "C - Index (95% CI)") +
  scale_fill_nejm() +
  theme_light() +
  xlim(0.4, 1) +
  geom_vline(
    xintercept = 0.5,
    linetype = "dotted",
    color = "black",
    size = 0.5
  ) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 5),
    legend.key.size = unit(0.5, 'cm'),
    legend.key.spacing = unit(0, 'cm'),
    legend.margin = margin(0, 10, 0, 0),
    legend.box.margin = margin(-10, -10, -10, -10),
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
  )

## save
pdf("<path>/Fig_2.pdf", height = 3.5, width = 7)
a
dev.off()

#############################################
####      Detailed unimodal (SFig 2)     ####
#############################################
## refactor modalities
cindex_individual_all$FeatureSet <-
  factor(cindex_individual_all$FeatureSet,
         levels = rev(levels(cindex_individual_all$FeatureSet)))

## plot extended unimodal performances
b <-
  ggplot(cindex_individual_all, aes(x = FeatureSet, y = Cindex.mean)) +
  coord_flip() +
  facet_wrap( ~ Outcome, ncol = 3) +
  geom_pointrange(
    aes(ymin = ci.low, ymax = ci.upp, fill = FeatureSet),
    linewidth = 0.25,
    shape = 21,
    size = .3
  ) +
  scale_fill_manual(values = reordered_colors) +
  geom_text(aes(label = round(Cindex.mean, 2)), hjust = -.5, size = 2) +
  labs(x = "", y = "C-index (95% CI)") +
  theme_light() +
  ylim(0.45, 1) +
  geom_hline(yintercept = 0.5,
             linetype = "dotted",
             color = "grey50") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    strip.text = element_text(color = "white"),
    strip.background = element_rect(
      color = "black",
      fill = "black",
      linetype = "solid"
    )
  )

## save
pdf("<path>/SuppFig_2.pdf",
    height = 6,
    width = 7)
b &
  theme(text = element_text(size = 10))
dev.off()

#############################################
####       Detailed unimodal (ST 4)      ####
#############################################
write_csv(cindex_individual_all, "results/IndividualCIndex.csv")

## ---> metrics for the text <--- ##
## median C-index General Health & Behaviour
round(median(cindex_individual_all[FeatureSet == "General health & lifestyle"]$Cindex.mean), 3)
round(range(cindex_individual_all[FeatureSet == "General health & lifestyle"]$Cindex.mean), 3)
# 0.713 (0.593 - 0.794)

## median C-index blood assays
round(median(cindex_individual_all[FeatureSet == "Clinical blood biomarkers"]$Cindex.mean), 3)
round(range(cindex_individual_all[FeatureSet == "Clinical blood biomarkers"]$Cindex.mean), 3)
# 0.695 (0.557 - 0.850)

## median PGS
round(median(cindex_individual_all[FeatureSet == "Polygenic scores"]$Cindex.mean), 3)
round(range(cindex_individual_all[FeatureSet == "Polygenic scores"]$Cindex.mean), 3)
# 0.564 (0.506 - 0.623)

#############################################
####    Top20 vs all-nonrandom (SFig 3)  ####
#############################################
cindex_individualnonrandom_list <- list()

## loop through each outcome
for (outc in outc_list) {
  ## load models
  feature.nonrandom.basic <-
    get(load(file = paste0("<path>", outc, ".RData")))
  feature.nonrandom.disdrug <-
    get(load(file = paste0("<path>", outc, ".RData")))
  feature.nonrandom.biomarker <-
    get(load(file = paste0("<path>", outc, ".RData")))
  feature.nonrandom.cardiopul <-
    get(load(file = paste0("<path>", outc, ".RData")))
  feature.nonrandom.bodycomp <-
    get(load(file = paste0("<path>", outc, ".RData")))
  feature.nonrandom.nmr <-
    get(load(file = paste0("<path>", outc, ".RData")))
  feature.nonrandom.pgs <-
    get(load(file = paste0("<path>", outc, ".RData")))
  
  ## create dataframe of performances
  cindex_data <- data.frame(
    Outcome = outc,
    FeatureSet = factor(
      c(
        "General health & lifestyle",
        "Diseases & drug intake",
        "Clinical blood biomarkers",
        "Cardiopulmonary parameters",
        "Body composition",
        "Plasma metabolites",
        "Polygenic scores"
      ),
      levels = c(
        "Polygenic scores",
        "Plasma metabolites",
        "Body composition",
        "Cardiopulmonary parameters",
        "Clinical blood biomarkers",
        "Diseases & drug intake",
        "General health & lifestyle"
      )
    ),
    Cindex.mean = c(
      feature.nonrandom.basic$Cindex.mean,
      feature.nonrandom.disdrug$Cindex.mean,
      feature.nonrandom.biomarker$Cindex.mean,
      feature.nonrandom.cardiopul$Cindex.mean,
      feature.nonrandom.bodycomp$Cindex.mean,
      feature.nonrandom.nmr$Cindex.mean,
      feature.nonrandom.pgs$Cindex.mean
    ),
    ci.low = c(
      feature.nonrandom.basic$ci.low,
      feature.nonrandom.disdrug$ci.low,
      feature.nonrandom.biomarker$ci.low,
      feature.nonrandom.cardiopul$ci.low,
      feature.nonrandom.bodycomp$ci.low,
      feature.nonrandom.nmr$ci.low,
      feature.nonrandom.pgs$ci.low
    ),
    ci.upp = c(
      feature.nonrandom.basic$ci.upp,
      feature.nonrandom.disdrug$ci.upp,
      feature.nonrandom.biomarker$ci.upp,
      feature.nonrandom.cardiopul$ci.upp,
      feature.nonrandom.bodycomp$ci.upp,
      feature.nonrandom.nmr$ci.upp,
      feature.nonrandom.pgs$ci.upp
    )
  )
  
  cindex_individualnonrandom_list[[outc]] <- cindex_data
  
}

## bind to dataset
cindex_individualnonrandom_all <-
  as.data.table(bind_rows(cindex_individualnonrandom_list))

## replace outcome labels
cindex_individualnonrandom_all$Outcome <-
  new_labels[cindex_individualnonrandom_all$Outcome]

## plot
a <-
  cindex_individual_all %>% left_join(
    cindex_individualnonrandom_all,
    join_by(Outcome == Outcome, FeatureSet == FeatureSet)
  ) %>%
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
      fill = as.factor(FeatureSet),
      xmin = ci.low.x,
      xmax = ci.upp.x,
      ymin = ci.low.y,
      ymax = ci.upp.y
    ),
    linewidth = 0.1,
    shape = 21,
    size = .4,
    alpha = .7
  ) +
  geom_pointrange(
    aes(
      fill = as.factor(FeatureSet),
      xmin = ci.low.x,
      xmax = ci.upp.x
    ),
    linewidth = 0.25,
    shape = 21,
    size = .5,
    alpha = .7,
    stroke = .1
  ) +
  labs(x = "C - Index (95% CI) 20 features", y = "C - Index (95% CI) all non-random features") +
  scale_fill_manual(values = reordered_colors) +
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
    legend.box.margin = margin(-10, -10, -10, -10),
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
  )

## save
pdf("<path>/SuppFig_3.pdf",
    height = 5 ,
    width = 6)
a & theme(text = element_text(size = 9))
dev.off()

### ---- metrics for the text ---- ###
## range difference non-random vs only top 20
range((
  cindex_individual_all %>% left_join(
    cindex_individualnonrandom_all,
    join_by(Outcome == Outcome, FeatureSet == FeatureSet)
  ) %>%
    filter(FeatureSet == "Diseases & drug intake") %>%
    mutate(delta = Cindex.mean.y - Cindex.mean.x)
)$delta)

## median
median((
  cindex_individual_all %>% left_join(
    cindex_individualnonrandom_all,
    join_by(Outcome == Outcome, FeatureSet == FeatureSet)
  ) %>%
    filter(FeatureSet == "Diseases & drug intake") %>%
    mutate(delta = Cindex.mean.y - Cindex.mean.x)
)$delta)
