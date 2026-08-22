#############################################################
####               Multimodal performances               ####
#############################################################

# This script generates:
## Fig 3: Performances of multimodal models
## ST 5 : Performances of multimodal models table

rm(list = ls())
setwd("<path>")

## --> packages needed <-- ##
require(data.table)
require(tidyverse)
require(dplyr)
require(magrittr)
require(ggrepel)
require(patchwork)
require(tidytext)
require(ggsci)

## define outcomes
outc_list <- c("<outcomes>")

## define outcome labels
new_labels <- c("<nicenames>")

#############################################
####   Multimodal performances (Fig 3)   ####
#############################################
## list to store results
cindex_data_list <- list()

# loop through each outcome
for (outc in outc_list) {
  # load models
  feature.top.basic <-
    get(load(
      file = paste0(
        "output/overweight.model.feature.top.Basic parameters.",
        outc,
        ".RData"
      )
    ))
  feature.top.disdrug <-
    get(load(
      file = paste0(
        "output/overweight.model.feature.top.+Diseases & Drugs.",
        outc,
        ".RData"
      )
    ))
  feature.top.biomarker <-
    get(load(
      file = paste0(
        "output/overweight.model.feature.top.+Biomarker.",
        outc,
        ".RData"
      )
    ))
  feature.top.cardiopul <-
    get(load(
      file = paste0(
        "output/overweight.model.feature.top.+Cardiopulmonary.",
        outc,
        ".RData"
      )
    ))
  feature.top.bodycomp <-
    get(load(
      file = paste0(
        "output/overweight.model.feature.top.+Body composition.",
        outc,
        ".RData"
      )
    ))
  feature.top.nmr <-
    get(load(
      file = paste0(
        "output/overweight.model.feature.top.+NMR.",
        outc,
        ".RData"
      )
    ))
  feature.top.pgs <-
    get(load(
      file = paste0(
        "output/overweight.model.feature.top.+PGS.",
        outc,
        ".RData"
      )
    ))
  
  ## create dataframe of performances
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
  
}

## combine all C-index data
cindex_data_all <- as.data.table(bind_rows(cindex_data_list))

## replace outcome names
cindex_data_all$Outcome <- new_labels[cindex_data_all$Outcome]

## define plot order based on decreasing mean C-index
ordered_outcomes <-
  cindex_data_all[FeatureSet == "+Polygenic scores"][order(-Cindex.mean), Outcome]

## reorder order of outcomes
cindex_data_all$Outcome <-
  factor(cindex_data_all$Outcome, levels = unique(ordered_outcomes))

## prepare ribbon shading dataset
cindex_data_ribbon <- cindex_data_all %>%
  group_by(Outcome) %>%
  mutate(ymin = min(Cindex.mean), ymax = max(Cindex.mean))

## prepare order of data modalities
cindex_data_all$FeatureSet <- factor(cindex_data_all$FeatureSet,
                                     levels = rev(levels(cindex_data_all$FeatureSet)))

## plot and save
pdf("<path>/Fig_3.pdf", height = 7, width = 7)
ggplot(cindex_data_all, aes(x = FeatureSet, y = Cindex.mean)) +
  coord_flip() +
  facet_wrap( ~ Outcome, ncol = 3) +
  geom_rect(
    data = cindex_data_ribbon,
    aes(
      xmin = -Inf,
      xmax = Inf,
      ymin = 0.5,
      ymax = ymax,
      group = Outcome
    ),
    fill = "#4DBBD5B2",
    alpha = .025
  ) +
  geom_pointrange(
    aes(ymin = ci.low, ymax = ci.upp, fill = FeatureSet),
    linewidth = 0.25,
    shape = 21,
    size = 0.35
  ) +
  scale_fill_nejm() +
  geom_text(aes(label = round(Cindex.mean, 2)), hjust = -0.5, size = 2) +
  labs(x = "", y = "C-index (95% CI)") +
  theme_light() +
  ylim(0.5, 1) +
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
dev.off()

#############################################
####    Multimodal performances (ST5)    ####
#############################################
write_csv(cindex_data_all, "<path>/StepwiseCIndex.csv")

### ---- metrics for the text ---- ###
## median delta +Clinical blood biomarkers
round(median((
  cindex_data_all[FeatureSet == "General health & lifestyle"] %>%
    left_join(cindex_data_all[FeatureSet == "+Clinical blood biomarkers"], by =
                "Outcome") %>%
    mutate(delta = Cindex.mean.y - Cindex.mean.x)
)$delta), 3)
# 0.025

## median delta +Clinical blood biomarkers to +Cardiopulmonary parameters
round(median((
  cindex_data_all[FeatureSet == "+Clinical blood biomarkers"] %>%
    left_join(cindex_data_all[FeatureSet == "+Cardiopulmonary parameters"], by =
                "Outcome") %>%
    mutate(delta = Cindex.mean.y - Cindex.mean.x)
)$delta), 4)
# -0.0006

## median delta +Cardiopulmonary parameters to +Body composition
round(median((
  cindex_data_all[FeatureSet == "+Cardiopulmonary parameters"] %>%
    left_join(cindex_data_all[FeatureSet == "+Body composition"], by =
                "Outcome") %>%
    mutate(delta = Cindex.mean.y - Cindex.mean.x)
)$delta), 3)
# 0.002

## median delta +Body composition to +NMR
round(median((
  cindex_data_all[FeatureSet == "+Body composition"] %>%
    left_join(cindex_data_all[FeatureSet == "+Plasma metabolites"], by =
                "Outcome") %>%
    mutate(delta = Cindex.mean.y - Cindex.mean.x)
)$delta), 4)
# -0.0003

## median delta +NMR to +PGS
round(median((
  cindex_data_all[FeatureSet == "+Plasma metabolites"] %>%
    left_join(cindex_data_all[FeatureSet == "+Polygenic scores"], by =
                "Outcome") %>%
    mutate(delta = Cindex.mean.y - Cindex.mean.x)
)$delta), 3)
# 0.003

## median (range) all modalities
round(median((cindex_data_all[FeatureSet == "+Polygenic scores"])$Cindex.mean), 3)
round(range((cindex_data_all[FeatureSet == "+Polygenic scores"])$Cindex.mean), 3)
# 0.75 (0.634, 0.864)
