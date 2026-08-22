#############################################################
####            Performance with less features           ####
#############################################################

# This script generates:
## Reviewer Fig 2: Performance of OBSCORE with less features

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
require(ggrepel)
require(stringr)

## define outcomes
outc_list <- c("<outcomes>")

## list to store results
cindex_data_list <- list()

## loop through each outcome
for (outc in outc_list) {
  
  ## load models
  m.20 <- get(load(file = paste0("<path>", outc, ".RData")))
  m.19 <- get(load(file = paste0("<path>", outc, ".20.RData")))
  m.18 <- get(load(file = paste0("<path>", outc, ".19.RData")))
  m.17 <- get(load(file = paste0("<path>", outc, ".18.RData")))
  m.16 <- get(load(file = paste0("<path>", outc, ".17.RData")))
  m.15 <- get(load(file = paste0("<path>", outc, ".16.RData")))
  m.14 <- get(load(file = paste0("<path>", outc, ".15.RData")))
  m.13 <- get(load(file = paste0("<path>", outc, ".14.RData")))
  m.12 <- get(load(file = paste0("<path>", outc, ".13.RData")))
  m.11 <- get(load(file = paste0("<path>", outc, ".12.RData")))
  m.10 <- get(load(file = paste0("<path>", outc, ".11.RData")))
  
  ## data frame with performances
  cindex_data <- data.frame(
    Outcome = outc,
    FeatureSet = factor(
      c(
        "20 Features",
        "19 Features",
        "18 Features",
        "17 Features",
        "16 Features",
        "15 Features",
        "14 Features",
        "13 Features",
        "12 Features",
        "11 Features",
        "10 Features"
      ),
      levels = c(
        "10 Features",
        "11 Features",
        "12 Features",
        "13 Features",
        "14 Features",
        "15 Features",
        "16 Features",
        "17 Features",
        "18 Features",
        "19 Features",
        "20 Features"
      )
    ),
    Cindex.mean = c(
      m.20$Cindex.mean,
      m.19$Cindex.mean,
      m.18$Cindex.mean,
      m.17$Cindex.mean,
      m.16$Cindex.mean,
      m.15$Cindex.mean,
      m.14$Cindex.mean,
      m.13$Cindex.mean,
      m.12$Cindex.mean,
      m.11$Cindex.mean,
      m.10$Cindex.mean
    ),
    ci.low = c(
      m.20$ci.low,
      m.19$ci.low,
      m.18$ci.low,
      m.17$ci.low,
      m.16$ci.low,
      m.15$ci.low,
      m.14$ci.low,
      m.13$ci.low,
      m.12$ci.low,
      m.11$ci.low,
      m.10$ci.low
    ),
    ci.upp = c(
      m.20$ci.upp,
      m.19$ci.upp,
      m.18$ci.upp,
      m.17$ci.upp,
      m.16$ci.upp,
      m.15$ci.upp,
      m.14$ci.upp,
      m.13$ci.upp,
      m.12$ci.upp,
      m.11$ci.upp,
      m.10$ci.upp
    )
  )
  
  ## append to the list
  cindex_data_list[[outc]] <- cindex_data
  
}

## combine all C-index data
cindex_data_all <- as.data.table(bind_rows(cindex_data_list))

## define outcome labels
new_labels <- c("<nicenames>")

## replace outcome names
cindex_data_all$Outcome <- new_labels[cindex_data_all$Outcome]

## colors for outcomes
cols <- c(
  "dodgerblue2",
  "#E31A1C",
  "green4",
  "#6A3D9A",
  "#FF7F00",
  "gold1",
  "skyblue2",
  "#FB9A99",
  "palegreen2",
  "#CAB2D6",
  "#FDBF6F",
  "gray70",
  "khaki2",
  "maroon",
  "orchid1",
  "deeppink1",
  "blue1",
  "steelblue4"
)
## levels
feat_levels <- paste(10:20, "Features")
## some refactoring
cindex_data_all2 <- cindex_data_all %>%
  mutate(
    nfeat = as.integer(str_extract(FeatureSet, "\\d+")),
    FeatureSet = factor(FeatureSet, levels = feat_levels)
  )
## baseline with 20
base20 <- cindex_data_all2 %>%
  filter(nfeat == 20) %>%
  select(Outcome, base = Cindex.mean)

med_labels <- cindex_data_all2 %>%
  filter(nfeat <= 20, nfeat >= 10) %>%
  left_join(base20, by = "Outcome") %>%
  mutate(diff_from20 = Cindex.mean - base) %>%
  filter(nfeat != 20) %>%
  group_by(FeatureSet, nfeat) %>%
  summarise(med_diff = median(diff_from20, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(label = paste0(round(med_diff, digits = 4)))

pdf("<path>/stepwiseexclude10features.pdf", width = 7, height = 2.75)
ggplot(cindex_data_all2,
       aes(
         x = FeatureSet,
         y = Cindex.mean,
         fill = Outcome,
         color = Outcome
       )) +
  coord_flip(clip = "off") +
  geom_line(linetype = "dashed",  aes(group = Outcome)) +
  geom_pointrange(
    aes(ymin = ci.low, ymax = ci.upp),
    linewidth = 0.25,
    shape = 21,
    size = 0.4,
    stroke = .1
  ) +
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  labs(x = "Exclusion of features", y = "C-index (95% CI)") +
  theme_light() +
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
    ),
    text = element_text(size = 7)
  ) +
  geom_text(
    data = med_labels,
    aes(x = FeatureSet, y = 0.9, label = label),
    inherit.aes = FALSE,
    size = 2,
    vjust = 1
  ) +
  geom_text(
    data = cindex_data_all2[FeatureSet == "20 Features"],
    aes(
      x = "20 Features",
      y = Cindex.mean,
      label = Outcome,
      color = Outcome
    ),
    inherit.aes = FALSE,
    angle = 90,
    size = 2,
    hjust = -0.2,
    vjust = 0.5
  )
dev.off()
