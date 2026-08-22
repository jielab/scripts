#############################################################
####                     t-ROC analyses                  ####
#############################################################
# This script generates:
## Time dependent ROC curves for reviewer figure

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
require(survminer)
require(risksetROC)

## import data
ukb.phe <- fread("<path>")
ukb.phe.original <- ukb.phe
## import labels
lab.phe <- fread("<path>")

## define outcomes
outc_list <- c("<outcomes>")

## define follow-up
cdat_list <- c("<followup>")

############################################
####      Time-dependent ROC curves     ####
############################################

## data frames to store results
td_auc_list <- list()

# loop across outcomes
for (i in seq_along(outc_list)) {
  outc <- outc_list[i]
  cdat <- cdat_list[i]
  
  ## set seed
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
  # convert to data frame
  ukb.phe  <- as.data.frame(ukb.phe)
  
  ## load models
  loaded_object_name <-
    load(file = paste0("<path>", outc, ".RData"))
  feature.com <- get(loaded_object_name)
  
  ukb.test <- ukb.phe %>% filter(f.eid %in% v.idx)
  ukb.test <- as.data.table(ukb.test)
  
  ## join linear predictor to dataset
  ukb.test[, lin.pred := feature.com$linear.predictor]
  
  ## run function from risksetROC package
  td_auc <- risksetAUC(
    Stime = ukb.test[[cdat]],
    status = ukb.test[[outc]],
    marker = ukb.test$lin.pred,
    tmax = 10,
    plot = F
  )
  
  td_auc <- data.frame(Outcome = outc,
                       AUC = td_auc$AUC,
                       t = td_auc$utimes)
  
  td_auc_list[[outc]] <- td_auc
}

td_auc_data <- as.data.table(bind_rows(td_auc_list))

## rename outcomes
new_labels <- c("<nicenames>")

## rename outcomes
td_auc_data$Outcome <- new_labels[td_auc_data$Outcome]

## get first AUC of all
firsts <- td_auc_data %>%
  group_by(Outcome) %>%
  slice_min(t, n = 1, with_ties = FALSE) %>%
  transmute(Outcome, first_auc = AUC)

## get 5y and 10y auc for each outcome
targets <- tibble(target = c(2.5 , 5 , 7.5 , 10))
closest <- td_auc_data %>%
  right_join(targets, by = character()) %>%
  group_by(Outcome, target) %>%
  slice_min(abs(t - target), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(Outcome,
            target,
            t_closest = t,
            auc_target = AUC)

## get changes in AUC
deltas <- closest %>%
  left_join(firsts, by = "Outcome") %>%
  mutate(
    delta = auc_target - first_auc,
    label = scales::number(delta, accuracy = 0.001, signed = TRUE),
    x = target - 0.9,
    xarrow = target,
    y = pmin(auc_target + 0.1, 0.99)
  )

## plot
pdf("graphics/timedepROC.pdf",
    width = 7,
    height = 7)
td_auc_data %>%
  ggplot(aes(t, AUC)) +
  facet_wrap( ~ Outcome, ncol = 3, scales = "free") +
  geom_line(size = .2, color = "#0072B5FF") +
  labs(x = "Follow-up time", y = "AUC (t)") +
  theme_light() + ylim(0.4, 1) + xlim(0, 10) +
  geom_hline(yintercept = 0.5,
             linetype = "dotted",
             color = "grey50") +
  geom_segment(
    data = deltas,
    aes(
      x = xarrow,
      xend = xarrow,
      yend = auc_target + 0.02,
      y = y
    ),
    inherit.aes = FALSE,
    arrow = arrow(length = unit(0.03, "inches")),
    color = "grey40",
    linewidth = 0.3
  ) +
  geom_label(
    data = deltas,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 2,
    fill = NA,
    label.size = 0,
    label.r = unit(0.1, "lines")
  ) +
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
    text = element_text(size = 10)
  )
dev.off()
