#############################################################
####       Cumulative incidences and predictor list      ####
#############################################################

# This script generates:
## Fig 1A:  A plot with cumulative disease incidences for 18 obesity outcomes stratified by BMI
## Fig 1B:  A plot displaying number of potential one hot encoded predictors per category
## ST 10 :  A predictor list table

rm(list = ls())
setwd("<path>")
options(stringsAsFactors = F)

## --> packages needed <-- ##
require(data.table)
require(readxl)
require(tidyverse)
require(magrittr)
require(ggrepel)
require(forcats)
require(survival)
require(survminer)
require(patchwork)
require(RColorBrewer)
require(ggsci)

## read data
ukb_cleaned <- fread("<path>")
## read labs
lab.ext <- fread("<path>")

##################################################
####    Cumulative incidence plot (Fig 1A)    ####
##################################################

## define outcomes
outcomes <- c("<outcomes>")

## define follow-up columns
followup <- c("<followup>")

## define custom outcome names
new_labels <- c("<nicenames>")

## custom labels match outcomes
names(new_labels) <- outcomes

## define BMI categories
ukb_cleaned <- ukb_cleaned %>%
  mutate(bmi_category = case_when(
    bmi >= 27 & bmi < 30 ~ "Overweight",
    bmi >= 30 ~ "Obese",
    TRUE ~ NA_character_
  ))

## list to store results
cuminc_bmi_list <- list()

## loop through BMI categories
for (bmi_cat in unique(ukb_cleaned$bmi_category)) {
  ## filter data for the current BMI category
  bmi_data <- ukb_cleaned %>%
    filter(bmi_category == bmi_cat)
  
  ## list to store results
  bmi_cuminc_list <- list()
  
  ## loop through outcomes
  for (i in seq_along(outcomes)) {
    outcome <- outcomes[i]
    followup_time <- followup[i]
    
    ## filter baseline and early cases
    subset_data <-
      bmi_data %>% filter(!(((
        eval(as.name(followup_time)) < 0.5
      )) &
        (eval(
          as.name(outcome)
        ) == 1)))
    
    ## fit Cox model
    surv_object <-
      Surv(subset_data[[followup_time]], subset_data[[outcome]])
    fit <- survfit(surv_object ~ 1, data = subset_data)
    
    ## get cumulative incidence at 10 years
    cuminc_10 <- fit$surv[fit$time <= 10]
    if (length(cuminc_10) > 0) {
      cuminc_value <- 1 - tail(cuminc_10, 1)
    } else {
      cuminc_value <- NA
    }
    
    ## append to the list
    bmi_cuminc_list[[outcome]] <- data.frame(
      outcome_label = new_labels[outcome],
      cuminc = cuminc_value,
      bmi_category = bmi_cat
    )
  }
  
  ## combine results for the current BMI category
  cuminc_bmi_list[[bmi_cat]] <- bind_rows(bmi_cuminc_list)
}

## combine results
cuminc_bmi_data <- bind_rows(cuminc_bmi_list)

## convert cumulative incidence to percentages
cuminc_bmi_data <- cuminc_bmi_data %>%
  mutate(cuminc_percent = cuminc * 100)

## create bar chart
a <-
  ggplot(cuminc_bmi_data,
         aes(
           x = cuminc_percent,
           y = reorder(outcome_label, cuminc_percent),
           fill = bmi_category
         )) +
  geom_bar(
    stat = "identity",
    position = "dodge",
    color = "black",
    size = .2,
    width = 1,
    alpha = 0.8
  ) +
  scale_x_continuous(
    limits = c(0, 17),
    breaks = seq(0, 17, by = 5),
    expand = expansion(mult = c(0.01, 0.07))
  ) +
  labs(x = "Cumulative 10-year incidence (%)", y = "", ) +
  scale_color_jama(guide = guide_legend(
    ncol = 1,
    byrow = TRUE,
    keyheight = unit(.01, "cm"),
    keywidth = unit(.01, "cm")
  )) +
  scale_fill_jama() +
  theme_light() +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  )

##################################################
####    Number of predictors plot (Fig 1B)    ####
##################################################
## renaming
lab.ext <- lab.ext %>%
  mutate(
    category = case_when(
      category == "Ancestry" ~ "Ancestries",
      category == "Basic demographics" ~ "Basic demographics",
      category == "Health" ~ "Health factors",
      category == "Family history" ~ "Family history of diseases",
      category == "Socioeconomic" ~ "Socioeconomic factors",
      category == "Diet" ~ "Diet preferences",
      category == "Diseases" ~ "Diseases",
      category == "Drugs" ~ "Drugs",
      category == "Biomarker" ~ "Clinical biomarkers",
      category == "Blood cell counts" ~ "Blood cell counts",
      category == "Cardiovascular" ~ "Cardiovascular parameters",
      category == "Pulmonary" ~ "Pulmonary parameters",
      category == "Body composition" ~ "Body composition parameters",
      category == "NMR" ~ "Plasma metabolites",
      category == "PGS" ~ "Polygenic scores",
      TRUE ~ category
    )
  )

## ordering
category_order <-
  rev(
    c(
      "Ancestries",
      "Basic demographics",
      "Health factors",
      "Family history of diseases",
      "Socioeconomic factors",
      "Diet preferences",
      "Diseases",
      "Drugs",
      "Clinical biomarkers",
      "Blood cell counts",
      "Cardiovascular parameters",
      "Pulmonary parameters",
      "Body composition parameters",
      "Plasma metabolites",
      "Polygenic scores"
    )
  )

## define order for plotting
lab.ext$category <-
  factor(lab.ext$category, levels = category_order)

## remove recruitment centre
category_counts <- lab.ext %>% filter(!is.na(category)) %>%
  filter(short_name != "centre") %>% count(category) %>% na.omit()

## color palette for labels
cat.col <- data.table(category = unique(c(lab.ext$category)),
                      cl = colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(16))
cat.col.vector <- setNames(cat.col$cl, cat.col$category)

## function to transform data (to break axis)
trans <- function(x) {
  pmin(x, 300) + 0.05 * pmax(x - 300, 0)
}

## define y-axis ticks
yticks <- c(0, 100, 200, 300, 1000)

## transform data onto the display scale
category_counts$n_t <- trans(category_counts$n)

## create bar chart
b <-
  ggplot(category_counts, aes(x = category, y = n_t, fill = category)) +
  geom_bar(stat = "identity",
           color = "black",
           size = 0.2) +
  coord_flip() +
  scale_fill_manual(values = cat.col.vector) +
  geom_label(
    aes(label = n, y = 0),
    hjust = 0,
    size = 1.7,
    alpha = 0.8,
    fill = "white",
    color = "black",
    label.size = NA
  ) +
  labs(x = "", y = "Number of predictors") +
  theme_bw() +
  scale_y_continuous(limits = c(0, NA),
                     breaks = trans(yticks),
                     labels = yticks) +
  geom_rect(aes(
    xmin = 0,
    xmax = Inf,
    ymin = trans(320),
    ymax = trans(680)
  ), fill = "white") +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  )

layout <- "
AB
"

jpeg(
  "<path>/Fig_1.jpg",
  width = 7.5,
  height = 3,
  units = "in",
  res = 1500
)
a + b + plot_layout(design = layout) &
  (theme(text = element_text(size = 7)))
dev.off()

##################################################
####      Predictor list (Supp Table 10)      ####
##################################################
tmp <- lab.ext %>% filter(!is.na(category)) %>%
  filter(short_name != "centre") %>% select(label_new, category)
write_csv(tmp, "<path>/Predictors.csv")
