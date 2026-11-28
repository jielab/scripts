source(file.path(Sys.getenv("FLT_CODE"), "runtime.R"))
library(dplyr)
library(ggplot2)
library(patchwork)

project <- Sys.getenv("FLT_OUT")
input <- Sys.getenv("FLT_COHORT")
outdir <- file.path(project, "01lm")
figdir <- file.path(Sys.getenv("FLT_PAPER"), "01lm/figures")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

# 只读取分析结果；三个热图共用色阶和图例。
r <- read.csv(file.path(outdir, "lm_results.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(r) > 0, !anyDuplicated(r[c("trait", "variable")]))
liking <- unique(r[c("variable", "food_liking")])
limit <- 0.20  # 仅限制显示范围，原始 beta 不变。
r <- r %>%
  mutate(category = ifelse(grepl("^Brain", category), "Brain behavioral phenotypes", category),
         variable = factor(variable, levels = liking$variable),
         valid = status == "ok" & is.finite(beta) & is.finite(se) & se > 0 & is.finite(FDR),
         display_beta = pmax(-limit, pmin(limit, beta)))
stopifnot(any(r$valid))
labels <- c(fat_pct="Body fat percentage", whole_body_fat_mass="Whole body fat mass", bmi="BMI",
 smoke_status="Smoking status", alcohol_status="Alcohol status", alcohol_freq="Alcohol frequency",
 sleep_duration="Sleep duration", phone_call_time="Phone call time", tv_time="TV time",
 computer_time="Computer time", pa_met_tot="Total physical activity (MET)", summed_mvpa="Moderate/vigorous activity",
 walking_duration="Walking duration", walk_time="Walking time", liking_gym="Liking gym",
 liking_stairs="Liking stairs", liking_sweat="Liking sweat", education_degree="Education degree",
 tdi="Townsend deprivation index", inc="Income", fluid_intelligence.z="Fluid intelligence",
 reaction_time.z="Reaction time", pairs_matching.z="Pairs matching", numeric_memory.z="Numeric memory",
 life_satisfaction.z="Life satisfaction", subjective_wellbeing.z="Subjective wellbeing",
 phq4.z="PHQ-4", phq9.z="PHQ-9", cidi.z="CIDI", gad7.z="GAD-7",
 prospective_memory.b="Prospective memory", general_mental_health.b="General mental health")
for (part in c("leg", "arm")) for (side in c("right", "left")) {
 for (measure in c("fat_pct", "fat_mass", "fat_free_mass", "predicted_mass")) {
  key <- paste(part, measure, side, sep="_")
  desc <- c(fat_pct="fat percentage", fat_mass="fat mass", fat_free_mass="fat-free mass", predicted_mass="predicted mass")
  labels[key] <- paste0(tools::toTitleCase(part), " ", desc[measure], " (", side, ")")
 }
}

stopifnot(all(r$trait %in% names(labels)))
panel <- function(group) {
  z <- r %>% filter(category == group) %>%
    mutate(trait = factor(trait, levels = rev(unique(trait))))
  ggplot(z, aes(variable, trait)) +
    geom_tile(fill = "white", colour = "grey80") +
    geom_point(data = filter(z, valid), aes(fill = display_beta, size = abs(display_beta)), shape = 22, colour = NA) +
    geom_point(data = filter(z, valid & FDR >= 0.05), shape = 4, size = 2) +
    geom_text(data = filter(z, !valid), label = "NA", size = 3, colour = "grey40") +
    scale_x_discrete(labels = setNames(liking$food_liking, liking$variable), drop = FALSE) +
    scale_y_discrete(labels = labels, drop = FALSE) +
    scale_fill_gradient2(low = "#29365F", mid = "white", high = "#852121",
                         limits = c(-limit, limit), name = "Beta (capped)") +
    scale_size_area(max_size = 10, limits = c(0, limit), guide = "none") +
    labs(title = group, x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(colour = "black"), plot.title = element_text(face = "bold", size = 12))
}
p <- wrap_plots(lapply(unique(r$category), panel), nrow = 1, guides = "collect") +
  plot_annotation(title = "Food liking and other phenotypes",
    caption = paste("Beta per 1-point higher FLT; outcome units differ. Colour and square area saturate at |beta| =", limit,
                    "\nx: FDR >= 0.05; NA: unavailable estimate. BH-FDR within each phenotype (3 FLTs)."))
p <- p & theme(legend.position = "bottom")
ggsave(file.path(figdir, "otherphe_lm_heatmap_beta.pdf"), p, width = 18, height = 10, bg = "white")
