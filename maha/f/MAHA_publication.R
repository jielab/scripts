dir0 <- ifelse(Sys.info()[["sysname"]] == "Windows", "D:", "/mnt/d")
outdir <- Sys.getenv("MAHA_OUTDIR", unset = file.path(dir0, "analysis", "maha"))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
.this_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.this_file <- if (length(.this_arg)) sub("^--file=", "", .this_arg[[1]]) else file.path(dir0, "scripts", "maha", "f", "MAHA_publication.R")
.this_dir <- dirname(normalizePath(.this_file, winslash = "/", mustWork = FALSE))

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(writexl)
  library(patchwork)
  library(cowplot)
  library(scales)
})
source(file.path(.this_dir, "comm.f.R"))
source(file.path(dir0, "scripts", "0f", "plot.f.R"))

args <- commandArgs(trailingOnly = TRUE)
cohort_arg <- grep("^--cohort=", args, value = TRUE)
scope <- if (length(cohort_arg)) sub("^--cohort=", "", cohort_arg[[length(cohort_arg)]]) else "all"
if (!scope %in% c("ukb", "nhanes", "chns", "all")) stop("--cohort must be ukb, nhanes, chns, or all", call. = FALSE)

pal_diet <- c(
  "MAHA" = "#F26D60", "DASH" = "#16B9C0", "MIND" = "#B97AF7", "MEDI" = "#7CAE00",
  "MAHA-balanced" = "#F4A261", "MAHA-strict" = "#A873E8",
  "MAHA-no dairy" = "#2A9D8F", "MAHA-no protein" = "#577590",
  "MAHA-density" = "#F4A261", "MAHA-residual" = "#2A9D8F"
)
pal_strata <- c("low" = "#66C2A5", "middle" = "#FC8D62", "high" = "#8DA0CB")
diet.lst <- c(
  dash = "DASH", medi = "MEDI", mind = "MIND", maha = "MAHA",
  balanced = "-balanced", strict = "-strict",
  no_dairy = "-no_dairy", no_protein = "-no_protein"
)
primary_diets <- c("DASH", "MEDI", "MIND", "MAHA")
maha_variants <- c("MAHA", "MAHA-balanced", "MAHA-strict", "MAHA-no dairy", "MAHA-no protein")
diet_display <- c(
  "DASH" = "DASH", "MEDI" = "MEDI", "MIND" = "MIND", "MAHA" = "MAHA",
  "MAHA-balanced" = "-balanced", "MAHA-strict" = "-strict",
  "MAHA-no dairy" = "-no_dairy", "MAHA-no protein" = "-no_protein"
)
diet_plot_order <- c("DASH", "MEDI", "MIND", "MAHA", "MAHA-balanced", "MAHA-strict", "MAHA-no dairy", "MAHA-no protein")

ukb_outcome_plot_levels <- c(
  "Coronary artery disease", "Ischemic stroke", "Heart failure",
  "Type 2 diabetes", "Chronic kidney disease", "All-cause mortality"
)

pub_theme <- function(base_size = 13) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, size = base_size * 1.05),
      plot.subtitle = element_text(size = base_size * .82, color = "#444444"),
      axis.title = element_text(face = "bold"),
      axis.title.y = element_text(margin = margin(r = 12)),
      axis.text = element_text(color = "black", face = "bold"),
      legend.text = element_text(face = "bold"),
      legend.position = "bottom",
      legend.title = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      panel.grid = element_blank(),
      plot.margin = margin(7, 8, 7, 8)
    )
}

aux_path <- function(cohort, filename) file.path(outdir, cohort, filename)
pub_path <- function(filename) file.path(outdir, filename)

wb_sheets <- function(path) {
  if (!file.exists(path)) stop("Required intermediate workbook does not exist: ", path, call. = FALSE)
  setNames(lapply(excel_sheets(path), function(sh) read_excel(path, sheet = sh)), excel_sheets(path))
}
get_sheet <- function(wb, name, required = TRUE) {
  if (name %in% names(wb)) return(as_tibble(wb[[name]]))
  if (required) stop("Required sheet missing: ", name, call. = FALSE)
  NULL
}
write_book <- function(x, filename) writexl::write_xlsx(x, pub_path(filename))
save_pub <- function(plot, filename, width, height, dpi = 500) {
  ggsave(pub_path(filename), plot = plot, width = width, height = height, dpi = dpi, bg = "white")
  message("[publication] wrote ", filename)
}
copy_if_exists <- function(from, to) {
  if (file.exists(from)) {
    file.copy(from, to, overwrite = TRUE)
    TRUE
  } else FALSE
}

supplement_renumbering <- tibble::tribble(
  ~old_number, ~new_number, ~content, ~tmp_stem, ~final_stem, ~obsolete_stem,
  10L,  6L, "CHNS fixed-wave mortality sensitivity",
  "FigS6_tmp.chns.mortality_sensitivity", "FigS6.chns.mortality_sensitivity", "FigS10.chns.mortality_sensitivity",
   9L,  7L, "CHNS high-event mortality sensitivity",
  "FigS7_tmp.chns.high_event_sensitivity", "FigS7.chns.high_event_sensitivity", "FigS9.chns.high_event_sensitivity",
   6L,  8L, "UKB construct stress tests",
  "FigS8_tmp.ukb.construct_stress_tests", "FigS8.ukb.construct_stress_tests", "FigS6.ukb.construct_stress_tests",
   7L,  9L, "NHANES criterion and construct validation",
  "FigS9_tmp.nhanes.validation", "FigS9.nhanes.validation", "FigS7.nhanes.validation",
   8L, 10L, "NHANES construct profiles and score concordance",
  "FigS10_tmp.nhanes.construct_concordance", "FigS10.nhanes.construct_concordance", "FigS8.nhanes.construct_concordance"
)

finalize_supplement_renumbering <- function(cohort = "all") {
  selected <- if (cohort == "all") supplement_renumbering else
    supplement_renumbering %>% filter(grepl(paste0(".",cohort,"."),tmp_stem,fixed=TRUE))
  jobs <- tidyr::crossing(selected, suffix = c(".png", ".out.xlsx")) %>%
    mutate(from = pub_path(paste0(tmp_stem, suffix)), to = pub_path(paste0(final_stem, suffix)))
  info <- file.info(jobs$from)
  missing <- jobs$from[is.na(info$size) | info$size < 1000]
  if (length(missing)) {
    stop(
      "Supplement renumbering was not finalized because temporary outputs are missing or incomplete: ",
      paste(basename(missing), collapse = ", "), call. = FALSE
    )
  }

  copied <- Map(function(from, to) file.copy(from, to, overwrite = TRUE), jobs$from, jobs$to)
  if (!all(unlist(copied, use.names = FALSE))) {
    stop("Could not copy every verified temporary supplement output to its final filename.", call. = FALSE)
  }
  source_md5 <- unname(tools::md5sum(jobs$from))
  final_md5 <- unname(tools::md5sum(jobs$to))
  if (!identical(source_md5, final_md5)) {
    stop("Supplement renumbering checksum verification failed; temporary files were retained.", call. = FALSE)
  }

  obsolete <- unlist(lapply(selected$obsolete_stem, function(stem) {
    pub_path(paste0(stem, c(".png", ".out.xlsx")))
  }), use.names = FALSE)
  obsolete <- setdiff(obsolete, jobs$to)
  unlink(obsolete[file.exists(obsolete)], force = TRUE)
  unlink(jobs$from, force = TRUE)
  if (any(file.exists(jobs$from))) stop("Temporary supplement outputs could not be removed after verification.", call. = FALSE)

  readr::write_tsv(
    supplement_renumbering %>% select(old_number, new_number, content, final_stem),
    pub_path("FigS_renumbering_manifest.tsv")
  )
  message("[publication] finalized verified supplement filenames for scope=",cohort)
}

align_publication_grid <- function(rows, rel_widths = NULL, rel_heights = NULL,
                                   row_gap = .10, side_pad = .02,
                                   axis_text_size = 10.5, axis_title_size = 11.5) {
  if (!is.list(rows) || !length(rows) || !all(vapply(rows, is.list, logical(1)))) {
    stop("rows must be a non-empty list of plot lists.", call. = FALSE)
  }
  ncols <- vapply(rows, length, integer(1))
  if (length(unique(ncols)) != 1L) stop("All publication grid rows must have the same number of panels.", call. = FALSE)
  widths_by_row <- if (is.null(rel_widths)) {
    rep(list(NULL), length(rows))
  } else if (is.list(rel_widths)) {
    if (length(rel_widths) != length(rows)) stop("rel_widths list must match rows.", call. = FALSE)
    rel_widths
  } else {
    rep(list(rel_widths), length(rows))
  }
  themed <- lapply(unlist(rows, recursive = FALSE), function(p) {
    if (inherits(p, "ggplot")) {
      p + theme(
        axis.title = element_text(face = "bold", size = axis_title_size),
        axis.text = element_text(face = "bold", size = axis_text_size, color = "black")
      )
    } else p
  })
  # Align left/right boundaries within each column across rows. Aligning plot
  # heights globally can push top annotations or titles outside the canvas when
  # another panel has a tall legend, so top/bottom alignment is applied only
  # while the panels are rebuilt within each row.
  row_count <- length(rows)
  col_count <- ncols[[1]]
  aligned_by_column <- lapply(seq_len(col_count), function(column) {
    indices <- column + (seq_len(row_count) - 1L) * col_count
    cowplot::align_plots(plotlist = themed[indices], align = "v", axis = "lr")
  })
  row_grobs <- Map(function(row_index, widths) {
    row_panels <- lapply(aligned_by_column, `[[`, row_index)
    args <- list(plotlist = row_panels, nrow = 1, align = "h", axis = "tb")
    if (!is.null(widths)) args$rel_widths <- widths
    row <- do.call(cowplot::plot_grid, args)
    if (is.finite(side_pad) && side_pad > 0 && side_pad < .25) {
      row <- cowplot::plot_grid(NULL, row, NULL, nrow = 1,
                               rel_widths = c(side_pad, 1 - 2 * side_pad, side_pad))
    }
    row
  }, seq_len(row_count), widths_by_row)

  base_heights <- if (is.null(rel_heights)) rep(1, length(rows)) else rel_heights
  if (length(base_heights) != length(rows)) stop("rel_heights must match rows.", call. = FALSE)
  if (length(row_grobs) > 1L && is.finite(row_gap) && row_gap > 0) {
    out <- vector("list", length(row_grobs) * 2L - 1L)
    out[seq(1L, length(out), 2L)] <- row_grobs
    out[seq(2L, length(out), 2L)] <- rep(list(NULL), length(row_grobs) - 1L)
    gap_height <- row_gap * mean(base_heights)
    out_heights <- as.numeric(rbind(head(base_heights, -1L), rep(gap_height, length(base_heights) - 1L)))
    row_grobs <- out
    base_heights <- c(out_heights, tail(base_heights, 1L))
  }
  cowplot::plot_grid(plotlist = row_grobs, ncol = 1, align = "none", rel_heights = base_heights)
}


assert_contiguous_figure_numbers <- function(prefix, expected) {
  escaped_prefix <- stringr::str_replace_all(prefix, "([\\W])", "\\\\\\1")
  pattern <- paste0("^", escaped_prefix, "([0-9]+)\\..+\\.png$")
  files <- list.files(outdir, pattern = pattern, full.names = FALSE)
  numbers <- sort(as.integer(stringr::str_match(files, pattern)[, 2]))
  if (!identical(numbers, as.integer(expected))) {
    stop(
      prefix, " figure numbering is not contiguous. Expected ",
      paste(expected, collapse = ", "), "; found ", paste(numbers, collapse = ", "),
      ". Files: ", paste(sort(files), collapse = ", "), call. = FALSE
    )
  }
  message("[publication] verified contiguous ", prefix, " numbering: ", min(expected), "-", max(expected))
}

parse_ci <- function(x) {
  m <- stringr::str_match(as.character(x), "([0-9.]+)\\s*\\(([0-9.]+),\\s*([0-9.]+)\\)")
  tibble(est = suppressWarnings(as.numeric(m[,2])), lo = suppressWarnings(as.numeric(m[,3])), hi = suppressWarnings(as.numeric(m[,4])))
}

normalize_chns_label <- function(x) {
  recode(as.character(x),
    "MAHA-CN foodcode-enhanced" = "MAHA",
    "MAHA-CN expanded proxy" = "MAHA-balanced",
    "MAHA-CN core" = "MAHA-strict",
    "DASH-CN proxy" = "DASH", "DASH-CN foodcode-complete" = "DASH",
    "MIND-CN proxy" = "MIND", "MIND-CN foodcode-complete" = "MIND",
    "MEDI-CN proxy" = "MEDI", "MEDI-CN foodcode-complete" = "MEDI",
    .default = as.character(x)
  )
}

forest_multi <- function(d, title, xlab, ref = 1, estimate = "HR", low = "LCI", high = "UCI",
                         outcome = "Outcome", diet = "Diet", diets = primary_diets,
                         hollow = character(), limits = NULL, show_text = FALSE) {
  d <- d %>% filter(.data[[diet]] %in% diets, is.finite(.data[[estimate]]), is.finite(.data[[low]]), is.finite(.data[[high]]))
  observed_outcomes <- unique(as.character(d[[outcome]]))
  if (length(observed_outcomes) > 0 && all(observed_outcomes %in% ukb_outcome_plot_levels)) {
    out_levels <- intersect(ukb_outcome_plot_levels, observed_outcomes)
  } else {
    out_levels <- rev(observed_outcomes)
  }
  d <- d %>% mutate(
    .out = factor(as.character(.data[[outcome]]), levels = out_levels),
    .diet = factor(as.character(.data[[diet]]), levels = diets),
    .y0 = as.numeric(.out),
    .off = scales::rescale(as.numeric(.diet), to = c(-.28, .28)),
    .y = .y0 + .off,
    .secondary = as.character(.diet) %in% hollow
  )
  p <- ggplot(d) +
    geom_vline(xintercept = ref, linetype = 2, color = "#777777", linewidth = .45) +
    geom_segment(aes(x = .data[[low]], xend = .data[[high]], y = .y, yend = .y, color = .diet), linewidth = .75) +
    geom_point(data = d %>% filter(!.secondary), aes(x = .data[[estimate]], y = .y, color = .diet), size = 2.8) +
    geom_point(data = d %>% filter(.secondary), aes(x = .data[[estimate]], y = .y, color = .diet), shape = 21, fill = "white", stroke = .9, size = 3.0) +
    scale_color_manual(values = pal_diet, breaks = diets, labels = unname(diet_display[diets]), drop = FALSE) +
    scale_y_continuous(breaks = seq_along(out_levels), labels = out_levels, expand = expansion(add = .65)) +
    labs(title = title, x = xlab, y = NULL) + pub_theme() +
    guides(color = guide_legend(
      nrow = if (length(diets) > 4) 2 else 1, byrow = TRUE,
      override.aes = list(
      shape = ifelse(diets %in% hollow, 21, 16),
      fill = ifelse(diets %in% hollow, "white", unname(pal_diet[diets]))
    )))
  if (!is.null(limits)) p <- p + coord_cartesian(xlim = limits)
  if (show_text) {
    d <- d %>% mutate(.lab = sprintf("%.2f (%.2f, %.2f)", .data[[estimate]], .data[[low]], .data[[high]]))
    p <- p + geom_text(data = d, aes(x = .data[[high]], y = .y, label = .lab, color = .diet), hjust = -.08, size = 2.7, show.legend = FALSE)
  }
  p
}

forest_single <- function(d, title, xlab, ref = 1, estimate = "HR", low = "CI_low", high = "CI_high",
                          label = "Diet", color = "Diet", hollow = character(), limits = NULL,
                          right_label = NULL, zone = NULL) {
  diet_order_full <- c(primary_diets, setdiff(maha_variants, "MAHA"))
  d <- d %>% filter(is.finite(.data[[estimate]]), is.finite(.data[[low]]), is.finite(.data[[high]])) %>%
    mutate(.lab_key = as.character(.data[[label]])) %>%
    arrange(match(.lab_key, diet_order_full)) %>%
    mutate(.lab = dplyr::coalesce(unname(diet_display[.lab_key]), .lab_key), .col = as.character(.data[[color]]), .y = rev(seq_len(n())))
  col_levels <- unique(d$.col)
  legend_order <- intersect(diet_order_full, col_levels)
  col_values <- unname(pal_diet[col_levels]); names(col_values) <- col_levels
  missing_col <- is.na(col_values)
  if (any(missing_col)) col_values[missing_col] <- scales::hue_pal(l = 55, c = 90)(sum(missing_col))
  p <- ggplot(d)
  if (!is.null(zone) && length(zone) == 2) p <- p + annotate("rect", xmin = zone[1], xmax = zone[2], ymin = -Inf, ymax = Inf, fill = "grey70", alpha = .14)
  p <- p +
    geom_vline(xintercept = ref, linetype = 2, color = "#777777", linewidth = .45) +
    geom_segment(aes(x = .data[[low]], xend = .data[[high]], y = .y, yend = .y, color = .col), linewidth = .8) +
    geom_point(data = d %>% filter(!.lab_key %in% hollow), aes(x = .data[[estimate]], y = .y, color = .col), size = 2.9) +
    geom_point(data = d %>% filter(.lab_key %in% hollow), aes(x = .data[[estimate]], y = .y, color = .col), shape = 21, fill = "white", stroke = .9, size = 3.1) +
    scale_color_manual(values = col_values, breaks = legend_order,
                       labels = unname(dplyr::coalesce(diet_display[legend_order], legend_order)),
                       drop = FALSE) +
    scale_y_continuous(breaks = d$.y, labels = d$.lab, expand = expansion(add = .78)) +
    labs(title = title, x = xlab, y = NULL) + pub_theme()
  if (!is.null(right_label) && right_label %in% names(d)) {
    xr <- range(c(d[[low]], d[[high]]), finite = TRUE); span <- diff(xr)
    p <- p + geom_text(aes(x = .data[[high]] + .04 * span, y = .y, label = .data[[right_label]], color = .col), hjust = 0, size = 3.15, fontface = "bold", show.legend = FALSE)
  }
  if (!is.null(limits)) {
    p <- p + coord_cartesian(xlim = limits, clip = "off")
  } else if (!is.null(right_label)) {
    p <- p + scale_x_continuous(expand = expansion(mult = c(.04, .62))) +
      theme(plot.margin = margin(7, 16, 7, 8))
  }
  p + guides(color = guide_legend(
    nrow = if (length(legend_order) > 4) 2 else 1, byrow = TRUE,
    override.aes = list(
      shape = ifelse(legend_order %in% hollow, 21, 16),
      fill = ifelse(legend_order %in% hollow, "white", unname(pal_diet[legend_order]))
    )
  ))
}

risk_bar <- function(d, title, xlab) {
  if ("Facet_group" %in% names(d)) {
    facet_col <- "Facet_group"
    line_col <- "Line_group"
  } else {
    wanted <- if (grepl("MAHA", xlab, ignore.case = TRUE)) "diet.maha.3c" else "diet.dash.3c"
    other <- if (wanted == "diet.maha.3c") "diet.dash.3c" else "diet.maha.3c"
    facet_col <- if (wanted %in% names(d)) wanted else intersect(c("diet.dash.3c", "diet.maha.3c"), names(d))[1]
    line_col <- if (other %in% names(d)) other else setdiff(intersect(c("diet.maha.3c", "diet.dash.3c"), names(d)), facet_col)[1]
  }
  if (is.na(facet_col) || is.na(line_col)) stop("Cannot identify risk-strata columns.", call. = FALSE)
  ycol <- if ("risk_pct" %in% names(d)) "risk_pct" else "mean"
  lcol <- if ("risk_lower_pct" %in% names(d)) "risk_lower_pct" else "lower_plot"
  ucol <- if ("risk_upper_pct" %in% names(d)) "risk_upper_pct" else "upper_plot"
  d <- d %>% transmute(
    Facet = factor(tolower(as.character(.data[[facet_col]])), levels = c("low", "middle", "high")),
    Line = factor(tolower(as.character(.data[[line_col]])), levels = c("low", "middle", "high")),
    Risk = as.numeric(.data[[ycol]]), Low = as.numeric(.data[[lcol]]), High = as.numeric(.data[[ucol]])
  ) %>% drop_na(Facet, Line, Risk)
  ggplot(d, aes(Facet, Risk, fill = Line)) +
    geom_col(position = position_dodge(width = .78), width = .68, alpha = .92) +
    geom_errorbar(aes(ymin = Low, ymax = High), position = position_dodge(width = .78), width = .12, linewidth = .55) +
    scale_fill_manual(values = pal_strata, drop = FALSE) +
    labs(title = title, x = xlab, y = "10-year risk (%)") + pub_theme()
}

incidence_panel <- function(d, title, xlab) {
  ggplot(d %>% mutate(Facet_group = factor(Facet_group, levels = c("low", "middle", "high")), Line_group = factor(Line_group, levels = c("low", "middle", "high"))),
         aes(time, risk_pct, color = Line_group, fill = Line_group)) +
    geom_ribbon(aes(ymin = risk_lower_pct, ymax = risk_upper_pct), alpha = .10, linewidth = 0) +
    geom_line(linewidth = .75) +
    facet_wrap(~Facet_group, nrow = 1) +
    scale_color_manual(values = pal_strata) + scale_fill_manual(values = pal_strata) +
    scale_x_continuous(limits = c(0, 10), breaks = seq(0, 10, 2.5), labels = c("0", "2.5", "5", "7.5", "10"), expand = expansion(mult = c(.01, .01))) +
    labs(title = title, x = xlab, y = "Cumulative incidence (%)") + pub_theme(11) +
    theme(plot.title = element_text(hjust = .5),
          panel.spacing.x = grid::unit(1.0, "lines"),
          panel.border = element_rect(color = "grey72", fill = NA, linewidth = .45))
}

score_heatmap <- function(m, title, rename_chns = FALSE) {
  m <- m[, vapply(m, function(x) all(is.na(x) | is.finite(suppressWarnings(as.numeric(x)))), logical(1)), drop = FALSE]
  if (rename_chns) names(m) <- normalize_chns_label(names(m))
  original_order <- names(m)
  keep <- intersect(diet_plot_order, original_order)
  row_index <- match(keep, original_order)
  m <- m[row_index, keep, drop = FALSE]
  rn <- keep
  mm <- as.matrix(dplyr::mutate(m, dplyr::across(dplyr::everything(), as.numeric)))
  if (nrow(mm) != ncol(mm)) stop("Concordance matrix is not square: ", title, call. = FALSE)
  dimnames(mm) <- list(rn, rn)
  d <- as_tibble(as.data.frame(as.table(mm)))
  names(d) <- c("Row", "Col", "rho")
  d$Row <- factor(as.character(d$Row), levels = rev(rn), labels = rev(unname(diet_display[rn])))
  d$Col <- factor(as.character(d$Col), levels = rn, labels = unname(diet_display[rn]))
  ggplot(d, aes(Col, Row, fill = rho)) +
    geom_tile(color = "white", linewidth = .35) +
    geom_text(aes(label = sprintf("%.2f", rho)), size = 2.8, fontface = "bold") +
    scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#D73027", midpoint = .5, limits = c(-1,1), oob = squish) +
    labs(title = title, x = NULL, y = NULL, fill = "Spearman\nrho") +
    guides(fill = guide_colorbar(barwidth = grid::unit(4.8, "cm"), barheight = grid::unit(.28, "cm"))) +
    pub_theme(11) + theme(axis.text.x = element_text(angle = 38, hjust = 1))
}

score_violin_panel <- function(curves, summary, title, chns = FALSE, show_y = TRUE) {
  curves <- curves %>% transmute(Diet = if (chns) normalize_chns_label(Diet) else as.character(Diet_score), score, density) %>%
    filter(Diet %in% diet_plot_order, is.finite(score), is.finite(density))
  sm <- if (chns) summary %>% transmute(Diet = normalize_chns_label(Diet), Center = mean, Low = pmax(0, mean-sd), High = pmin(100, mean+sd)) else
    summary %>% transmute(Diet = as.character(Diet_score), Center = Median, Low = P25, High = P75)
  ykey <- tibble(Diet = diet_plot_order, y0 = rev(seq_along(diet_plot_order)), label = unname(diet_display[diet_plot_order]))
  poly <- curves %>% inner_join(ykey, by = "Diet") %>% group_by(Diet) %>% mutate(v = .38*density/max(density, na.rm=TRUE)) %>%
    arrange(score, .by_group=TRUE) %>% group_modify(~ bind_rows(transmute(.x, score, yy=y0+v), transmute(arrange(.x, desc(score)), score, yy=y0-v))) %>% ungroup()
  sm <- sm %>% inner_join(ykey, by="Diet")
  ggplot(poly, aes(score, yy, group=Diet, fill=Diet)) +
    geom_polygon(alpha=.76, color=NA) +
    geom_segment(data=sm, aes(x=Low, xend=High, y=y0, yend=y0), inherit.aes=FALSE, linewidth=1.25, color="grey25") +
    geom_point(data=sm, aes(Center, y0), inherit.aes=FALSE, shape=21, size=2.1, stroke=.65, fill="white", color="grey20") +
    scale_fill_manual(values=pal_diet, guide="none") + scale_x_continuous(limits=c(0,100), breaks=seq(0,100,20), expand=expansion(mult=c(.01,.01))) +
    scale_y_continuous(breaks=ykey$y0, labels=if (show_y) ykey$label else NULL, expand=expansion(add=.55)) +
    labs(title=title, x="Harmonized dietary score (0–100)", y=NULL) + pub_theme(10) +
    theme(plot.title=element_text(hjust=.5), axis.text.y=if(show_y) element_text(face="bold") else element_blank(), axis.ticks.y=if(show_y) element_line() else element_blank())
}

score_range_panel <- function(d, title, chns = FALSE) {
  if (chns) {
    d <- d %>% transmute(Diet_score = normalize_chns_label(Diet), Center = mean,
                         Low = pmax(0, mean - sd), High = pmin(100, mean + sd))
    xlab <- "Harmonized score (mean +/- SD)"
  } else {
    d <- d %>% transmute(Diet_score, Center = Median, Low = P25, High = P75)
    xlab <- "Harmonized score (median and IQR)"
  }
  d <- d %>% filter(Diet_score %in% names(pal_diet)) %>%
    mutate(Diet_score = factor(Diet_score, levels = rev(names(pal_diet))))
  ggplot(d, aes(Center, Diet_score, color = Diet_score)) +
    geom_segment(aes(x = Low, xend = High, yend = Diet_score), linewidth = 1.15) +
    geom_point(size = 2.5) + scale_color_manual(values = pal_diet, drop = FALSE) +
    scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    labs(title = title, x = xlab, y = NULL) + pub_theme(10) + guides(color = "none")
}

joint_risk_panel <- function(d, title) {
  d <- d %>% mutate(
    dash = factor(tolower(dash), levels = c("low", "middle", "high")),
    maha = factor(tolower(maha), levels = c("low", "middle", "high"))
  )
  ggplot(d, aes(maha, dash, fill = risk10)) +
    geom_tile(color = "white", linewidth = .65) +
    geom_text(aes(label = sprintf("%.2f%%", risk10)), size = 4.3) +
    scale_x_discrete(labels = c(low = "Low", middle = "Middle", high = "High")) +
    scale_y_discrete(labels = c(low = "Low", middle = "Middle", high = "High")) +
    scale_fill_gradient(low = "#F7FBFF", high = "#08519C") +
    coord_fixed() +
    labs(title = title, x = "MAHA adherence", y = "DASH adherence", fill = "10-year risk (%)") +
    pub_theme(11) +
    theme(plot.title = element_text(face = "bold", hjust = .5, size = 13),
          axis.text = element_text(size = 10), axis.title = element_text(size = 11),
          legend.text = element_text(size = 9), legend.title = element_text(size = 9))
}

leaveout_panel <- function(d, title, xlab, estimate, low, high, label) {
  d <- d %>% filter(is.finite(.data[[estimate]]), is.finite(.data[[low]]), is.finite(.data[[high]])) %>%
    mutate(.lab = .data[[label]], .primary = grepl("Primary", .lab), .y = rev(seq_len(n())))
  ggplot(d) + geom_vline(xintercept = 1, linetype = 2, color = "#777777", linewidth = .45) +
    geom_segment(aes(x = .data[[low]], xend = .data[[high]], y = .y, yend = .y, color = .primary), linewidth = .8) +
    geom_point(aes(x = .data[[estimate]], y = .y, color = .primary), size = 2.6) +
    scale_color_manual(values = c(`FALSE` = "#4C78A8", `TRUE` = pal_diet[["MAHA"]]), guide = "none") +
    scale_y_continuous(breaks = d$.y, labels = d$.lab, expand = expansion(add = .55)) +
    labs(title = title, x = xlab, y = NULL) + pub_theme(11)
}

make_shared_score_figure <- function() {
  u <- wb_sheets(aux_path("ukb", "FigS3.score_distribution_concordance.out.xlsx"))
  n <- wb_sheets(aux_path("nhanes", "FigS3.score_distribution_concordance.out.xlsx"))
  c <- wb_sheets(aux_path("chns", "FigS7.score_distribution_concordance.out.xlsx"))
  pu <- score_violin_panel(get_sheet(u,"density_curve"), get_sheet(u,"score_distribution"), "a. Score distributions (UKB)", show_y=TRUE)
  pn <- score_violin_panel(get_sheet(n,"density_curve"), get_sheet(n,"score_distribution"), "b. Score distributions (NHANES)", show_y=FALSE)
  pc <- score_violin_panel(get_sheet(c,"density_curve"), get_sheet(c,"summary"), "c. Score distributions (CHNS)", chns=TRUE, show_y=FALSE)
  hu <- score_heatmap(get_sheet(u, "spearman_matrix"), "d. Concordance (UKB)") + theme(axis.text.y=element_text(face="bold"))
  hn <- score_heatmap(get_sheet(n, "spearman_matrix"), "e. Concordance (NHANES)") + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
  hc <- score_heatmap(get_sheet(c, "spearman"), "f. Concordance (CHNS)", rename_chns = TRUE) + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
  figS1 <- align_publication_grid(list(list(pu, pn, pc), list(hu, hn, hc)),
                                  rel_widths = c(1, .8, .8), rel_heights = c(.82, 1.18), row_gap = .14)
  save_pub(figS1, "FigS1.score_distribution_concordance.png", 14, 12.4)
  write_book(list(ukb_distribution = get_sheet(u, "score_distribution"), ukb_spearman = get_sheet(u, "spearman_matrix"),
                  nhanes_distribution = get_sheet(n, "score_distribution"), nhanes_spearman = get_sheet(n, "spearman_matrix"),
                  chns_distribution = get_sheet(c, "summary"), chns_spearman = get_sheet(c, "spearman")),
             "FigS1.score_distribution_concordance.out.xlsx")
}

profile_panel <- function(d, title) {
  ggplot(d, aes(decile, mean_score, color = Component, group = Component)) +
    geom_line(linewidth = .75) + geom_point(size = 1.6) +
    scale_x_continuous(expand = expansion(mult = c(.02, .04))) +
    labs(title = title, x = "MAHA decile", y = "Mean component score") +
    pub_theme(11) + theme(legend.position = "right", legend.text = element_text(size = 7, face = "bold"), legend.key.width = grid::unit(10, "pt"))
}

matched_null_panel <- function(direction_null, summary_row, metric, observed_col, title, xlab) {
  null <- direction_null %>% filter(!is_prespecified_MAHA) %>% filter(is.finite(.data[[metric]]))
  obs <- as.numeric(summary_row[[observed_col]][1])
  ggplot(null, aes(.data[[metric]])) +
    geom_histogram(bins = 36, fill = "#9ECAE1", color = "white", linewidth = .2) +
    geom_vline(xintercept = obs, color = pal_diet[["MAHA"]], linewidth = 1.1) +
    annotate("text", x = obs, y = Inf, label = "prespecified MAHA", color = pal_diet[["MAHA"]],
             hjust = 1.05, vjust = 1.15, fontface = "bold", size = 3.0) +
    coord_cartesian(clip = "off") +
    labs(title = title, x = xlab, y = "Matched arbitrary scores") + pub_theme(11)
}

make_fig1 <- function() {
  cohort_box <- tibble(
    x = c(1, 2, 3), y = 3,
    label = c("UK Biobank\nMain cohort\nIncident disease + mortality + geography",
              "NHANES\nExternal validation\nMortality",
              "CHNS\nExternal validation\nMortality + construct checks")
  )
  mid_box <- tibble(x = c(1, 2, 3), y = 2,
                    label = c("MAHA\nPrimary + sensitivity specifications", "Comparators\nDASH • MIND • MEDI", "Validation\nConvergence • robustness • QC"))
  low_box <- tibble(x = c(1, 2, 3), y = 1,
                    label = c("Fig 2\nUKB precision & risk", "Fig 3\nNHANES mortality", "Fig 4\nCHNS mortality"))
  boxes <- bind_rows(cohort_box, mid_box, low_box)
  p <- ggplot() +
    geom_tile(data = boxes, aes(x, y), width = .78, height = .55, fill = "white", color = "#303030", linewidth = .6) +
    geom_text(data = boxes, aes(x, y, label = label), size = 4, lineheight = .98) +
    geom_segment(data = tibble(x = rep(1:3, 2), y = rep(c(2.70, 1.70), each = 3), yend = rep(c(2.31, 1.31), each = 3)),
                 aes(x = x, xend = x, y = y, yend = yend), arrow = grid::arrow(length = unit(.16, "cm")), linewidth = .55) +
    annotate("segment", x = 2.62, xend = 3.00, y = 2, yend = 2, arrow = grid::arrow(length = unit(.15, "cm")), linewidth = .5) +
    annotate("text", x = 2, y = .28,
             label = "Primary question: does a higher MAHA score show the same consistent protection as established dietary scores?",
             fontface = "bold", size = 4.5) +
    coord_cartesian(xlim = c(.45, 3.55), ylim = c(.05, 3.45), clip = "off") +
    labs(title = "Fig 1. Overall study design",
         subtitle = "Three cohorts; one prespecified comparison framework; cohort-specific outcome validation") +
    theme_void(base_size = 14) +
    theme(plot.title = element_text(face = "bold", size = 18), plot.subtitle = element_text(size = 12))
}

make_ukb <- function() {
  wb2 <- wb_sheets(aux_path("ukb", "Fig2.prospective.out.xlsx"))
  wb3 <- wb_sheets(aux_path("ukb", "Fig3.precision_attenuation.out.xlsx"))

  main <- get_sheet(wb3, "Fig3_main_summary") %>%
    bind_cols(parse_ci(.$HR_95CI) %>% rename(HR = est, LCI = lo, UCI = hi))
  maha <- get_sheet(wb3, "Fig3_MAHA_summary") %>%
    bind_cols(parse_ci(.$HR_90CI) %>% rename(HR = est, LCI90 = lo, UCI90 = hi))
  h2h <- get_sheet(wb3, "Fig3_head2head") %>% mutate(Comp2 = sub("^MAHA vs ", "", Comparison))
  mde <- stringr::str_match(as.character(maha$MDE80_HR_label), "([0-9.]+)[^0-9]+([0-9.]+)")
  maha <- maha %>% mutate(MDE_HR_lower = suppressWarnings(as.numeric(mde[,2])), MDE_HR_upper = suppressWarnings(as.numeric(mde[,3])))

  p1 <- forest_multi(main, "a. Main associations", "Hazard ratio per 1 SD healthier score",
                     estimate = "HR", low = "LCI", high = "UCI", diets = primary_diets)
  p2 <- forest_multi(h2h %>% transmute(Outcome, Diet = Comp2, HR = HR_ratio_maha_vs_trad, LCI = LCI_ratio, UCI = UCI_ratio),
                     "b. Direct attenuation versus established scores", "HR ratio: MAHA / comparator",
                     estimate = "HR", low = "LCI", high = "UCI", diets = c("DASH", "MIND", "MEDI"))

  eq <- maha %>% transmute(Outcome, Diet = "MAHA", HR, LCI = LCI90, UCI = UCI90,
                           lab = paste0("TOST P=", ap_p_lab(TOST_p), "; BF₀₁",
                                        if_else(!is.na(BF01) & BF01 < .001, " < .001", paste0("=", ap_bf_lab(BF01)))))
  p3 <- forest_single(eq, "c. MAHA equivalence and Bayes evidence", "MAHA hazard ratio with 90% CI",
                      estimate = "HR", low = "LCI", high = "UCI", label = "Outcome", color = "Outcome",
                      zone = c(1/1.05, 1.05), right_label = "lab") + guides(color = "none")
  pow <- maha %>% transmute(Outcome, Diet = "MAHA", MDE = MDE_HR_upper, LCI = 1, UCI = MDE_HR_upper,
                            lab = sprintf("%.2f–%.2f", MDE_HR_lower, MDE_HR_upper))
  p4 <- forest_single(pow, "d. MAHA detectable effect", "Minimum detectable HR at 80% power",
                      estimate = "MDE", low = "LCI", high = "UCI", label = "Outcome", color = "Outcome",
                      ref = 1.05, right_label = "lab") + guides(color = "none")
  p5 <- risk_bar(get_sheet(wb2, "Fig2_10y_risk_by_DASH"), "e. 10-year mortality risk by DASH strata", "DASH strata")
  p6 <- risk_bar(get_sheet(wb2, "Fig2_10y_risk_by_MAHA"), "f. 10-year mortality risk by MAHA strata", "MAHA strata")
  fig2 <- align_publication_grid(list(list(p1, p2), list(p3, p4), list(p5, p6)),
                                 rel_widths = c(1, 1), rel_heights = c(1, 1, .95))
  save_pub(fig2, "Fig2.ukb.precision.png", 16, 16.5)
  write_book(list(main_summary = main, head2head = h2h, MAHA_precision = maha,
                  risk_by_DASH = get_sheet(wb2, "Fig2_10y_risk_by_DASH"),
                  risk_by_MAHA = get_sheet(wb2, "Fig2_10y_risk_by_MAHA")),
             "Fig2.ukb.precision.out.xlsx")

  if (copy_if_exists(aux_path("ukb", "Fig1.phewas.png"), pub_path("Fig1.ukb.phewas.png"))) {
    copy_if_exists(aux_path("ukb", "Fig1.phewas.out.xlsx"), pub_path("Fig1.ukb.phewas.out.xlsx"))
  }

  alt <- get_sheet(wb2, "Alternative_MAHA_assoc") %>% mutate(Diet = recode(Diet,
    balanced = "MAHA-balanced", strict = "MAHA-strict", `no dairy` = "MAHA-no dairy", `no protein` = "MAHA-no protein", .default = Diet))
  ukbmain <- get_sheet(wb2, "Fig2_main_assoc")
  pS2a <- forest_multi(ukbmain, "a. Prospective associations", "Hazard ratio per 1 SD", diets = primary_diets)
  pS2b <- forest_multi(alt, "b. Alternative MAHA specifications", "Hazard ratio per 1 SD",
                       diets = c("MAHA-balanced", "MAHA-strict", "MAHA-no dairy", "MAHA-no protein"),
                       hollow = c("MAHA-balanced", "MAHA-strict", "MAHA-no dairy", "MAHA-no protein"))
  pS2c <- incidence_panel(get_sheet(wb2, "Fig2_incidence_by_DASH"), "c. Cumulative incidence by DASH strata", "Years")
  pS2d <- incidence_panel(get_sheet(wb2, "Fig2_incidence_by_MAHA"), "d. Cumulative incidence by MAHA strata", "Years")
  figS2 <- align_publication_grid(list(list(pS2a, pS2b), list(pS2c, pS2d)),
                                  rel_widths = c(1, 1), rel_heights = c(1, 1))
  save_pub(figS2, "FigS3.ukb.prospective_sensitivity.png", 16, 12.4)
  write_book(list(main_assoc = ukbmain, alternative_MAHA = alt,
                  discordant = get_sheet(wb2, "Fig2_main_discordant"), AIC = get_sheet(wb2, "Fig2_main_aic"),
                  incidence_DASH = get_sheet(wb2, "Fig2_incidence_by_DASH"), incidence_MAHA = get_sheet(wb2, "Fig2_incidence_by_MAHA"),
                  strata_counts = get_sheet(wb2, "Fig2_strata_counts")), "FigS3.ukb.prospective_sensitivity.out.xlsx")

  copy_if_exists(aux_path("ukb", "FigS2.phewas_head2head.png"), pub_path("FigS2.ukb.phewas_head2head.png"))
  copy_if_exists(aux_path("ukb", "FigS2.phewas_head2head.out.xlsx"), pub_path("FigS2.ukb.phewas_head2head.out.xlsx"))
  jr_path <- aux_path("ukb", "FigS4.joint_risk.out.xlsx")
  jr <- read_excel(jr_path) %>% as_tibble()
  jr_titles <- intersect(rev(ukb_outcome_plot_levels), unique(as.character(jr$Outcome)))
  jr_panels <- lapply(jr_titles, function(z) joint_risk_panel(jr %>% filter(Outcome == z), z))
  fig_joint <- wrap_plots(jr_panels, ncol = 3) +
    plot_annotation(title = "Joint 10-year risk across MAHA and DASH adherence strata",
                    theme = theme(plot.title = element_text(face="bold", hjust=.5, size=16, margin=margin(b=10))))
  save_pub(fig_joint, "FigS4.ukb.joint_risk.png", 13.5, 11.8, dpi = 400)
  write_book(list(joint_risk = jr), "FigS4.ukb.joint_risk.out.xlsx")

  diag_png <- aux_path("ukb", "FigS5.survival_diagnostics.png")
  diag_xlsx <- aux_path("ukb", "FigS5.survival_diagnostics.out.xlsx")
  if (!file.exists(diag_png) || !file.exists(diag_xlsx)) {
    stop("UKB survival diagnostics are missing. Run ./maha.sh ukb --steps survdiag first.", call. = FALSE)
  }
  copy_if_exists(diag_png, pub_path("FigS5.ukb.survival_diagnostics.png"))
  copy_if_exists(diag_xlsx, pub_path("FigS5.ukb.survival_diagnostics.out.xlsx"))

  # These auxiliary files used the pre-renumbering figure numbers.  Remove them
  # only after the replacement inputs and publication outputs have been checked.
  obsolete_ukb_aux <- aux_path("ukb", c(
    "FigS4.phewas_head2head.png", "FigS4.phewas_head2head.out.xlsx",
    "FigS5.joint_risk.png", "FigS5.joint_risk.out.xlsx",
    "UKB_survival_diagnostics.png", "UKB_survival_diagnostics.out.xlsx"
  ))
  replacement_outputs <- pub_path(c(
    "FigS2.ukb.phewas_head2head.png", "FigS2.ukb.phewas_head2head.out.xlsx",
    "FigS4.ukb.joint_risk.png", "FigS4.ukb.joint_risk.out.xlsx",
    "FigS5.ukb.survival_diagnostics.png", "FigS5.ukb.survival_diagnostics.out.xlsx"
  ))
  if (all(file.exists(replacement_outputs)) && all(file.info(replacement_outputs)$size > 1000)) {
    unlink(obsolete_ukb_aux[file.exists(obsolete_ukb_aux)])
  }

  v <- wb_sheets(aux_path("ukb", "FigS1.validate.out.xlsx")); cp <- wb_sheets(aux_path("ukb", "FigS2.construct_profile.out.xlsx")); sc <- wb_sheets(aux_path("ukb", "FigS3.score_distribution_concordance.out.xlsx"))
  vs <- get_sheet(v, "summary"); vn <- get_sheet(v, "direction_null"); vl <- get_sheet(v, "leave_one_component_out")
  pprimarya <- matched_null_panel(vn, vs, "protective_Z", "MAHA_protective_Z", "a. Mortality criterion validity", "Protective Z statistic (-beta/SE)")
  pprimaryb <- matched_null_panel(vn, vs, "rho_established_anchor", "MAHA_anchor_Spearman", "b. Convergent construct validity", "Spearman rho with established-diet anchor")
  pprimaryc <- leaveout_panel(vl, "c. Leave-one-component-out mortality", "HR per 1 SD healthier score",
                         "HR", "CI_low", "CI_high", "label")
  pprimaryd <- profile_panel(get_sheet(cp, "decile_profile"), "d. Component profile across MAHA deciles")
  figprimary <- align_publication_grid(list(list(pprimarya, pprimaryb), list(pprimaryc, pprimaryd)), rel_widths = c(1, 1))
  save_pub(figprimary, "FigS8_tmp.ukb.construct_stress_tests.png", 16, 11.6)
  write_book(list(validation_summary = vs, direction_null = vn, leave_one_out = vl,
                  decile_profile = get_sheet(cp, "decile_profile"), component_metrics = get_sheet(cp, "component_metrics"),
                  convergent_validity = get_sheet(cp, "convergent_validity")), "FigS8_tmp.ukb.construct_stress_tests.out.xlsx")

  copy_if_exists(aux_path("ukb", "Fig5.geography.png"), pub_path("Fig5.ukb.geography.png"))
  copy_if_exists(aux_path("ukb", "Fig5.geography.out.xlsx"), pub_path("Fig5.ukb.geography.out.xlsx"))

  copy_if_exists(aux_path("ukb", "Table1.out.xlsx"), pub_path("Table1.ukb.out.xlsx"))
}

make_nhanes <- function() {
  wb <- wb_sheets(aux_path("nhanes", "Fig4.mortality_validation.out.xlsx"))
  a <- get_sheet(wb, "Fig4A_assoc")
  h <- get_sheet(wb, "Fig4G_head2head")
  eq <- get_sheet(wb, "Fig4D_equiv") %>% bind_cols(parse_ci(.$HR_90CI) %>% rename(HR = est, LCI = lo, UCI = hi))
  bf <- get_sheet(wb, "Fig4E_BF")
  pw <- get_sheet(wb, "Fig4F_power")
  eqbf <- eq %>% left_join(bf %>% select(Diet, BF01), by = "Diet") %>% mutate(lab = paste0("TOST P=", ap_p_lab(TOST_p), "; BF₀₁=", ap_bf_lab(BF01)))

  p1 <- forest_single(a %>% select(Diet, estimate, conf.low, conf.high) %>% mutate(Outcome = "All-cause mortality"),
                      "a. Mortality associations", "HR per 1 SD healthier score", estimate = "estimate", low = "conf.low", high = "conf.high",
                      label = "Diet", color = "Diet", hollow = setdiff(maha_variants, "MAHA"))
  p2 <- forest_single(h %>% transmute(Diet = Comp2, Comparison, HR = HR_ratio_maha_vs_trad, LCI = LCI_ratio, UCI = UCI_ratio),
                      "b. Direct attenuation versus established scores", "HR ratio: MAHA / comparator", estimate = "HR", low = "LCI", high = "UCI", label = "Comparison", color = "Diet")
  p3 <- forest_single(eqbf, "c. MAHA equivalence and Bayes evidence", "Hazard ratio with 90% CI", estimate = "HR", low = "CI90_low_HR", high = "CI90_high_HR",
                      label = "Diet", color = "Diet", hollow = setdiff(maha_variants, "MAHA"), right_label = "lab", zone = c(1/1.05,1.05))
  p4 <- forest_single(pw %>% mutate(MDE = MDE_HR_upper, LCI = 1, UCI = MDE_HR_upper, lab = MDE80_HR_label),
                      "d. Detectable effect", "Minimum detectable HR at 80% power", ref = 1.05,
                      estimate = "MDE", low = "LCI", high = "UCI", label = "Diet", color = "Diet",
                      hollow = setdiff(maha_variants, "MAHA"), right_label = "lab")
  p5 <- risk_bar(get_sheet(wb, "Fig4B_risk_DASH"), "e. 10-year mortality risk by DASH strata", "DASH strata")
  p6 <- risk_bar(get_sheet(wb, "Fig4C_risk_MAHA"), "f. 10-year mortality risk by MAHA strata", "MAHA strata")
  fig3 <- align_publication_grid(list(list(p1, p2), list(p3, p4), list(p5, p6)), rel_widths = c(1, 1))
  save_pub(fig3, "Fig3.nhanes.mortality.png", 16, 16.5)
  write_book(list(associations = a, head2head = h, equivalence = eq, bayes = bf, power = pw,
                  risk_by_DASH = get_sheet(wb, "Fig4B_risk_DASH"), risk_by_MAHA = get_sheet(wb, "Fig4C_risk_MAHA")),
             "Fig3.nhanes.mortality.out.xlsx")

  v <- wb_sheets(aux_path("nhanes", "FigS1.validate.out.xlsx")); cp <- wb_sheets(aux_path("nhanes", "FigS2.construct_profile.out.xlsx")); sc <- wb_sheets(aux_path("nhanes", "FigS3.score_distribution_concordance.out.xlsx"))
  vs <- get_sheet(v, "summary"); vn <- get_sheet(v, "direction_null"); vl <- get_sheet(v, "leave_one_component_out")
  pS4a <- matched_null_panel(vn, vs, "protective_Z", "MAHA_protective_Z", "a. Mortality criterion validity", "Protective Z statistic (-beta/SE)")
  pS4b <- matched_null_panel(vn, vs, "rho_established_anchor", "MAHA_anchor_Spearman", "b. Convergent construct validity", "Spearman rho with established-diet anchor")
  pS4c <- leaveout_panel(vl, "c. Leave-one-component-out mortality", "Survey-weighted HR per 1 SD", "HR", "CI_low", "CI_high", "label")
  pS4d <- ggplot(vl, aes(rho_primary, HR, color = omitted)) +
    geom_hline(yintercept = 1, linetype = 2, color = "grey55") +
    geom_vline(xintercept = 1, linetype = 2, color = "grey55") +
    geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = .008, linewidth = .65) +
    geom_point(size = 2.8) +
    labs(title = "d. Rank stability and mortality effect after omission",
         x = "Spearman rho with primary MAHA", y = "Mortality HR per 1 SD", color = "Omitted domain") +
    pub_theme(11) + theme(legend.position = "right", legend.text = element_text(size = 7, face = "bold"))
  figS4 <- align_publication_grid(list(list(pS4a, pS4b), list(pS4c, pS4d)), rel_widths = c(1, 1))
  save_pub(figS4, "FigS9_tmp.nhanes.validation.png", 16, 11.6)
  write_book(list(summary = vs, direction_null = vn, leave_one_out = vl), "FigS9_tmp.nhanes.validation.out.xlsx")

  pS5a <- profile_panel(get_sheet(cp, "decile_profile"), "a. Component profile across MAHA deciles")
  cm <- get_sheet(cp, "component_metrics")
  pS5b <- ggplot(cm, aes(D10_minus_D1, reorder(Component, D10_minus_D1), color = Component)) + geom_segment(aes(x = 0, xend = D10_minus_D1, yend = reorder(Component, D10_minus_D1)), linewidth = .7) + geom_point(size = 2.6) +
    labs(title = "b. Component separation", x = "Decile 10 – decile 1", y = NULL) + pub_theme(11) + guides(color = "none")
  sd <- get_sheet(sc, "score_distribution")
  pS5c <- ggplot(sd, aes(Median, factor(Diet_score, levels = rev(Diet_score)), color = Diet_score)) + geom_segment(aes(x = P25, xend = P75, yend = factor(Diet_score, levels = rev(Diet_score))), linewidth = 1.2) + geom_point(size = 2.7) +
    scale_color_manual(values = pal_diet) + labs(title = "c. Score distributions", x = "Harmonized score (median and IQR)", y = NULL) + pub_theme(11)
  pS5d <- score_heatmap(get_sheet(sc, "spearman_matrix"), "d. Spearman concordance")
  figS5 <- align_publication_grid(list(list(pS5a, pS5b), list(pS5c, pS5d)), rel_widths = c(1, 1))
  save_pub(figS5, "FigS10_tmp.nhanes.construct_concordance.png", 16, 12.1)
  write_book(list(decile_profile = get_sheet(cp, "decile_profile"), component_metrics = cm, convergent_validity = get_sheet(cp, "convergent_validity"),
                  score_distribution = sd, spearman_matrix = get_sheet(sc, "spearman_matrix")), "FigS10_tmp.nhanes.construct_concordance.out.xlsx")
  copy_if_exists(aux_path("nhanes", "Table1.out.xlsx"), pub_path("Table1.nhanes.out.xlsx"))
}

make_chns <- function() {
  # Remove the superseded construct-validation FigS9 so a subsequent
  # `all --pub-only` run has exactly one numbered CHNS FigS9.
  stale_root <- pub_path(c("FigS9.chns.validation.png","FigS9.chns.validation.out.xlsx"))
  unlink(stale_root[file.exists(stale_root)],force=TRUE)

  pubin_path <- aux_path("chns","publication_inputs.out.xlsx")
  if(!file.exists(pubin_path)) stop("CHNS publication inputs are missing. Run ./maha.sh chns first: ",pubin_path,call.=FALSE)
  wb <- wb_sheets(pubin_path)
  sens <- get_sheet(wb,"sensitivity") %>% mutate(Diet=normalize_chns_label(Diet))
  benchmark <- get_sheet(wb,"benchmark",required=FALSE)
  meta <- get_sheet(wb,"analysis_metadata",required=FALSE)

  focused_path <- aux_path("chns","mortality_models_final.out.xlsx")
  if(!file.exists(focused_path))
    stop("Focused CHNS high-event results are missing. Run ./maha.sh chns first: ",focused_path,call.=FALSE)
  fw <- wb_sheets(focused_path)
  selected <- get_sheet(fw,"primary_results")
  selected_h2h <- get_sheet(fw,"head2head")
  selected_defs <- get_sheet(fw,"model_definitions")
  selected_pattern <- get_sheet(fw,"pattern_diagnostic")
  last_contact <- get_sheet(fw,"last_contact_results")
  maha_algorithms <- get_sheet(fw,"maha_algorithm_results")
  model1_id <- "Model1_default_reference"
  model4_id <- "Model4_socioeconomic"
  model5_id <- "Model5_bmi"
  expected_models <- c(model1_id,model4_id,model5_id)
  found_models <- sort(unique(as.character(selected$Model)))
  if(!setequal(found_models,expected_models) || length(found_models)>6L)
    stop("Focused CHNS workbook must contain only Models 1, 4 and 5; found: ",
         paste(found_models,collapse=", "),call.=FALSE)
  model_events <- selected %>% group_by(Model) %>% summarise(Events=first(Events),.groups="drop")
  event_n <- function(model_id) {
    z <- model_events$Events[model_events$Model==model_id]
    if(length(z)==1L&&is.finite(z)) format(z,big.mark=",") else "unknown"
  }

  p1 <- forest_single(
    selected %>% filter(Model==model1_id,Diet%in%primary_diets) %>%
      mutate(Diet=factor(Diet,levels=rev(primary_diets))),
    paste0("a. Primary fixed-wave analysis (",event_n(model1_id)," deaths)"),
    "Adjusted HR per 1 SD healthier score",
    estimate="HR",low="CI_low",high="CI_high",label="Diet",color="Diet"
  )

  p2 <- forest_single(
    selected %>% filter(Model==model4_id,Diet%in%primary_diets) %>%
      mutate(Diet=factor(Diet,levels=rev(primary_diets))),
    paste0("b. High-event cumulative-average analysis (",event_n(model4_id)," deaths)"),
    "Adjusted HR per 1 SD healthier score",
    estimate="HR",low="CI_low",high="CI_high",label="Diet",color="Diet"
  )

  p3 <- forest_single(
    selected_h2h %>% filter(Model==model1_id,Comparator%in%c("DASH","MIND","MEDI")) %>%
      transmute(Diet=sub("-reference$","",Comparator),Comparison,
                Ratio=HR_ratio_MAHA_vs_comparator,LCI=CI_low_ratio,UCI=CI_high_ratio),
    "c. Primary direct coefficient comparison","HR ratio: MAHA / comparator",
    estimate="Ratio",low="LCI",high="UCI",label="Comparison",color="Diet"
  )

  p4 <- forest_single(
    selected_h2h %>% filter(Model==model4_id,Comparator%in%c("DASH-reference","MIND","MEDI")) %>%
      transmute(Diet=sub("-reference$","",Comparator),Comparison=sub("-reference", "", Comparison),
                Ratio=HR_ratio_MAHA_vs_comparator,LCI=CI_low_ratio,UCI=CI_high_ratio),
    "d. High-event direct coefficient comparison","HR ratio: MAHA / comparator",
    estimate="Ratio",low="LCI",high="UCI",label="Comparison",color="Diet"
  )

  fig4 <- align_publication_grid(list(list(p1,p2),list(p3,p4)),rel_widths=c(1,1))
  save_pub(fig4,"Fig4.chns.mortality.png",16,11.8)
  out4 <- list(primary_results=selected%>%filter(Model==model1_id),
               high_event_results=selected%>%filter(Model==model4_id),
               primary_head2head=selected_h2h%>%filter(Model==model1_id),
               high_event_head2head=selected_h2h%>%filter(Model==model4_id),
               model_definitions=selected_defs,pattern_diagnostic=selected_pattern)
  if(!is.null(benchmark)) out4$benchmark <- benchmark
  if(!is.null(meta)) out4$analysis_metadata <- meta
  write_book(out4,"Fig4.chns.mortality.out.xlsx")

  # CHNS construct validation.
  v <- wb_sheets(aux_path("chns","FigS1.validate.out.xlsx"))
  cp <- wb_sheets(aux_path("chns","FigS2.construct_profile.out.xlsx"))
  vs <- get_sheet(v,"summary")
  vn <- get_sheet(v,"direction_null")
  vl <- get_sheet(v,"leave_one_out")

  pS9a <- ggplot(vn %>% filter(!is_prespecified,is.finite(rho_established)),aes(rho_established))+
    geom_histogram(bins=36,fill="#9ECAE1",color="white")+
    geom_vline(xintercept=vs$Established_anchor_rho[1],color=pal_diet[["MAHA"]],linewidth=1.1)+
    labs(title="a. Established-score convergence",x="Spearman rho with established-score anchor",y="Matched arbitrary scores")+pub_theme(11)

  pS9b <- ggplot(vn %>% filter(!is_prespecified,is.finite(rho_preference)),aes(rho_preference))+
    geom_histogram(bins=36,fill="#FDD0A2",color="white")+
    geom_vline(xintercept=vs$Preference_rho[1],color=pal_diet[["MAHA"]],linewidth=1.1)+
    labs(title="b. Independent food-preference criterion",x="Spearman rho with healthy preference index",y="Matched arbitrary scores")+pub_theme(11)

  vl <- vl %>% mutate(Label=ifelse(is.na(omitted),"Primary MAHA",paste0("Without ",omitted)))
  pS9c <- leaveout_panel(vl,"c. Leave-one-domain-out mortality","Adjusted mortality HR per 1 SD",
                         "mortality_HR","mortality_CI_low","mortality_CI_high","Label")
  pS9d <- profile_panel(get_sheet(cp,"decile_profile"),"d. Component profile across MAHA deciles")
  figS9 <- align_publication_grid(list(list(pS9a,pS9b),list(pS9c,pS9d)),rel_widths=c(1,1))
  save_pub(figS9,"Exploratory_CHNS_construct_validation.png",16,12.1)
  write_book(list(summary=vs,direction_null=vn,leave_one_out=vl,
                  decile_profile=get_sheet(cp,"decile_profile"),
                  component_metrics=get_sheet(cp,"component_metrics"),
                  convergent_validity=get_sheet(cp,"convergent_validity")),
             "Exploratory_CHNS_construct_validation.out.xlsx")

  # External construct / clinical validation.
  fp <- wb_sheets(aux_path("chns","Figprimary.food_preference_validation.out.xlsx"))
  cv <- wb_sheets(aux_path("chns","FigS6.cardiometabolic_validation.out.xlsx"))
  sc <- wb_sheets(aux_path("chns","FigS7.score_distribution_concordance.out.xlsx"))
  fa <- get_sheet(fp,"adjusted_associations") %>% mutate(CI_low=beta-1.96*se,CI_high=beta+1.96*se,Diet="MAHA")
  pEa <- forest_single(fa,"a. Adjusted food-preference associations","Change in preference per 1 SD MAHA",
                       ref=0,estimate="beta",low="CI_low",high="CI_high",label="Outcome",color="Diet")+
    scale_color_manual(values=c(MAHA="#4C78A8"))+guides(color="none")
  fc <- get_sheet(fp,"correlations")
  pEb <- ggplot(fc,aes(Spearman_rho,reorder(Criterion,Spearman_rho)))+
    geom_vline(xintercept=0,linetype=2)+
    geom_segment(aes(x=0,xend=Spearman_rho,yend=reorder(Criterion,Spearman_rho)),color="#E69F00",linewidth=.75)+
    geom_point(color="#E69F00",size=2.6)+
    labs(title="b. Independent preference convergence",x="Spearman rho",y=NULL)+pub_theme(11)
  pi <- get_sheet(cv,"prospective_incident") %>% mutate(Diet=normalize_chns_label(Score)) %>% filter(Diet %in% primary_diets)
  pEc <- forest_multi(pi,"c. Prospective cardiometabolic validation","Adjusted OR per 1 SD healthier score",
                       estimate="OR",low="CI_low",high="CI_high",outcome="Outcome",diet="Diet",diets=primary_diets)
  sp <- get_sheet(sc,"spearman")
  pEd <- score_heatmap(sp,"d. Score concordance",rename_chns=TRUE)
  figE <- align_publication_grid(list(list(pEa,pEb),list(pEc,pEd)),rel_widths=c(1,1))
  save_pub(figE,"Exploratory_CHNS_external_validation.png",16,12.6)
  write_book(list(preference_adjusted=fa,preference_correlations=fc,prospective_incident=pi,spearman=sp),
             "Exploratory_CHNS_external_validation.out.xlsx")

  # Focused high-event sensitivity figure: only the retained Model 4/5 family.
  pH1 <- forest_single(
    selected %>% filter(Model==model5_id,Diet%in%primary_diets) %>%
      mutate(Diet=factor(Diet,levels=rev(primary_diets))),
    "a. High-event model with BMI adjustment","Adjusted HR per 1 SD healthier score",
    estimate="HR",low="CI_low",high="CI_high",label="Diet",color="Diet"
  )
  pH2 <- forest_single(
    last_contact %>% filter(Adjustment==model4_id,Diet%in%primary_diets) %>%
      mutate(Diet=factor(Diet,levels=rev(primary_diets))),
    "b. Socioeconomic model: censor at last contact","Adjusted HR per 1 SD healthier score",
    estimate="HR",low="CI_low",high="CI_high",label="Diet",color="Diet"
  )
  pH3 <- forest_single(
    last_contact %>% filter(Adjustment==model5_id,Diet%in%primary_diets) %>%
      mutate(Diet=factor(Diet,levels=rev(primary_diets))),
    "c. BMI model: censor at last contact","Adjusted HR per 1 SD healthier score",
    estimate="HR",low="CI_low",high="CI_high",label="Diet",color="Diet"
  )
  alg_plot <- maha_algorithms %>% filter(Model%in%c(model4_id,model5_id)) %>%
    mutate(
      Model_label=recode(Model,
        Model4_socioeconomic="Socioeconomic model",
        Model5_bmi="+ BMI model"),
      Score_algorithm=factor(Score_algorithm,levels=c("MAHA","MAHA-density","MAHA-residual"))
    )
  pH4 <- ggplot(alg_plot,aes(HR,Model_label,color=Score_algorithm))+
    geom_vline(xintercept=1,linetype=2,color="grey55")+
    geom_errorbar(aes(xmin=CI_low,xmax=CI_high),orientation="y",width=0,
                  position=position_dodge(width=.62),linewidth=.65)+
    geom_point(position=position_dodge(width=.62),size=2.7)+
    scale_color_manual(values=pal_diet[c("MAHA","MAHA-density","MAHA-residual")],drop=FALSE)+
    labs(title="d. MAHA energy-adjustment algorithms",x="Adjusted HR per 1 SD",y=NULL,color=NULL)+
    pub_theme(11)
  figS9_focused <- align_publication_grid(list(list(pH1,pH2),list(pH3,pH4)),rel_widths=c(1,1))
  save_pub(figS9_focused,"FigS7_tmp.chns.high_event_sensitivity.png",16,11.8)
  write_book(list(
    model5_end_wave=selected%>%filter(Model==model5_id),
    model4_last_contact=last_contact%>%filter(Adjustment==model4_id),
    model5_last_contact=last_contact%>%filter(Adjustment==model5_id),
    maha_algorithms=maha_algorithms,
    model_definitions=selected_defs
  ),"FigS7_tmp.chns.high_event_sensitivity.out.xlsx")

  # Same-cohort mortality sensitivities only.
  counts <- sens %>% filter(Diet=="MAHA") %>% distinct(Analysis,N,Events)
  analysis_labels <- counts %>%
    transmute(Analysis,Outcome=paste0(Analysis,"\n(mortality N = ",Events,")"))
  splot <- sens %>%
    filter(Diet %in% primary_diets,is.finite(HR)) %>%
    left_join(analysis_labels,by="Analysis")
  figS10 <- forest_multi(
    splot,"Same-cohort mortality sensitivity analyses","Adjusted HR per 1 SD healthier score",
    estimate="HR",low="CI_low",high="CI_high",outcome="Outcome",diet="Diet",diets=primary_diets
  )
  save_pub(figS10,"FigS6_tmp.chns.mortality_sensitivity.png",10.5,7.5)

  qc <- wb_sheets(aux_path("chns","FigS8.food_mapping_qc.out.xlsx"))
  fu <- wb_sheets(aux_path("chns","FigS9.followup_qc.out.xlsx"))
  tb <- wb_sheets(aux_path("chns","Table1.out.xlsx"))
  write_book(list(
    sensitivity=sens,
    followup_overall=get_sheet(fu,"overall"),
    cohort_table=get_sheet(tb,"Table1"),
    foodcode_field_coverage=get_sheet(qc,"foodcode_field_coverage",required=FALSE),
    metadata=get_sheet(tb,"metadata"),
    input_sources=get_sheet(tb,"input_sources")
  ),"FigS6_tmp.chns.mortality_sensitivity.out.xlsx")

  copy_if_exists(aux_path("chns","Table1.out.xlsx"),pub_path("Table1.chns.out.xlsx"))
}

if (scope == "all") make_shared_score_figure()
if (scope %in% c("ukb", "all")) make_ukb()
if (scope %in% c("nhanes", "all")) make_nhanes()
if (scope %in% c("chns", "all")) make_chns()
finalize_supplement_renumbering(scope)
if (scope == "all") {
  assert_contiguous_figure_numbers("Fig", 1:5)
  assert_contiguous_figure_numbers("FigS", 1:10)
}

message("[publication] complete: ", outdir)
