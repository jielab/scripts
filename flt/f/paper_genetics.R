source(file.path(Sys.getenv("FLT_CODE"), "runtime.R"))
library(data.table)
library(ggplot2)
library(patchwork)
out <- Sys.getenv("FLT_PAPER")
meta <- as.data.table(readxl::read_excel(Sys.getenv("FLT_DICTIONARY"), sheet = "item182"))
setnames(meta, "gwas_id", "trait1")
registry <- as.data.table(flt_panel())
rg_path <- file.path(flt_run("ldsc"), "ldsc_hdl_results.csv")
cox_path <- file.path(flt_root, "01cox", "cox_results_table.csv")
uv_path <- file.path(flt_run("uvmr"), "uvmr_results.csv")
mv_path <- file.path(flt_run("mvmr"), "mvmr_conditional_F.csv")
coloc_path <- file.path(flt_root, "03placo", "coloc-results.csv")
read_optional <- function(path) if (file.exists(path)) fread(path) else data.table()
save <- function(plot, name, width, height)
    ggsave(file.path(out, paste0(name, ".png")), plot, width = width, height = height,
           dpi = 300, limitsize = FALSE, bg = "white")

rg <- read_optional(rg_path)
if (nrow(rg)) {
    d <- merge(rg[method == "ldsc"], meta[, .(trait1, item_raw, category)], by = "trait1")
    d[, label := factor(label, levels = registry$label)]
    d[, order := match(trait1, meta$trait1)]
    d[, category := factor(category, levels = unique(meta$category))]
    d[, low_precision := is.finite(rg) & !(h2_precision %in% "both_h2_Z_ge_4")]
    d[, display := fifelse(low_precision, NA_real_, rg)]
    pos <- d[, .(order = mean(range(order))), by = category]
    rng <- max(.5, max(abs(d$display), na.rm = TRUE))
    main <- d[panel %in% c("primary", "control")]
    heatmap <- ggplot(main, aes(order, label, fill = display)) +
        geom_tile(width = .95, height = .90) +
        geom_point(data = main[is.finite(FDR) & FDR < .05 & !low_precision],
                   shape = 21, fill = "black", colour = "white", size = 1.4, stroke = .2) +
        geom_point(data = main[low_precision], shape = 4, colour = "#697480", size = .9) +
        scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                             limits = c(-rng, rng), na.value = "#E1E5EA") +
        scale_x_continuous(breaks = pos$order, labels = pos$category, expand = expansion(add = 1)) +
        labs(x = NULL, y = NULL, fill = expression(r[g]),
             title = "Genome-wide genetic correlations",
             caption = "Dots: panel-wide BH q < 0.05. Crosses: at least one h² Z < 4 or unavailable. Grey: unavailable or low precision.") +
        theme_minimal(base_size = 12) +
        theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1))
    save(heatmap, "Genetics-LDSC", 14, 4.8)
    if (any(d$panel == "extended")) {
        alt <- heatmap %+% d[panel == "extended"]
        alt$layers[[2]]$data <- d[panel == "extended" & is.finite(FDR) & FDR < .05 & !low_precision]
        alt$layers[[3]]$data <- d[panel == "extended" & low_precision]
        save(alt + labs(title = "Extended GWAS panel · exploratory"), "Genetics-LDSC-extended", 15, 10)
    }
}
if (file.exists(cox_path) && file.exists(uv_path)) {
    cx <- flt_cox()
    cx <- merge(as.data.table(cx), meta[, .(name, trait1)], by = "name")
    cx[, trait2 := names(flt_labels)[match(disease, flt_labels)]]
    setnames(cx, c("FDR", "P"), c("Cox_q", "Cox_P"))
    raw <- read_optional(uv_path)
    pick <- function(dir) {
        z <- raw[direction == dir & method_id == "mr_ivw"]
        if (!nrow(z)) return(data.table(trait1 = character(), trait2 = character()))
        if (dir == "forward") {
            dat <- z[, .(trait1 = exposure, trait2 = outcome, MR_b = b, MR_SE = se,
                         MR_P = pval, MR_q = FDR, MR_n = nsnp, MR_F = F_mean)]
        } else {
            dat <- z[, .(trait1 = outcome, trait2 = exposure, reverse_b = b, reverse_SE = se,
                         reverse_P = pval, reverse_q = FDR, reverse_n = nsnp, reverse_F = F_mean)]
        }
        stopifnot(!anyDuplicated(dat[, .(trait1, trait2)]))
        dat
    }
    ev <- merge(cx, pick("forward"), by = c("trait1", "trait2"), all.x = TRUE, sort = FALSE)
    ev <- merge(ev, pick("reverse"), by = c("trait1", "trait2"), all.x = TRUE, sort = FALSE)
    if (nrow(rg)) ev <- merge(ev, rg[method == "ldsc", .(trait1, trait2, rg, rg_se, rg_p,
        LDSC_q = FDR, h2_precision)], by = c("trait1", "trait2"), all.x = TRUE, sort = FALSE)
    mv <- read_optional(mv_path)
    if (nrow(mv)) {
        mv <- mv[food == exposure, .(trait1 = food, trait2 = outcome, MVMR_conditional_F = conditional_F,
                                      MVMR_covariance_assumption = covariance_assumption)]
        ev <- merge(ev, mv, by = c("trait1", "trait2"), all.x = TRUE, sort = FALSE)
    }
    fwrite(ev, file.path(out, "Evidence-by-pair.csv"))
    fwrite(ev[trait1 == "G732" & trait2 == "mdd_2025"],
           file.path(out, "Chicken-MDD-evidence.csv"))
}
sensitivity <- read_optional(file.path(flt_root, "01cox", "cox-sensitivity.csv"))
if (nrow(sensitivity) && file.exists(cox_path)) {
    d <- sensitivity[name == "p20638" & outcome == "depression",
                     .(model, N, events, HR, lower, upper, PH_P)]
    cx <- as.data.table(flt_cox())[name == "p20638" & disease == "Depression"]
    if (nrow(cx) == 1L) d <- rbindlist(list(data.table(model = "FPQ primary", N = cx$N, events = cx$events,
        HR = cx$HR, lower = cx$CI_low, upper = cx$CI_high, PH_P = NA_real_), d))
    fwrite(d, file.path(out, "Chicken-MDD-cox.csv"))
    d <- d[is.finite(HR) & is.finite(lower) & lower > 0 & is.finite(upper)]
    if (nrow(d)) {
        d[, model := factor(paste0(model, "  (", events, " cases)"),
                            levels = rev(paste0(model, "  (", events, " cases)")))]
        plot <- ggplot(d, aes(HR, model)) +
            geom_vline(xintercept = 1, linetype = 2, colour = "#8E9AA4") +
            geom_segment(aes(x = lower, xend = upper, yend = model), colour = "#285D73") +
            geom_point(colour = "#285D73", size = 2.6) +
            scale_x_log10() +
            labs(title = "Chicken liking and incident depression", x = "HR per one-point higher liking (95% CI)", y = NULL,
                 caption = "FPQ = YGB food preferences questionnaire. Common-sample comparisons use the same participant subset.") +
            theme_minimal(base_size = 12)
        save(plot, "Chicken-MDD-cox", 8, 4.1)
    }
}
coloc <- read_optional(coloc_path)
if (nrow(coloc)) {
    d <- coloc[status == "ok" & p12 == 5e-6][order(-PP_H4, PLACO_P)]
    if (nrow(d)) {
        d <- head(d, 24L)
        d[, region := factor(paste0("Chr ", chr, " · ", lead), levels = rev(paste0("Chr ", chr, " · ", lead)))]
        long <- melt(d, id.vars = "region", measure.vars = paste0("PP_H", 0:4),
                     variable.name = "hypothesis", value.name = "posterior")
        plot <- ggplot(long, aes(posterior, region, fill = hypothesis)) +
            geom_col() +
            scale_fill_manual(values = c("#EAEFF3", "#B5C4D5", "#829EB8", "#DA9266", "#2F7E7A"),
                              labels = c("H0 Neither trait", "H1 FLT only", "H2 MDD only",
                                         "H3 Distinct variants", "H4 Shared variant")) +
            labs(title = "Colocalization in PLACO candidate regions", subtitle = "Single-signal model · top 24 by PP.H4",
                 x = "Posterior probability", y = NULL, fill = NULL,
                 caption = "p12 = 5 × 10⁻⁶ per SNP. All windows and prior sensitivities are in the source table.") +
            theme_minimal(base_size = 12) + theme(legend.position = "bottom")
        save(plot, "Genetics-colocalization", 9, max(4, 2.5 + .23 * nrow(d)))
    }
}
if (!any(file.exists(c(rg_path, cox_path, uv_path, coloc_path)))) {
    writeLines("No completed Cox, LDSC, UVMR or coloc input was available.",
               file.path(Sys.getenv("FLT_PAPER_REPORTS"), "genetics-coverage.txt"))
}
