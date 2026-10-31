#!/usr/bin/env Rscript
# mechanism: protein, pathways, network, pathway_plot, metabolite, metabolite_plot, inflammation, inflammation_plot, mri, mri_refit, mri_wmh, sem, sem_plot, sem_sensitivity, summary
# Launch with dose.sh; source this file to register definitions for reuse.
if (!exists("dose_step_code", mode="function")) {
    source(file.path(Sys.getenv("DOSE_SCRIPT_DIR", if (.Platform$OS.type == "windows")
        "D:/scripts/dose/f" else "/mnt/d/scripts/dose/f"), "0.common.R"), local=TRUE)
}

# ===========================================================================
# mechanism / protein
# ===========================================================================
dose_register_step("mechanism", "protein", expression({
# DOSE -> protein -> incident ischemic stroke (IS)
# Proteomics-specific import/QC; shared models are in 0.common.R.

options(stringsAsFactors = FALSE)

DIR_PROJECT <- dose_output_dir()
DIR_UKB <- dose_ukb_dir()
DIR_DATA <- file.path(DIR_PROJECT, "ukb/00data")
DIR_SCORE <- file.path(DIR_PROJECT, "ukb/04score")
DIR_OUT <- file.path(DIR_PROJECT, "ukb/05mechanism/protein/M2")

FILE_FUNCTIONS <- file.path(dose_script_dir(), "0.common.R")
FILE_COHORT <- file.path(DIR_DATA, "dat.base.imp.rds")
FILE_CONFIG <- file.path(DIR_DATA, "analysis_config.rds")
FILE_SCORE <- file.path(DIR_SCORE, "dose_score.rds")
FILE_PROT <- file.path(DIR_UKB, "Rdata/prot.rds")

OUTCOME <- "cvd_stroke_i"
TIME_VAR <- paste0(OUTCOME, ".t2e")
EVENT_VAR <- paste0(OUTCOME, ".Yt2e")
MIN_N <- 1000L
MIN_EVENTS <- 20L
MEDIATION_SIMS <- 100L
SEED <- 12345L

if (!file.exists(FILE_FUNCTIONS)) stop("Missing shared functions: ", FILE_FUNCTIONS)

s4_load_packages()
dir.create(DIR_OUT, recursive = TRUE, showWarnings = FALSE)

input_files <- c(functions = FILE_FUNCTIONS, cohort = FILE_COHORT,
                 config = FILE_CONFIG, dose = FILE_SCORE, protein = FILE_PROT)
if (any(!file.exists(input_files))) {
    stop("Missing input file(s): ",
         paste(input_files[!file.exists(input_files)], collapse = "; "))
}

cfg <- readRDS(FILE_CONFIG)
dat <- s4_validate_ids(as.data.frame(readRDS(FILE_COHORT)), "dat")
dose_checked <- s4_validate_dose(readRDS(FILE_SCORE))
dose <- dose_checked$data
dose_meta <- dose_checked$meta
prot <- s4_validate_ids(as.data.frame(readRDS(FILE_PROT)), "prot")

if (!OUTCOME %in% cfg$Ys ||
    !identical(cfg$primary_development_outcome, OUTCOME) ||
    !identical(cfg$primary_model, "M2")) {
    stop("analysis_config.rds does not define cvd_stroke_i/M2 as primary.")
}
m2_covariates <- s4_get_m2_covariates(cfg)

# Protein values are used as supplied by the UKB-PPP prepared file and are
# Z-standardized directly within each protein-specific complete-case sample.
protein_vars0 <- setdiff(names(prot), "eid")
prepared <- s4_prepare_base(
    dat = dat, dose = dose, mediator_data = prot,
    mediator_vars = protein_vars0, time_var = TIME_VAR,
    event_var = EVENT_VAR, covariates = m2_covariates,
    cohort_label = "proteomic"
)
datp0 <- prepared$data

protein_qc <- data.frame(
    protein = protein_vars0,
    is_numeric = vapply(datp0[protein_vars0],
                        function(x) is.numeric(x) || is.integer(x), logical(1)),
    N_nonmissing = vapply(datp0[protein_vars0],
                          function(x) sum(!is.na(x)), integer(1)),
    percent_missing = 100 * vapply(datp0[protein_vars0],
                                    function(x) mean(is.na(x)), numeric(1)),
    N_unique = vapply(datp0[protein_vars0],
                      function(x) dplyr::n_distinct(x, na.rm = TRUE), integer(1)),
    SD = vapply(datp0[protein_vars0], function(x) {
        if (is.numeric(x) || is.integer(x)) stats::sd(as.numeric(x), na.rm = TRUE)
        else NA_real_
    }, numeric(1))
)
protein_vars <- protein_qc %>% dplyr::filter(
    is_numeric, N_nonmissing >= MIN_N, N_unique > 1L, is.finite(SD), SD > 0
) %>% dplyr::pull(protein)
if (!length(protein_vars)) stop("No protein passed QC.")

qc_sample <- data.frame(
    main_cohort_N = nrow(dat), proteomics_file_N = nrow(prot),
    overlap_N = prepared$overlap_N, complete_M2_N = nrow(datp0),
    complete_M2_events = sum(datp0[[EVENT_VAR]]),
    event_percent = 100 * mean(datp0[[EVENT_VAR]]),
    dose_mean = mean(datp0$dose), dose_SD = stats::sd(datp0$dose),
    dose_min = min(datp0$dose), dose_max = max(datp0$dose),
    proteins_raw_N = length(protein_vars0), proteins_QC_N = length(protein_vars),
    transformation = "direct Z-standardization"
)
s4_write_csv(qc_sample, file.path(DIR_OUT, "00_QC_sample.csv"))
s4_write_csv(protein_qc, file.path(DIR_OUT, "00_QC_protein_variables.csv"))
s4_write_csv(data.frame(covariate = m2_covariates),
             file.path(DIR_OUT, "00_M2_covariates_used.csv"))

message("Proteomic complete M2 cohort: N=", nrow(datp0),
        "; IS events=", sum(datp0[[EVENT_VAR]]),
        "; proteins=", length(protein_vars))
rm(dat, prot, dose, dose_checked, prepared)
invisible(gc())

result <- s4_run_analysis(
    dat = datp0, mediator_vars = protein_vars,
    transform_mediator = s4_transform_direct_z,
    mediator_label = "protein", time_var = TIME_VAR,
    event_var = EVENT_VAR, covariates = m2_covariates,
    min_n = MIN_N, min_events = MIN_EVENTS,
    sims = MEDIATION_SIMS, seed = SEED,
    fdr_family_n = length(protein_vars)
)
result$summary <- cbind(
    data.frame(N_main_cohort = qc_sample$main_cohort_N,
               N_proteomic_complete_M2 = nrow(datp0),
               N_proteins_raw = length(protein_vars0),
               N_proteins_QC = length(protein_vars)),
    result$summary
)

s4_export_results(
    result, DIR_OUT, mediator_label = "protein",
    workbook_name = "DOSE_protein_IS_three_method_mediation_results.xlsx",
    candidate_prefix = "PPI_protein",
    extra_sheets = list("Sample QC" = qc_sample,
                        "Protein QC" = protein_qc,
                        "M2 covariates" = data.frame(covariate = m2_covariates)),
    notes = c(
        "Proteins are used as supplied and standardized directly; no log transformation is applied.",
        "Protein measurements may precede WebQ dietary assessment; interpret pathways as exploratory."
    )
)
s4_write_manifest(input_files, DIR_OUT)
print(result$summary)
message("Completed protein mediation analysis. Results: ", DIR_OUT)
}))

# ===========================================================================
# mechanism / pathways
# ===========================================================================
dose_register_step("mechanism", "pathways", expression({
# DOSE-associated and IS-associated protein candidates:
# pathway enrichment and STRING PPI preparation.
#
# Primary gene set:
#   X-M BH-FDR < 0.05 AND M-Y Cox P < 0.05, as generated by
#   mechanism.R.
# Sensitivity gene set:
#   candidates whose product beta_XM * beta_MY has the same sign as the
#   overall DOSE -> IS Cox log-HR.

options(stringsAsFactors = FALSE)

DIR_PROJECT <- dose_output_dir()
DIR_MED <- file.path(
    DIR_PROJECT, "ukb/05mechanism/protein/M2"
)
DIR_OUT <- file.path(DIR_MED, "PPI_pathway")

FILE_CANDIDATES <- file.path(
    DIR_MED, "04_candidate_protein_XM_FDR05_MY_P05.csv"
)
FILE_TOTAL <- file.path(DIR_MED, "00_overall_DOSE_IS_total_effect.csv")

TOP_N_PER_DATABASE <- 12L
STRING_SPECIES <- 9606L

required_packages <- c(
    "dplyr", "stringr", "ggplot2", "enrichR", "openxlsx"
)
missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
    stop(
        "Install required packages first: ",
        paste(missing_packages, collapse = ", ")
    )
}
suppressPackageStartupMessages({
    library(dplyr)
    library(stringr)
    library(ggplot2)
    # Attaching enrichR initializes the service URL used by its API calls.
    library(enrichR)
})

dir.create(DIR_OUT, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(FILE_CANDIDATES)) {
    stop(
        "Candidate file not found: ", FILE_CANDIDATES,
        "\nRun mechanism.R first."
    )
}
if (!file.exists(FILE_TOTAL)) {
    stop(
        "Overall DOSE-IS result not found: ", FILE_TOTAL,
        "\nRun mechanism.R first."
    )
}

write_csv <- function(x, filename) {
    utils::write.csv(
        x, file.path(DIR_OUT, filename), row.names = FALSE, na = ""
    )
}

write_gene_list <- function(x, filename) {
    x <- unique(as.character(x))
    x <- x[!is.na(x) & nzchar(x)]
    writeLines(x, file.path(DIR_OUT, filename), useBytes = TRUE)
}

assert_columns <- function(x, required, object_name) {
    missing <- setdiff(required, names(x))
    if (length(missing)) {
        stop(
            object_name, " is missing required column(s): ",
            paste(missing, collapse = ", ")
        )
    }
}

# Protein column names in the UKB-PPP prepared file are expected to be gene
# symbols. Preserve the original protein name and use an upper-case symbol for
# enrichment/STRING. A small explicit correction handles R-safe HLA names.
protein_to_gene <- function(x) {
    gene <- toupper(trimws(as.character(x)))
    gene <- sub("^HLA_([A-Z0-9]+)$", "HLA-\\1", gene)
    gene[!nzchar(gene)] <- NA_character_
    gene
}

candidates <- utils::read.csv(
    FILE_CANDIDATES, check.names = FALSE, stringsAsFactors = FALSE
)
if (!"protein" %in% names(candidates) && "mediator" %in% names(candidates)) {
    names(candidates)[names(candidates) == "mediator"] <- "protein"
}
assert_columns(
    candidates,
    c("protein", "beta_XM", "P_XM", "FDR_XM", "beta_MY", "HR_MY",
      "P_MY", "FDR_MY"),
    "Protein candidate table"
)

total_effect <- utils::read.csv(
    FILE_TOTAL, check.names = FALSE, stringsAsFactors = FALSE
)
assert_columns(total_effect, c("HR_DOSE_per_SD", "P"), "Total-effect table")
if (nrow(total_effect) != 1L ||
    !is.finite(total_effect$HR_DOSE_per_SD[1]) ||
    total_effect$HR_DOSE_per_SD[1] <= 0) {
    stop("The total-effect file must contain one finite positive HR.")
}

total_hr <- total_effect$HR_DOSE_per_SD[1]
total_log_hr <- log(total_hr)
total_p <- total_effect$P[1]

candidate_nodes <- candidates %>%
    mutate(
        Gene = protein_to_gene(protein),
        indirect_product_screen = beta_XM * beta_MY,
        direction = case_when(
            beta_XM < 0 & beta_MY > 0 ~ "DOSE lowers risk-associated protein",
            beta_XM > 0 & beta_MY < 0 ~ "DOSE raises protective protein",
            beta_XM > 0 & beta_MY > 0 ~ "DOSE raises risk-associated protein",
            beta_XM < 0 & beta_MY < 0 ~ "DOSE lowers protective protein",
            TRUE ~ "Zero or unclear"
        ),
        same_direction_as_total = case_when(
            !is.finite(indirect_product_screen) ~ NA,
            !is.finite(total_log_hr) | total_log_hr == 0 ~ NA,
            TRUE ~ sign(indirect_product_screen) == sign(total_log_hr)
        ),
        directionally_consistent = same_direction_as_total %in% TRUE,
        overall_DOSE_IS_HR = total_hr,
        overall_DOSE_IS_logHR = total_log_hr,
        overall_DOSE_IS_P = total_p
    ) %>%
    filter(!is.na(Gene), Gene != "") %>%
    arrange(desc(directionally_consistent), FDR_XM, P_MY)

duplicate_gene_qc <- candidate_nodes %>%
    count(Gene, name = "N_protein_columns") %>%
    filter(N_protein_columns > 1L) %>%
    arrange(desc(N_protein_columns), Gene)

# STRING/enrichment needs one row per gene; retain the strongest candidate row
# when duplicate protein columns map to the same gene.
candidate_nodes_unique <- candidate_nodes %>%
    arrange(Gene, FDR_XM, P_MY) %>%
    distinct(Gene, .keep_all = TRUE)
direction_nodes <- candidate_nodes_unique %>%
    filter(directionally_consistent)

genes_overlap <- candidate_nodes_unique$Gene
genes_direction <- direction_nodes$Gene
if (!length(genes_overlap)) stop("No valid gene symbol remained after QC.")

write_csv(candidate_nodes, "Protein_overlap_PPI_nodes_all_rows.csv")
write_csv(candidate_nodes_unique, "Protein_overlap_PPI_nodes_unique_genes.csv")
write_csv(direction_nodes,
          "Protein_directionally_consistent_PPI_nodes.csv")
write_csv(duplicate_gene_qc, "Protein_gene_symbol_duplicate_QC.csv")
write_gene_list(genes_overlap, "Protein_overlap_gene_list.txt")
write_gene_list(genes_overlap, "Protein_overlap_STRING_input.txt")
write_gene_list(
    genes_direction, "Protein_directionally_consistent_gene_list.txt"
)
write_gene_list(
    genes_direction, "Protein_directionally_consistent_STRING_input.txt"
)

make_string_url <- function(genes, species = STRING_SPECIES) {
    if (!length(genes)) return(NA_character_)
    identifiers <- utils::URLencode(
        paste(genes, collapse = "\r"), reserved = TRUE
    )
    paste0(
        "https://string-db.org/cgi/network?species=", species,
        "&identifiers=", identifiers
    )
}

string_urls <- data.frame(
    gene_set = c("Overlap_candidates", "Directionally_consistent"),
    N_genes = c(length(genes_overlap), length(genes_direction)),
    STRING_url = c(
        make_string_url(genes_overlap), make_string_url(genes_direction)
    )
)
write_csv(string_urls, "Protein_STRING_urls.csv")
writeLines(
    string_urls$STRING_url[!is.na(string_urls$STRING_url)],
    file.path(DIR_OUT, "Protein_STRING_urls.txt"), useBytes = TRUE
)

message(
    "Protein sets prepared: overlap candidates=", length(genes_overlap),
    "; directionally consistent=", length(genes_direction),
    "; total DOSE-IS HR=", signif(total_hr, 4),
    " (P=", signif(total_p, 4), ")."
)

# ---------------------------------------------------------------------------
# Enrichr: GO Biological Process / Reactome / KEGG / WikiPathways
# ---------------------------------------------------------------------------

available_dbs <- enrichR::listEnrichrDbs()$libraryName

choose_database <- function(preferred, regex, available = available_dbs) {
    hit <- preferred[preferred %in% available]
    if (length(hit)) return(hit[1])
    hit <- grep(regex, available, value = TRUE)
    if (!length(hit)) {
        stop(
            "No Enrichr database matched: ", regex,
            ". Check enrichR::listEnrichrDbs()."
        )
    }
    # Prefer the latest year when the requested historical database is absent.
    sort(hit, decreasing = TRUE)[1]
}

db_map <- c(
    GO_Biological_Process = choose_database(
        c("GO_Biological_Process_2023", "GO_Biological_Process_2025",
          "GO_Biological_Process_2024"),
        "^GO_Biological_Process_[0-9]{4}$"
    ),
    Reactome = choose_database(
        c("Reactome_2022", "Reactome_2024", "Reactome_2023"),
        "^Reactome_[0-9]{4}$"
    ),
    KEGG_Human = choose_database(
        c("KEGG_2021_Human", "KEGG_2024_Human", "KEGG_2023_Human"),
        "^KEGG_[0-9]{4}_Human$"
    ),
    WikiPathways_Human = choose_database(
        c("WikiPathway_2023_Human", "WikiPathways_2024_Human",
          "WikiPathways_2023_Human"),
        "^WikiPathways?_[0-9]{4}_Human$"
    )
)
dbs <- unname(db_map)

run_enrichr <- function(genes, gene_set) {
    if (!length(genes)) return(data.frame())
    message("Running Enrichr for ", gene_set, ": ", length(genes), " genes.")
    ans <- enrichR::enrichr(genes, dbs)
    bind_rows(lapply(names(ans), function(db) {
        x <- ans[[db]]
        if (!nrow(x)) return(NULL)
        x$Database <- db
        x
    })) %>%
        mutate(
            Gene_set = gene_set,
            Database_group = names(db_map)[match(Database, db_map)],
            Term2 = Term %>%
                str_remove(" \\(GO:\\d+\\)$") %>%
                str_remove(" R-HSA-\\d+$"),
            log10FDR = -log10(pmax(Adjusted.P.value, 1e-300)),
            Overlap_n = suppressWarnings(
                as.numeric(str_extract(Overlap, "^[0-9]+"))
            ),
            Database_group = factor(
                Database_group,
                levels = names(db_map)
            )
        ) %>%
        arrange(Adjusted.P.value)
}

enrichment_overlap <- run_enrichr(genes_overlap, "Overlap_candidates")
enrichment_direction <- run_enrichr(
    genes_direction, "Directionally_consistent"
)
enrichment_all <- bind_rows(enrichment_overlap, enrichment_direction)
write_csv(enrichment_all, "Protein_pathway_enrichment_all.csv")
write_csv(
    enrichment_all %>% filter(Adjusted.P.value < 0.05),
    "Protein_pathway_enrichment_FDR05.csv"
)

plot_enrichment <- function(dat, gene_set, n = TOP_N_PER_DATABASE) {
    d <- dat %>%
        filter(Gene_set == gene_set, is.finite(log10FDR)) %>%
        group_by(Database_group) %>%
        slice_min(Adjusted.P.value, n = n, with_ties = FALSE) %>%
        ungroup() %>%
        mutate(
            Term_plot = str_wrap(Term2, width = 55),
            panel_term = paste(Database_group, Term_plot, sep = "___")
        )
    if (!nrow(d)) return(NULL)
    d$panel_term <- factor(
        d$panel_term,
        levels = rev(unique(d$panel_term[order(d$Database_group,
                                               d$Adjusted.P.value)]))
    )
    ggplot(d, aes(log10FDR, panel_term, size = Overlap_n,
                  color = log10FDR)) +
        geom_point(alpha = 0.9) +
        facet_grid(Database_group ~ ., scales = "free_y", space = "free_y") +
        scale_y_discrete(labels = function(x) sub("^.*___", "", x)) +
        scale_color_gradient(low = "#80B1D3", high = "#D95F5F") +
        labs(
            title = paste0(gene_set, " pathway enrichment"),
            subtitle = paste0(
                "Top ", n, " terms per database; dashed line: FDR=0.05"
            ),
            x = expression(-log[10](FDR)), y = NULL,
            size = "Overlap", color = expression(-log[10](FDR))
        ) +
        geom_vline(
            xintercept = -log10(0.05), linetype = 2, color = "grey50"
        ) +
        theme_bw(base_size = 11) +
        theme(
            plot.title = element_text(face = "bold", hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            axis.text.y = element_text(size = 8),
            strip.text.y = element_text(face = "bold")
        )
}

p_overlap <- plot_enrichment(enrichment_all, "Overlap_candidates")
if (!is.null(p_overlap)) {
    bounded_ggsave(
        file.path(DIR_OUT, "Fig_protein_overlap_pathway_enrichment.pdf"),
        p_overlap, width = 10.5, height = 12, device = grDevices::cairo_pdf
    )
}
p_direction <- plot_enrichment(
    enrichment_all, "Directionally_consistent"
)
if (!is.null(p_direction)) {
    bounded_ggsave(
        file.path(
            DIR_OUT,
            "Fig_protein_directionally_consistent_pathway_enrichment.pdf"
        ),
        p_direction, width = 10.5, height = 12,
        device = grDevices::cairo_pdf
    )
}

# ---------------------------------------------------------------------------
# Optional: import enrichment table downloaded from STRING and plot its top
# terms. Put the downloaded TSV in DIR_OUT under one of the names below, then
# rerun this script. No manually chosen old terms are imposed on the new data.
# ---------------------------------------------------------------------------

string_enrichment_candidates <- file.path(
    DIR_OUT,
    c("PPI.enrichment.all.tsv", "PPI.nrichment.all.tsv",
      "enrichment.enrichment.tsv", "enrichment.tsv")
)
string_enrichment_file <- string_enrichment_candidates[
    file.exists(string_enrichment_candidates)
][1]

standardize_string_enrichment <- function(path) {
    x <- utils::read.delim(
        path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE
    )
    aliases <- list(
        category = c("#category", "category"),
        term_id = c("term ID", "term_id", "term"),
        term_description = c("term description", "description"),
        observed_gene_count = c("observed gene count", "number_of_genes"),
        background_gene_count = c("background gene count",
                                  "number_of_genes_in_background"),
        FDR = c("false discovery rate", "fdr"),
        matching_proteins = c(
            "matching proteins in your network (labels)",
            "preferredNames", "inputGenes"
        )
    )
    for (new_name in names(aliases)) {
        old_name <- aliases[[new_name]][aliases[[new_name]] %in% names(x)][1]
        if (!is.na(old_name)) names(x)[names(x) == old_name] <- new_name
    }
    assert_columns(
        x,
        c("category", "term_id", "term_description",
          "observed_gene_count", "FDR"),
        "STRING enrichment table"
    )
    x %>%
        mutate(
            observed_gene_count = as.numeric(observed_gene_count),
            FDR = as.numeric(FDR),
            log10FDR = -log10(pmax(FDR, 1e-300)),
            Facet = case_when(
                category == "GO Process" ~ "GO Biological Process",
                category == "GO Function" ~ "GO Molecular Function",
                category == "GO Component" ~ "GO Cellular Component",
                category == "STRING clusters" ~ "STRING Cluster",
                category == "Reactome" ~ "Reactome",
                category == "KEGG" ~ "KEGG",
                str_detect(category, regex("Wiki", ignore_case = TRUE)) ~
                    "WikiPathways",
                TRUE ~ category
            )
        ) %>%
        filter(is.finite(FDR), is.finite(log10FDR)) %>%
        arrange(FDR)
}

if (!is.na(string_enrichment_file)) {
    string_enrichment <- standardize_string_enrichment(
        string_enrichment_file
    )
    write_csv(string_enrichment, "STRING_enrichment_all_standardized.csv")

    string_plot_data <- string_enrichment %>%
        filter(FDR < 0.05) %>%
        group_by(Facet) %>%
        slice_min(FDR, n = TOP_N_PER_DATABASE, with_ties = FALSE) %>%
        ungroup() %>%
        mutate(
            Term2 = str_wrap(term_description, width = 55),
            panel_term = paste(Facet, Term2, sep = "___")
        )
    if (nrow(string_plot_data)) {
        string_plot_data$panel_term <- factor(
            string_plot_data$panel_term,
            levels = rev(unique(string_plot_data$panel_term[
                order(string_plot_data$Facet, string_plot_data$FDR)
            ]))
        )
        p_string <- ggplot(
            string_plot_data,
            aes(log10FDR, panel_term, size = observed_gene_count,
                color = Facet)
        ) +
            geom_segment(
                aes(x = 0, xend = log10FDR, yend = panel_term),
                color = "grey82", linewidth = 2.5, lineend = "round",
                show.legend = FALSE
            ) +
            geom_point(alpha = 0.95) +
            facet_grid(Facet ~ ., scales = "free_y", space = "free_y") +
            scale_y_discrete(labels = function(x) sub("^.*___", "", x)) +
            scale_size_continuous(range = c(2.8, 6.5)) +
            labs(
                x = expression(-log[10](FDR)), y = NULL,
                color = NULL, size = "Overlap",
                title = "STRING functional enrichment"
            ) +
            theme_classic(base_size = 11) +
            theme(
                legend.position = "right",
                strip.background = element_rect(fill = "grey92", color = NA),
                strip.text.y = element_text(face = "bold"),
                axis.text.y = element_text(size = 8)
            )
        bounded_ggsave(
            file.path(DIR_OUT, "Fig_STRING_enrichment.pdf"),
            p_string, width = 9.5,
            height = max(7, 2 + 0.22 * nrow(string_plot_data)),
            device = grDevices::cairo_pdf, limitsize = FALSE
        )
    }
} else {
    writeLines(
        c(
            "STRING enrichment table was not found.",
            "1. Open Protein_STRING_urls.csv or upload Protein_overlap_STRING_input.txt to STRING.",
            "2. Select Homo sapiens (taxon 9606).",
            "3. Export the STRING enrichment table as TSV.",
            "4. Save it here as PPI.enrichment.all.tsv and rerun this script."
        ),
        file.path(DIR_OUT, "README_STRING_manual_step.txt"),
        useBytes = TRUE
    )
    message(
        "STRING enrichment TSV not found. Input files and URLs were created; ",
        "see README_STRING_manual_step.txt."
    )
}

summary_counts <- data.frame(
    N_candidate_protein_rows = nrow(candidate_nodes),
    N_unique_candidate_genes = length(genes_overlap),
    N_directionally_consistent_genes = length(genes_direction),
    N_duplicate_gene_symbols = nrow(duplicate_gene_qc),
    overall_DOSE_IS_HR = total_hr,
    overall_DOSE_IS_P = total_p,
    N_enrichment_rows = nrow(enrichment_all),
    N_enrichment_FDR05 = sum(
        enrichment_all$Adjusted.P.value < 0.05, na.rm = TRUE
    )
)
write_csv(summary_counts, "00_PPI_pathway_summary.csv")
write_csv(
    data.frame(Database_group = names(db_map), Enrichr_database = db_map),
    "00_Enrichr_databases_used.csv"
)

workbook_sheets <- list(
    "Summary" = summary_counts,
    "Candidate nodes" = candidate_nodes_unique,
    "Directional nodes" = direction_nodes,
    "Duplicate gene QC" = if (nrow(duplicate_gene_qc)) duplicate_gene_qc else
        data.frame(note = "No duplicated gene symbol"),
    "Enrichr databases" = data.frame(
        Database_group = names(db_map), Enrichr_database = db_map
    ),
    "Pathway enrichment" = if (nrow(enrichment_all)) enrichment_all else
        data.frame(note = "No enrichment result"),
    "Pathway FDR05" = if (any(enrichment_all$Adjusted.P.value < 0.05,
                               na.rm = TRUE)) {
        enrichment_all %>% filter(Adjusted.P.value < 0.05)
    } else data.frame(note = "No pathway with FDR < 0.05")
)
openxlsx::write.xlsx(
    workbook_sheets,
    file.path(DIR_OUT, "DOSE_protein_PPI_pathway_enrichment.xlsx"),
    overwrite = TRUE, asTable = TRUE
)

message("Completed protein PPI/pathway analysis. Results: ", DIR_OUT)
print(summary_counts)
}))

# ===========================================================================
# mechanism / network
# ===========================================================================
dose_register_step("mechanism", "network", expression({
# STRING v12.0 network for the full DOSE/IS association intersection.
# Run: D:/software/win/R-4.5.1/bin/Rscript.exe f/mechanism.R
# Cached API responses make subsequent figure generation reproducible.
options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
score_cutoff <- if (length(args)) as.numeric(args[1]) else 0.4
stopifnot(length(score_cutoff) == 1L, is.finite(score_cutoff), score_cutoff >= 0, score_cutoff < 1)
required_score <- as.integer(round(score_cutoff * 1000)) + 1L
score_tag <- sprintf("%03d", round(score_cutoff * 100))
suppressPackageStartupMessages({library(httr); library(jsonlite); library(igraph); library(ggplot2)})
root <- dose_output_dir()
input <- file.path(root, "ukb/05mechanism/protein/M2/04_candidate_protein_XM_FDR05_MY_P05.csv")
out <- file.path(root, "ukb/05mechanism/protein/M2/PPI_pathway/STRING_network")
if (score_cutoff != 0.4) out <- file.path(out, paste0("score_gt_", score_tag))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
api <- "https://version-12-0.string-db.org/api/json/"
cand <- read.csv(input, check.names = FALSE)
stopifnot(nrow(cand) > 0L, !anyDuplicated(cand$protein), all(cand$FDR_XM < .05), all(cand$P_MY < .05))
cand$Gene <- sub("^HLA_([A-Z0-9]+)$", "HLA-\\1", toupper(trimws(cand$protein)))
genes <- sort(unique(cand$Gene))

fetch <- function(method, args, filename) {
    path <- file.path(out, filename)
    params <- c(list(species = 9606L, caller_identity = "DOSE_IS_protein_network"), args)
    manifest <- list(url = paste0(api, method), parameters = params)
    request_file <- paste0(path, ".request.json")
    if (file.exists(path)) {
        previous <- fromJSON(request_file, simplifyVector = FALSE)
        if (!identical(toJSON(previous, auto_unbox=TRUE), toJSON(manifest, auto_unbox=TRUE))) stop("Cached request changed: ", path)
    } else {
        Sys.sleep(1)
        response <- RETRY("POST", paste0(api, method), body = params, encode = "form", timeout(60), times = 3, pause_min = 2)
        stop_for_status(response)
        txt <- content(response, as = "text", encoding = "UTF-8")
        parsed <- fromJSON(txt)
        if (!is.data.frame(parsed) || !nrow(parsed)) stop("Empty or invalid API result: ", method)
        writeLines(txt, path, useBytes = TRUE)
        write_json(manifest, request_file, auto_unbox = TRUE, pretty = TRUE)
    }
    fromJSON(path)
}
mapping <- fetch("get_string_ids", list(identifiers = paste(genes, collapse = "\r"), echo_query = 1), "01_STRING_mapping_raw.json")
mapping$Gene <- genes[mapping$queryIndex + 1L]
stopifnot(!anyDuplicated(mapping$Gene), all(mapping$ncbiTaxonId == 9606))
write.csv(mapping, file.path(out, "01_STRING_mapping.csv"), row.names = FALSE)
unmapped <- setdiff(genes, mapping$Gene)
write.csv(data.frame(Gene = unmapped), file.path(out, "01_unmapped_genes.csv"), row.names = FALSE)
if (anyDuplicated(mapping$stringId)) stop("Multiple candidate genes map to one STRING protein; inspect mapping before proceeding.")
message("Mapped: ", nrow(mapping), "/", length(genes), "; aliases: ", paste(mapping$Gene[mapping$Gene != mapping$preferredName], collapse=", "))
identifiers <- paste(sort(mapping$stringId), collapse = "\r")
edges <- fetch("network", list(identifiers = identifiers, required_score = required_score, network_type = "functional", add_nodes = 0L), "02_STRING_network_raw.json")
edges <- edges[edges$score > score_cutoff & edges$stringId_A != edges$stringId_B, ]
edges$from <- mapping$Gene[match(edges$stringId_A, mapping$stringId)]
edges$to <- mapping$Gene[match(edges$stringId_B, mapping$stringId)]
stopifnot(!anyNA(edges$from), !anyNA(edges$to))
key <- apply(edges[c("from", "to")], 1, function(x) paste(sort(x), collapse="|"))
edges <- edges[!duplicated(key), ]
edges <- edges[order(edges$score, edges$from, edges$to), ]
write.csv(edges, file.path(out, paste0("02_STRING_edges_score_gt_", score_tag, ".csv")), row.names = FALSE)
ppi <- fetch("ppi_enrichment", list(identifiers = identifiers, required_score = required_score, network_type = "functional"), "03_STRING_ppi_enrichment_raw.json")
write.csv(ppi, file.path(out, "03_STRING_ppi_enrichment.csv"), row.names = FALSE)

g <- graph_from_data_frame(edges[c("from", "to", "score")], directed = FALSE, vertices = data.frame(name=genes))
stopifnot(vcount(g) == length(genes), is_simple(g), ecount(g) == nrow(edges))
nodes <- cand[match(genes, cand$Gene), ]
nodes$STRING_ID <- mapping$stringId[match(nodes$Gene, mapping$Gene)]
nodes$STRING_preferred_name <- mapping$preferredName[match(nodes$Gene, mapping$Gene)]
nodes$mapped <- !is.na(nodes$STRING_ID)
nodes$degree <- degree(g)[nodes$Gene]
nodes$strength <- strength(g, weights=E(g)$score)[nodes$Gene]
nodes$betweenness_unweighted_normalized <- betweenness(g, directed=FALSE, normalized=TRUE, weights=NA)[nodes$Gene]
nodes$harmonic_centrality_unweighted_normalized <- harmonic_centrality(g, normalized=TRUE, weights=NA)[nodes$Gene]
nodes$component <- components(g)$membership[nodes$Gene]
write.csv(nodes[order(-nodes$degree, -nodes$strength, nodes$Gene), ], file.path(out, "04_nodes_and_centrality.csv"), row.names=FALSE)
write.csv(nodes[nodes$degree == 0, ], file.path(out, "05_isolated_nodes.csv"), row.names=FALSE)
write_graph(g, file.path(out, "DOSE_IS_STRING_network.graphml"), format="graphml")

# Alphabetical circular layout, all input proteins retained, including isolates.
# A fixed ordering avoids suggesting unsupported biological clusters.
nodes <- nodes[order(nodes$Gene), ]
theta <- pi/2 - 2*pi*(seq_len(nrow(nodes))-1)/nrow(nodes)
nodes$x <- cos(theta); nodes$y <- sin(theta)
nodes$label_x <- 1.055*nodes$x; nodes$label_y <- 1.055*nodes$y
nodes$angle <- theta*180/pi
nodes$hjust <- ifelse(nodes$x < -1e-10, 1, 0)
nodes$angle[nodes$x < -1e-10] <- nodes$angle[nodes$x < -1e-10] + 180
nodes$angle <- (nodes$angle + 180) %% 360 - 180
nodes$display_area <- nodes$degree + 5
# Cubic Bezier curves; control points radially contracted from endpoints.
t <- seq(0, 1, length.out = 60)
curves <- do.call(rbind, lapply(seq_len(nrow(edges)), function(i) {
    a <- nodes[match(edges$from[i], nodes$Gene), ]; b <- nodes[match(edges$to[i], nodes$Gene), ]
    data.frame(edge=i, score=edges$score[i],
      x=(1-t)^3*a$x+3*(1-t)^2*t*.25*a$x+3*(1-t)*t^2*.25*b$x+t^3*b$x,
      y=(1-t)^3*a$y+3*(1-t)^2*t*.25*a$y+3*(1-t)*t^2*.25*b$y+t^3*b$y)
}))
breaks <- unique(c(0, 10, 20, 40, max(nodes$degree)))
breaks <- breaks[breaks <= max(nodes$degree)]
p <- ggplot() +
    geom_path(data=curves, aes(x,y,group=edge,color=score), linewidth=.28, alpha=.43) +
    geom_point(data=nodes, aes(x,y,size=display_area), shape=21, fill="#CCBEDD", color="white", stroke=.3, alpha=.92) +
    geom_text(data=nodes, aes(label_x,label_y,label=Gene,angle=angle,hjust=hjust), size=3.3, family="Arial", color="#202020") +
    scale_color_gradientn(colors=c("#FFF1B5", "#FEC66B", "#EF763C", "#CA200D"), limits=c(.4,1), breaks=c(.4,.6,.8,1), name="STRING confidence") +
    scale_size_area(max_size=12, limits=c(0,max(nodes$display_area)), breaks=breaks+5, labels=breaks, name="Degree\n(number of partners)") +
    coord_equal(xlim=c(-1.32,1.32),ylim=c(-1.29,1.29),clip="off") +
    labs(title="DOSE-associated proteins and ischemic stroke", subtitle=sprintf("%d input proteins  |  %d interactions  |  STRING v12.0, score > %g", length(genes), nrow(edges), score_cutoff),
      caption=sprintf("All %d candidates shown; %d have no retained connection. No additional proteins.\nEdges represent STRING functional associations, including indirect evidence; they do not establish causality.", length(genes), sum(nodes$degree==0))) +
    theme_void(base_family="Arial", base_size=12) +
    theme(plot.title=element_text(size=20,face="bold",hjust=0), plot.subtitle=element_text(size=12,color="#555555",margin=margin(t=9,b=15)),
      plot.caption=element_text(size=10,color="#555555",hjust=0,margin=margin(t=12)),
      legend.position="right", legend.box="vertical", legend.spacing.y=grid::unit(1.2,"cm"),
      legend.title=element_text(size=12), legend.text=element_text(size=11),
      plot.margin=margin(25,30,22,25)) +
    guides(color=guide_colorbar(order=1,barheight=grid::unit(3.5,"cm"),barwidth=grid::unit(.5,"cm")),
      size=guide_legend(order=2,override.aes=list(fill="#CCBEDD",color="white",alpha=1)))
base <- file.path(out, "Fig_PPI_circular")
bounded_ggsave(paste0(base,".pdf"),p,width=18,height=16,device=grDevices::cairo_pdf,bg="white")
bounded_ggsave(paste0(base,".svg"),p,width=18,height=16,device=svglite::svglite,bg="white")
bounded_ggsave(paste0(base,".png"),p,width=18,height=16,dpi=300,bg="white")
bounded_ggsave(file.path(out,"preview.png"),p,width=18,height=16,dpi=100,bg="white")
summary <- data.frame(input_proteins=length(genes),mapped=nrow(mapping),unmapped=length(unmapped),connected=sum(nodes$degree>0),isolated=sum(nodes$degree==0),edges=ecount(g),components=components(g)$no,largest_component=max(components(g)$csize),average_degree=mean(nodes$degree),density=edge_density(g))
write.csv(summary,file.path(out,"00_network_summary.csv"),row.names=FALSE)
writeLines(c(
 "DOSE-IS candidate-protein STRING network", "",
 "Input: current candidate proteins with DOSE-protein BH-FDR<0.05 and protein-IS Cox P<0.05. No mediation-significance filter.",
 "Reference styling/method: https://www.nature.com/articles/s41538-026-01048-3 ; Figure 4E.",
 sprintf("STRING v12.0; Homo sapiens (9606); functional network; combined score strictly >%g (API required_score=%d); add_nodes=0.", score_cutoff, required_score),
 "All input nodes are displayed, including isolates. Node labels use submitted gene symbols; canonical STRING names and IDs are in the mapping file.",
 "Node area scales with degree + 5 to keep isolated nodes visible; legend labels report actual degree. Edge color is combined confidence; all retained edges are shown.",
 "Degree is the number of retained partners among the input proteins. Betweenness and harmonic centrality are unweighted, normalized and computed on the full network.",
 "Harmonic centrality is used because the network can be disconnected. Degree/hub status is topological, not proof of causality.",
 "STRING functional associations include direct and indirect evidence, including text mining; this is not an experimentally verified physical-interaction-only network.",
 "STRING PPI enrichment uses its default proteome null, NOT the 2,910-protein measurement panel. Network connectivity enrichment is not pathway enrichment or mediation testing.",
 paste("Generated:",format(Sys.time(),tz="UTC"),"UTC"),
 paste("Input MD5:",tools::md5sum(input)),
 "Cached raw JSON and request manifests are retained. Delete the relevant cache intentionally if the input/request changes."
), file.path(out,"README_methods.txt"))
capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
print(summary); print(ppi); print(head(nodes[order(-nodes$degree),c("Gene","degree","strength")],10))
message("Completed: ",out)
}))

# ===========================================================================
# mechanism / pathway_plot
# ===========================================================================
dose_register_step("mechanism", "pathway_plot", expression({
# Restyle the previously displayed default-background enrichment results.
# This does not alter the separate measured-background analysis or its entry point.
suppressPackageStartupMessages({library(dplyr);library(stringr);library(ggplot2);library(patchwork)})
root <- file.path(dose_output_dir(), "ukb/05mechanism/protein/M2/PPI_pathway/Top_significant_pathways")
out <- file.path(root,"Bubble_KEGG_Reactome_GO_BP")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
groups <- c("KEGG_Human","Reactome","GO_Biological_Process")
titles <- c(KEGG_Human="A  KEGG",Reactome="B  Reactome",GO_Biological_Process="C  GO Biological Process")
colors <- c(KEGG_Human="#A58BBC",Reactome="#70A7AE",GO_Biological_Process="#8BB79F")
all <- read.csv(file.path(dirname(root),"Protein_pathway_enrichment_all.csv"),check.names=FALSE) %>%
    filter(Gene_set == "Overlap_candidates") %>%
    mutate(Count=Overlap_n, label=Term2, significance=log10FDR)
write.csv(all,file.path(root,"01_all_five_databases.csv"),row.names=FALSE)
d <- all %>% filter(Database_group %in% groups,Adjusted.P.value<.05) %>%
 group_by(Database_group) %>% arrange(Adjusted.P.value,P.value,desc(Count),Term,.by_group=TRUE) %>% slice_head(n=6) %>% ungroup()
stopifnot(nrow(d)>0, all(d$Adjusted.P.value<.05))
write.csv(d,file.path(out,"Displayed_terms.csv"),row.names=FALSE)
colors <- c(KEGG_Human="#4786B1",Reactome="#51A46B",GO_Biological_Process="#ED943A")
d <- d %>% mutate(Category=factor(Database_group,levels=groups)) %>%
 arrange(Category,Adjusted.P.value,P.value,desc(Count),Term)
d$term_id <- factor(d$Term,levels=rev(d$Term))
labels <- str_to_sentence(str_to_lower(d$label))
for (token in c("MAPK","ERK1","ERK2","IGF","IGFBPs","DAP12","TNFs")) {
 labels <- str_replace_all(labels,regex(paste0("\\b",token,"\\b"),ignore_case=TRUE),token)
}
labels <- str_wrap(labels,62)
category_labels <- c(KEGG_Human="KEGG",Reactome="Reactome",GO_Biological_Process="GO\nBiological Process")
p <- ggplot(d,aes(significance,term_id))+
 geom_segment(aes(x=0,xend=significance,yend=term_id),linewidth=3.1,color="#C3C5C7",lineend="butt")+
 geom_point(aes(size=Count,color=Category),alpha=1)+
 facet_grid(Category~.,scales="free_y",space="free_y",switch="y",labeller=labeller(Category=category_labels))+
 scale_y_discrete(labels=setNames(labels,d$Term),expand=expansion(add=.65))+
 scale_x_continuous(limits=c(0,max(13.4,max(d$significance)*1.05)),expand=expansion(mult=c(0,0)))+
 scale_color_manual(values=setNames(unname(colors),groups),guide="none")+
 scale_size_area(max_size=5.8,limits=c(0,max(d$Count)),breaks=c(5,10,20,30,39),name="Overlap")+
 labs(x=expression(-log[10]~italic(FDR)),y=NULL)+
 theme_classic(base_family="Arial",base_size=10)+
 theme(axis.text.y=element_text(size=9.5,color="black",lineheight=1.05,margin=margin(r=3)),
 axis.text.x=element_text(size=9,color="black"),axis.title.x=element_text(size=11,margin=margin(t=7)),
 axis.ticks.y=element_blank(),axis.line=element_line(linewidth=.5),
 strip.placement="outside",strip.background=element_rect(fill="#E5E5E5",color=NA),
 strip.text.y.left=element_text(angle=0,size=10.5,margin=margin(l=7,r=7)),
 panel.spacing.y=grid::unit(.17,"cm"),legend.position="right",legend.justification="top",
 legend.title=element_text(size=9),legend.text=element_text(size=8.5),legend.key.height=grid::unit(.42,"cm"),
 legend.margin=margin(0,0,0,2),plot.margin=margin(12,9,9,7))+
 guides(size=guide_legend(override.aes=list(color="black",alpha=1)))
name <- file.path(out,"Fig_pathway_bubble_KEGG_Reactome_GO_BP")
bounded_ggsave(paste0(name,".pdf"),p,width=10,height=6.5,device=cairo_pdf,bg="white")
bounded_ggsave(paste0(name,".svg"),p,width=10,height=6.5,device=svglite::svglite,bg="white")
bounded_ggsave(paste0(name,".png"),p,width=10,height=6.5,dpi=600,bg="white")
bounded_ggsave(file.path(out,"preview.png"),p,width=10,height=6.5,dpi=150,bg="white")
writeLines(c("Current overlap proteins; Enrichr default background.","Up to six significant terms each for KEGG, Reactome and GO BP. FDR<0.05, ranked by FDR and raw P, then decreasing overlap count and term.","Reference-style grouped bubble plot: gray left category strips, gray stems from zero, colored bubbles; x=-log10(FDR), area=overlap protein count. All categories share one axis and size scale.","This figure is NOT based on the separate Measured_background_2910 results."),file.path(out,"README.txt"))
print(d[,c("Database_group","Term","Count","Adjusted.P.value")],n=20)
}))

# ===========================================================================
# mechanism / metabolite
# ===========================================================================
dose_register_step("mechanism", "metabolite", expression({
# DOSE -> metabolite -> incident ischemic stroke (IS)
# Metabolomics-specific import/QC; shared models are in 0.common.R.

options(stringsAsFactors = FALSE)

DIR_PROJECT <- dose_output_dir()
DIR_UKB <- dose_ukb_dir()
DIR_DATA <- file.path(DIR_PROJECT, "ukb/00data")
DIR_SCORE <- file.path(DIR_PROJECT, "ukb/04score")
DIR_OUT <- file.path(DIR_PROJECT, "ukb/05mechanism/metabolite/M2")

FILE_FUNCTIONS <- file.path(dose_script_dir(), "0.common.R")
FILE_COHORT <- file.path(DIR_DATA, "dat.base.imp.rds")
FILE_CONFIG <- file.path(DIR_DATA, "analysis_config.rds")
FILE_SCORE <- file.path(DIR_SCORE, "dose_score.rds")
FILE_MET <- file.path(DIR_UKB, "Rdata/met.rds")
FILE_MET_LIST <- file.path(DIR_UKB, "common/met.lst")

OUTCOME <- "cvd_stroke_i"
TIME_VAR <- paste0(OUTCOME, ".t2e")
EVENT_VAR <- paste0(OUTCOME, ".Yt2e")
MIN_N <- 1000L
MIN_EVENTS <- 20L
MEDIATION_SIMS <- 100L
SEED <- 12345L

if (!file.exists(FILE_FUNCTIONS)) stop("Missing shared functions: ", FILE_FUNCTIONS)

s4_load_packages(extra = "data.table")
dir.create(DIR_OUT, recursive = TRUE, showWarnings = FALSE)

input_files <- c(functions = FILE_FUNCTIONS, cohort = FILE_COHORT,
                 config = FILE_CONFIG, dose = FILE_SCORE,
                 metabolite = FILE_MET, metabolite_list = FILE_MET_LIST)
if (any(!file.exists(input_files))) {
    stop("Missing input file(s): ",
         paste(input_files[!file.exists(input_files)], collapse = "; "))
}

cfg <- readRDS(FILE_CONFIG)
dat <- s4_validate_ids(as.data.frame(readRDS(FILE_COHORT)), "dat")
dose_checked <- s4_validate_dose(readRDS(FILE_SCORE))
dose <- dose_checked$data
dose_meta <- dose_checked$meta
met <- s4_validate_ids(as.data.frame(readRDS(FILE_MET)), "met")

if (!OUTCOME %in% cfg$Ys ||
    !identical(cfg$primary_development_outcome, OUTCOME) ||
    !identical(cfg$primary_model, "M2")) {
    stop("analysis_config.rds does not define cvd_stroke_i/M2 as primary.")
}
m2_covariates <- s4_get_m2_covariates(cfg)

# met.rds has already passed clean_biom() missingness/variance QC and
# imputation. clean_biom() does not transform concentrations. Analyze only
# direct absolute-concentration metabolites from met.lst; exclude expressions
# containing '/', then apply log(x+1) followed by Z standardization.
met_list <- data.table::fread(
    FILE_MET_LIST, sep = "\t", select = 1:2, header = FALSE,
    fill = TRUE, data.table = FALSE
)
names(met_list) <- c("field_id", "var_name")
met_list$field_id <- as.character(met_list$field_id)
met_list$var_name <- as.character(met_list$var_name)
direct_metabolites <- unique(met_list$var_name[
    !grepl("/", met_list$field_id, fixed = TRUE) &
        !is.na(met_list$var_name) & nzchar(met_list$var_name)
])
metabolite_vars0 <- intersect(direct_metabolites, setdiff(names(met), "eid"))
if (!length(metabolite_vars0)) {
    stop("No direct metabolite from common/met.lst was found in met.rds.")
}
metabolite_inventory <- data.frame(
    metabolite = direct_metabolites,
    present_in_met_rds = direct_metabolites %in% names(met)
)

prepared <- s4_prepare_base(
    dat = dat, dose = dose, mediator_data = met,
    mediator_vars = metabolite_vars0, time_var = TIME_VAR,
    event_var = EVENT_VAR, covariates = m2_covariates,
    cohort_label = "metabolomic"
)
datp0 <- prepared$data

metabolite_qc <- data.frame(
    metabolite = metabolite_vars0,
    is_numeric = vapply(datp0[metabolite_vars0],
                        function(x) is.numeric(x) || is.integer(x), logical(1)),
    N_nonmissing = vapply(datp0[metabolite_vars0],
                          function(x) sum(!is.na(x)), integer(1)),
    percent_missing = 100 * vapply(datp0[metabolite_vars0],
                                    function(x) mean(is.na(x)), numeric(1)),
    N_unique = vapply(datp0[metabolite_vars0],
                      function(x) dplyr::n_distinct(x, na.rm = TRUE), integer(1)),
    N_negative = vapply(datp0[metabolite_vars0], function(x) {
        if (is.numeric(x) || is.integer(x)) sum(as.numeric(x) < 0, na.rm = TRUE)
        else NA_integer_
    }, integer(1)),
    raw_min = vapply(datp0[metabolite_vars0], function(x) {
        if (is.numeric(x) || is.integer(x)) min(as.numeric(x), na.rm = TRUE)
        else NA_real_
    }, numeric(1)),
    SD = vapply(datp0[metabolite_vars0], function(x) {
        if (is.numeric(x) || is.integer(x)) stats::sd(as.numeric(x), na.rm = TRUE)
        else NA_real_
    }, numeric(1))
)
metabolite_vars <- metabolite_qc %>% dplyr::filter(
    is_numeric, N_nonmissing >= MIN_N, N_unique > 1L,
    N_negative == 0L, is.finite(SD), SD > 0
) %>% dplyr::pull(metabolite)
if (!length(metabolite_vars)) stop("No metabolite passed QC.")

qc_sample <- data.frame(
    main_cohort_N = nrow(dat), metabolomics_file_N = nrow(met),
    overlap_N = prepared$overlap_N, complete_M2_N = nrow(datp0),
    complete_M2_events = sum(datp0[[EVENT_VAR]]),
    event_percent = 100 * mean(datp0[[EVENT_VAR]]),
    dose_mean = mean(datp0$dose), dose_SD = stats::sd(datp0$dose),
    dose_min = min(datp0$dose), dose_max = max(datp0$dose),
    direct_metabolites_listed_N = length(direct_metabolites),
    direct_metabolites_present_N = length(metabolite_vars0),
    metabolites_QC_N = length(metabolite_vars),
    transformation = "log(x+1), then Z-standardization"
)
s4_write_csv(qc_sample, file.path(DIR_OUT, "00_QC_sample.csv"))
s4_write_csv(metabolite_qc,
             file.path(DIR_OUT, "00_QC_metabolite_variables.csv"))
s4_write_csv(metabolite_inventory,
             file.path(DIR_OUT, "00_metabolite_inventory_met_list.csv"))
s4_write_csv(data.frame(covariate = m2_covariates),
             file.path(DIR_OUT, "00_M2_covariates_used.csv"))

message("Metabolomic complete M2 cohort: N=", nrow(datp0),
        "; IS events=", sum(datp0[[EVENT_VAR]]),
        "; metabolites=", length(metabolite_vars))
rm(dat, met, dose, dose_checked, prepared)
invisible(gc())

result <- s4_run_analysis(
    dat = datp0, mediator_vars = metabolite_vars,
    transform_mediator = s4_transform_log1p_z,
    mediator_label = "metabolite", time_var = TIME_VAR,
    event_var = EVENT_VAR, covariates = m2_covariates,
    min_n = MIN_N, min_events = MIN_EVENTS,
    sims = MEDIATION_SIMS, seed = SEED,
    fdr_family_n = length(metabolite_vars)
)
result$summary <- cbind(
    data.frame(N_main_cohort = qc_sample$main_cohort_N,
               N_metabolomic_complete_M2 = nrow(datp0),
               N_direct_metabolites_listed = length(direct_metabolites),
               N_direct_metabolites_present = length(metabolite_vars0),
               N_metabolites_QC = length(metabolite_vars)),
    result$summary
)

s4_export_results(
    result, DIR_OUT, mediator_label = "metabolite",
    workbook_name = "DOSE_metabolite_IS_three_method_mediation_results.xlsx",
    candidate_prefix = "metabolite",
    extra_sheets = list("Sample QC" = qc_sample,
                        "Metabolite inventory" = metabolite_inventory,
                        "Metabolite QC" = metabolite_qc,
                        "M2 covariates" = data.frame(covariate = m2_covariates)),
    notes = c(
        "met.rds was previously cleaned by clean_biom(), which does not transform concentrations.",
        "Only direct absolute-concentration metabolites listed in met.lst are analyzed.",
        "Every analyzed metabolite is transformed by log(x+1) and then Z-standardized within its complete-case sample.",
        "Metabolite sampling may precede WebQ dietary assessment; interpret pathways as exploratory."
    )
)
s4_write_manifest(input_files, DIR_OUT)
print(result$summary)
message("Completed metabolite mediation analysis. Results: ", DIR_OUT)
}))

# ===========================================================================
# mechanism / metabolite_plot
# ===========================================================================
dose_register_step("mechanism", "metabolite_plot", expression({
# Bubble plot from the upstream metabolite association and Cox mediation results.
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({library(ggplot2); library(ggrepel)})
ROOT <- dose_output_dir()
OUT <- file.path(ROOT, "ukb/05mechanism/metabolite/M2")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
INPUT <- c(file.path(OUT, "01_all_metabolite_associations.csv"),
           file.path(OUT, "07_mediation_Cox_product.csv"))
d <- s4_read_mediation_table(OUT, "metabolite")
if (!nrow(d)) stop("No metabolite meets positive proportion and P_proportion < 0.05; no mediation bubble plot can be drawn.")
d$metabolite <- d$display_name
d$x <- abs(d$beta_a); d$y <- abs(d$beta_b)
d$direction <- factor(ifelse(d$beta_a<0,"Negative","Positive"), levels=c("Negative","Positive"))
# Use all selected rows. Labels identify the eight largest mediation proportions.
# The shared table builder converts the upstream proportion to a percentage once.
d <- d[order(-d$proportion_percent,d$metabolite), ]
labels <- head(d,8)
colors <- c(Negative="#F3AA83",Positive="#8CBED5")
p <- ggplot(d,aes(x,y)) +
    geom_point(aes(size=proportion_percent,fill=direction),shape=21,color="white",stroke=.25,alpha=.75) +
    geom_text_repel(data=labels,aes(label=metabolite),color="#303030",size=3.15,family="Arial",
       seed=12345,box.padding=.45,point.padding=.75,min.segment.length=0,
       segment.color="#AAAAAA",segment.size=.3,max.overlaps=Inf,max.iter=30000,max.time=5) +
    scale_fill_manual(values=colors,drop=FALSE,name=expression("DOSE-metabolite association ("*beta[a]*")")) +
    scale_size_area(max_size=13,breaks=c(.5,1,2,4),limits=c(0,max(4.5, d$proportion_percent)),
       name="Proportion mediated (%)") +
    scale_x_continuous(limits=c(0,NA),breaks=seq(0,.1,.02),expand=expansion(mult=c(.02,.15))) +
    scale_y_continuous(limits=c(0,NA),expand=expansion(mult=c(.02,.13))) +
    labs(x=expression("Absolute value of "*beta[a]*" coef."),
         y=expression("Absolute value of "*beta[b]*" coef.")) +
    guides(fill=guide_legend(order=1,override.aes=list(size=4,alpha=1)),
           size=guide_legend(order=2,override.aes=list(fill="#333333",alpha=1))) +
    theme_classic(base_size=12,base_family="Arial") +
    theme(axis.line=element_line(linewidth=.65),axis.ticks=element_line(linewidth=.5),
          axis.text=element_text(color="#202020",size=11),axis.title=element_text(size=13),
          axis.title.x=element_text(margin=margin(t=10)),axis.title.y=element_text(margin=margin(r=10)),
          legend.position="right",legend.title=element_text(size=10.5),legend.text=element_text(size=10),
          legend.spacing.y=grid::unit(.6,"cm"),plot.margin=margin(15,16,10,10))
stem <- file.path(OUT,"Fig_DOSE_metabolite_mediation_beta_bubble")
bounded_ggsave(paste0(stem,".pdf"),p,width=8.6,height=5.8,device=grDevices::cairo_pdf,bg="white")
bounded_ggsave(paste0(stem,".svg"),p,width=8.6,height=5.8,device=svglite::svglite,bg="white")
bounded_ggsave(paste0(stem,".png"),p,width=8.6,height=5.8,dpi=600,bg="white")
write.csv(d,paste0(stem,"_plot_data.csv"),row.names=FALSE)
writeLines(c("Metabolite mediation bubble plot", "",
    paste0("Sources: ", paste(INPUT, collapse="; ")),
    "X: absolute beta_a (DOSE-metabolite association); Y: absolute beta_b (metabolite-IS association), from a_estimate and b_adjusted_estimate in the Cox mediation results.",
    "Color reflects the SIGN OF beta_a, not beta_b or the indirect effect: negative orange #F3AA83; positive blue #8CBED5.",
    "Bubble area is proportional to mediation percentage. The upstream Proportion_mediated is multiplied by 100 once.",
    "All selected metabolites have positive mediation proportions and unadjusted proportion P<0.05. This figure does not imply that they all pass FDR correction.",
    "Eight largest mediation proportions labelled. No axis break. The plot displays coefficients, not hazard ratios.",
    paste("Number of plotted metabolites:",nrow(d)),paste("Source MD5:",tools::md5sum(INPUT))
),paste0(stem,"_README.txt"))
print(table(d$direction)); print(d[c("metabolite","proportion_percent","beta_a","beta_b")]);message("Saved: ",stem)
}))

# ===========================================================================
# mechanism / inflammation
# ===========================================================================
dose_register_step("mechanism", "inflammation", expression({
# DOSE -> inflammatory marker -> incident ischemic stroke (IS)
# Inflammation-specific import/QC; shared models are in 0.common.R.

options(stringsAsFactors = FALSE)

DIR_PROJECT <- dose_output_dir()
DIR_UKB <- dose_ukb_dir()
DIR_DATA <- file.path(DIR_PROJECT, "ukb/00data")
DIR_SCORE <- file.path(DIR_PROJECT, "ukb/04score")
DIR_OUT <- file.path(DIR_PROJECT, "ukb/05mechanism/inflammation/M2")

FILE_FUNCTIONS <- file.path(dose_script_dir(), "0.common.R")
FILE_COHORT <- file.path(DIR_DATA, "dat.base.imp.rds")
FILE_CONFIG <- file.path(DIR_DATA, "analysis_config.rds")
FILE_SCORE <- file.path(DIR_SCORE, "dose_score.rds")
FILE_INFLAM <- file.path(DIR_UKB, "Rdata/bbc.rds")

OUTCOME <- "cvd_stroke_i"
TIME_VAR <- paste0(OUTCOME, ".t2e")
EVENT_VAR <- paste0(OUTCOME, ".Yt2e")
MIN_N <- 1000L
MIN_EVENTS <- 20L
MEDIATION_SIMS <- 100L
SEED <- 12345L

if (!file.exists(FILE_FUNCTIONS)) stop("Missing shared functions: ", FILE_FUNCTIONS)

s4_load_packages()
dir.create(DIR_OUT, recursive = TRUE, showWarnings = FALSE)

input_files <- c(functions = FILE_FUNCTIONS, cohort = FILE_COHORT,
                 config = FILE_CONFIG, dose = FILE_SCORE)
if (any(!file.exists(input_files))) {
    stop("Missing input file(s): ",
         paste(input_files[!file.exists(input_files)], collapse = "; "))
}

cfg <- readRDS(FILE_CONFIG)
dat <- s4_validate_ids(as.data.frame(readRDS(FILE_COHORT)), "dat")
dose_checked <- s4_validate_dose(readRDS(FILE_SCORE))
dose <- dose_checked$data
dose_meta <- dose_checked$meta

if (!OUTCOME %in% cfg$Ys ||
    !identical(cfg$primary_development_outcome, OUTCOME) ||
    !identical(cfg$primary_model, "M2")) {
    stop("analysis_config.rds does not define cvd_stroke_i/M2 as primary.")
}
m2_covariates <- s4_get_m2_covariates(cfg)

# Five measured markers plus four derived ratios. Prefer the prepared cohort.
# If a whole measured column is absent, supplement only that column from bbc.rds;
# do not fill individual missing observations from another source.
base_inflammation_vars <- c(
    "bc_NEUT", "bc_PLT", "bc_MONO", "bc_LYMPH", "bb_CRP"
)
inflammation_vars0 <- c(
    base_inflammation_vars, "LMR", "SII", "NLR", "PLR"
)
inflam <- data.frame(eid = dat$eid)
source_qc <- data.frame(
    inflammation = base_inflammation_vars,
    source = "prepared cohort", stringsAsFactors = FALSE
)
missing_fields <- setdiff(base_inflammation_vars, names(dat))
if (length(missing_fields)) {
    if (!file.exists(FILE_INFLAM)) {
        stop("Missing cohort inflammation fields: ",
             paste(missing_fields, collapse = ", "),
             "; fallback file missing: ", FILE_INFLAM)
    }
    bbc <- s4_validate_ids(as.data.frame(readRDS(FILE_INFLAM)), "bbc")
    if (!all(missing_fields %in% names(bbc))) {
        stop("bbc.rds lacks: ",
             paste(setdiff(missing_fields, names(bbc)), collapse = ", "))
    }
    idx <- match(inflam$eid, bbc$eid)
    for (v in missing_fields) inflam[[v]] <- bbc[[v]][idx]
    source_qc$source[source_qc$inflammation %in% missing_fields] <- FILE_INFLAM
    input_files <- c(input_files, inflammation_fallback = FILE_INFLAM)
    rm(bbc, idx)
}
for (v in intersect(base_inflammation_vars, names(dat))) inflam[[v]] <- dat[[v]]
if (any(!vapply(inflam[base_inflammation_vars],
                function(x) is.numeric(x) || is.integer(x), logical(1)))) {
    stop("Inflammation fields must be numeric.")
}

invalid_qc <- dplyr::bind_rows(lapply(base_inflammation_vars, function(v) {
    x <- suppressWarnings(as.numeric(inflam[[v]]))
    data.frame(inflammation = v, N_missing = sum(is.na(x)),
               N_infinite = sum(is.infinite(x)),
               N_negative = sum(is.finite(x) & x < 0),
               N_zero = sum(is.finite(x) & x == 0))
}))
for (v in base_inflammation_vars) {
    x <- suppressWarnings(as.numeric(inflam[[v]]))
    x[!is.finite(x) | x < 0] <- NA_real_
    inflam[[v]] <- x
}
s4_safe_ratio <- function(num, den) {
    out <- rep(NA_real_, length(num))
    ok <- is.finite(num) & is.finite(den) & num >= 0 & den > 0
    out[ok] <- num[ok] / den[ok]
    out[!is.finite(out)] <- NA_real_
    out
}
inflam$LMR <- s4_safe_ratio(inflam$bc_LYMPH, inflam$bc_MONO)
inflam$NLR <- s4_safe_ratio(inflam$bc_NEUT, inflam$bc_LYMPH)
inflam$PLR <- s4_safe_ratio(inflam$bc_PLT, inflam$bc_LYMPH)
sii_num <- with(inflam, ifelse(
    is.finite(bc_NEUT) & bc_NEUT >= 0 & is.finite(bc_PLT) & bc_PLT >= 0,
    bc_NEUT * bc_PLT, NA_real_
))
inflam$SII <- s4_safe_ratio(sii_num, inflam$bc_LYMPH)
inflam <- s4_validate_ids(inflam[c("eid", inflammation_vars0)], "inflam")
rm(sii_num)

definitions <- data.frame(
    inflammation = inflammation_vars0,
    definition = c(
        "Neutrophils", "Platelets", "Monocytes", "Lymphocytes",
        "C-reactive protein", "Lymphocytes / Monocytes",
        "Neutrophils * Platelets / Lymphocytes",
        "Neutrophils / Lymphocytes", "Platelets / Lymphocytes"
    ),
    transformation = "direct Z-standardization; no log transformation",
    stringsAsFactors = FALSE
)
inventory <- definitions %>% dplyr::left_join(source_qc, by = "inflammation") %>%
    dplyr::mutate(source = ifelse(is.na(source), "derived", source))

prepared <- s4_prepare_base(
    dat = dat, dose = dose, mediator_data = inflam,
    mediator_vars = inflammation_vars0, time_var = TIME_VAR,
    event_var = EVENT_VAR, covariates = m2_covariates,
    cohort_label = "inflammation"
)
datp0 <- prepared$data

inflammation_qc <- data.frame(
    inflammation = inflammation_vars0,
    is_numeric = vapply(datp0[inflammation_vars0],
                        function(x) is.numeric(x) || is.integer(x), logical(1)),
    N_nonmissing = vapply(datp0[inflammation_vars0],
                          function(x) sum(!is.na(x)), integer(1)),
    percent_missing = 100 * vapply(datp0[inflammation_vars0],
                                    function(x) mean(is.na(x)), numeric(1)),
    N_unique = vapply(datp0[inflammation_vars0],
                      function(x) dplyr::n_distinct(x, na.rm = TRUE), integer(1)),
    SD = vapply(datp0[inflammation_vars0], function(x) {
        if (is.numeric(x) || is.integer(x)) stats::sd(as.numeric(x), na.rm = TRUE)
        else NA_real_
    }, numeric(1))
)
inflammation_vars <- inflammation_qc %>% dplyr::filter(
    is_numeric, N_nonmissing >= MIN_N, N_unique > 1L, is.finite(SD), SD > 0
) %>% dplyr::pull(inflammation)
if (!length(inflammation_vars)) stop("No inflammatory marker passed QC.")

qc_sample <- data.frame(
    main_cohort_N = nrow(dat), inflammation_source_N = nrow(inflam),
    overlap_N = prepared$overlap_N, complete_M2_N = nrow(datp0),
    complete_M2_events = sum(datp0[[EVENT_VAR]]),
    event_percent = 100 * mean(datp0[[EVENT_VAR]]),
    dose_mean = mean(datp0$dose), dose_SD = stats::sd(datp0$dose),
    dose_min = min(datp0$dose), dose_max = max(datp0$dose),
    measured_inflammation_N = length(base_inflammation_vars),
    derived_inflammation_N = length(setdiff(inflammation_vars0,
                                             base_inflammation_vars)),
    inflammation_QC_N = length(inflammation_vars),
    transformation = "direct Z-standardization; no log transformation"
)
s4_write_csv(qc_sample, file.path(DIR_OUT, "00_QC_sample.csv"))
s4_write_csv(inflammation_qc,
             file.path(DIR_OUT, "00_QC_inflammation_variables.csv"))
s4_write_csv(inventory, file.path(DIR_OUT, "00_inflammation_inventory.csv"))
s4_write_csv(source_qc,
             file.path(DIR_OUT, "00_inflammation_field_sources.csv"))
s4_write_csv(invalid_qc,
             file.path(DIR_OUT, "00_inflammation_invalid_values.csv"))
s4_write_csv(definitions,
             file.path(DIR_OUT, "00_inflammation_definitions.csv"))
capture.output(summary(inflam[inflammation_vars0]),
               file = file.path(DIR_OUT, "00_inflammation_summary.txt"))
dose_save_rds(inflam, file.path(DIR_OUT, "inflammation_9.rds"))
s4_write_csv(data.frame(covariate = m2_covariates),
             file.path(DIR_OUT, "00_M2_covariates_used.csv"))

message("Inflammation complete M2 cohort: N=", nrow(datp0),
        "; IS events=", sum(datp0[[EVENT_VAR]]),
        "; inflammatory markers=", length(inflammation_vars))
rm(dat, inflam, dose, dose_checked, prepared)
invisible(gc())

result <- s4_run_analysis(
    dat = datp0, mediator_vars = inflammation_vars,
    transform_mediator = s4_transform_direct_z,
    mediator_label = "inflammation", time_var = TIME_VAR,
    event_var = EVENT_VAR, covariates = m2_covariates,
    min_n = MIN_N, min_events = MIN_EVENTS,
    sims = MEDIATION_SIMS, seed = SEED,
    fdr_family_n = length(inflammation_vars)
)
result$summary <- cbind(
    data.frame(N_main_cohort = qc_sample$main_cohort_N,
               N_inflammation_complete_M2 = nrow(datp0),
               N_measured_inflammation = length(base_inflammation_vars),
               N_derived_inflammation = length(setdiff(
                   inflammation_vars0, base_inflammation_vars)),
               N_inflammation_QC = length(inflammation_vars)),
    result$summary
)

s4_export_results(
    result, DIR_OUT, mediator_label = "inflammation",
    workbook_name = "DOSE_inflammation_IS_three_method_mediation_results.xlsx",
    candidate_prefix = "inflammation",
    extra_sheets = list("Sample QC" = qc_sample,
                        "Inflammation inventory" = inventory,
                        "Inflammation QC" = inflammation_qc,
                        "Invalid values" = invalid_qc,
                        "M2 covariates" = data.frame(covariate = m2_covariates)),
    notes = c(
        "Five measured markers plus LMR, NLR, PLR and SII are analyzed.",
        "Ratios are derived in original units before standardization; denominators must be >0.",
        "All nine markers are standardized directly; no log(x+1) transformation is applied.",
        "Inflammatory-marker sampling may precede WebQ dietary assessment; interpret pathways as exploratory."
    )
)
s4_write_manifest(input_files, DIR_OUT)
print(result$summary)
message("Completed inflammation mediation analysis. Results: ", DIR_OUT)
}))

# ===========================================================================
# mechanism / inflammation_plot
# ===========================================================================
dose_register_step("mechanism", "inflammation_plot", expression({
# Bubble plot from the upstream inflammatory-marker association and Cox mediation results.
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({library(ggplot2); library(ggrepel)})
ROOT <- dose_output_dir()
OUT <- file.path(ROOT, "ukb/05mechanism/inflammation/M2")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
INPUT <- c(file.path(OUT, "01_all_inflammation_associations.csv"),
           file.path(OUT, "07_mediation_Cox_product.csv"))
d <- s4_read_mediation_table(OUT, "inflammation")
if (!nrow(d)) stop("No inflammation meets positive proportion and P_proportion < 0.05; no mediation bubble plot can be drawn.")
d$inflammation <- d$display_name
d$x <- abs(d$beta_a); d$y <- abs(d$beta_b)
d$direction <- factor(ifelse(d$beta_a<0,"Negative","Positive"), levels=c("Negative","Positive"))
# Use all selected rows and label every inflammatory marker.
# The shared table builder converts the upstream proportion to a percentage once.
d <- d[order(-d$proportion_percent,d$inflammation), ]
labels <- d
colors <- c(Negative="#F3AA83",Positive="#8CBED5")
p <- ggplot(d,aes(x,y)) +
    geom_point(aes(size=proportion_percent,fill=direction),shape=21,color="white",stroke=.25,alpha=.75) +
    geom_text_repel(data=labels,aes(label=inflammation),color="#303030",size=3.15,family="Arial",
       seed=12345,box.padding=.45,point.padding=.75,min.segment.length=0,
       segment.color="#AAAAAA",segment.size=.3,max.overlaps=Inf,max.iter=30000,max.time=5) +
    scale_fill_manual(values=colors,drop=TRUE,name=expression("DOSE-marker association ("*beta[a]*")")) +
    scale_size_area(max_size=13,breaks=c(.5,1,2,4),limits=c(0,max(4.5, d$proportion_percent)),
       name="Proportion mediated (%)") +
    scale_x_continuous(limits=c(0,NA),breaks=seq(0,.1,.01),expand=expansion(mult=c(.02,.15))) +
    scale_y_continuous(limits=c(0,NA),expand=expansion(mult=c(.02,.13))) +
    labs(x=expression("Absolute value of "*beta[a]*" coef."),
         y=expression("Absolute value of "*beta[b]*" coef.")) +
    guides(fill=guide_legend(order=1,override.aes=list(size=4,alpha=1)),
           size=guide_legend(order=2,override.aes=list(fill="#333333",alpha=1))) +
    theme_classic(base_size=12,base_family="Arial") +
    theme(aspect.ratio=1,axis.line=element_line(linewidth=.65),axis.ticks=element_line(linewidth=.5),
          axis.text=element_text(color="#202020",size=11),axis.title=element_text(size=13),
          axis.title.x=element_text(margin=margin(t=10)),axis.title.y=element_text(margin=margin(r=10)),
          legend.position="bottom",legend.box="vertical",legend.title=element_text(size=10.5),legend.text=element_text(size=10),
          legend.spacing.y=grid::unit(.6,"cm"),plot.margin=margin(15,16,10,10))
stem <- file.path(OUT,"Fig_DOSE_inflammation_mediation_beta_bubble")
bounded_ggsave(paste0(stem,".pdf"),p,width=7.8,height=7.8,device=grDevices::cairo_pdf,bg="white")
bounded_ggsave(paste0(stem,".svg"),p,width=7.8,height=7.8,device=svglite::svglite,bg="white")
bounded_ggsave(paste0(stem,".png"),p,width=7.8,height=7.8,dpi=600,bg="white")
write.csv(d,paste0(stem,"_plot_data.csv"),row.names=FALSE)
writeLines(c("Inflammatory-marker mediation bubble plot", "",
    paste0("Sources: ", paste(INPUT, collapse="; ")),
    "X: absolute beta_a (DOSE-marker association); Y: absolute beta_b (inflammation-IS association), from a_estimate and b_adjusted_estimate in the Cox mediation results.",
    "Color reflects the SIGN OF beta_a, not beta_b or the indirect effect: negative orange #F3AA83; positive blue #8CBED5.",
    "Bubble area is proportional to mediation percentage. The upstream Proportion_mediated is multiplied by 100 once.",
    "All selected inflammatory markers have positive mediation proportions and unadjusted proportion P<0.05. This figure does not imply that they all pass FDR correction.",
    "All inflammatory markers labelled. No axis break. The plot displays coefficients, not hazard ratios.",
    "Square canvas (7.8 x 7.8 inches) and square plotting panel (aspect.ratio=1); axis display lengths are equal, not numeric units.",
    paste("Number of plotted inflammatory markers:",nrow(d)),paste("Source MD5:",tools::md5sum(INPUT))
),paste0(stem,"_README.txt"))
print(table(d$direction)); print(d[c("inflammation","proportion_percent","beta_a","beta_b")]);message("Saved: ",stem)
}))

# ===========================================================================
# mechanism / mri
# ===========================================================================
dose_register_step("mechanism", "mri", expression({
# DOSE MRI: observed outcomes, explicit age/centre source, separate test families.
# Defaults are absolute and do not depend on the working directory.
# Rewrites current MRI results; the run entry point archives previous outputs.
options(stringsAsFactors=FALSE)
PROJECT <- dose_output_dir()
UKB <- dose_ukb_dir()
RAW <- Sys.getenv("DOSE_MRI_RAW", unset=file.path(UKB,"rap/img.tab.gz"))
OUT <- file.path(PROJECT,"ukb/05mechanism/04mri")
# A cache must be explicitly selected after prepare.py verifies its provenance.
# The run entry point sets DOSE_MRI_RAW only after preparation succeeds.
DICT <- dose_mri_dictionary()
VISIT <- Sys.getenv("DOSE_MRI_VISIT", unset=file.path(OUT,"inputs/mri_visit_covariates.csv"))
MODEL <- "M2"
VISIT_MODE <- Sys.getenv("DOSE_MRI_VISIT_MODE", unset="baseline")
stopifnot(VISIT_MODE %in% c("baseline","imaging"))
MIN_N <- 100L
# User requested baseline age/centre. Set DOSE_MRI_VISIT_MODE=imaging only
# when intentionally switching to verified imaging-visit covariates.
# Remove age/site before adding the selected versions once, plus age squared.
# Audit 00_covariates_used.csv; edit this list if your actual names differ.
REPLACE_BASELINE_COVS <- c("age","age0","age_baseline","assessment_age", "centre", "center", "assessment_centre")
# FALSE matches the conventional lm SE used by the reference. TRUE uses HC3.
ROBUST_HC3 <- FALSE
for (p in c("data.table",if(ROBUST_HC3) "sandwich"))
    if (!requireNamespace(p,quietly=TRUE)) stop("Install package: ",p)
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
write_tab <- function(x,n) write.csv(x,file.path(OUT,n),row.names=FALSE,na="")
ids <- function(x,label) {
    x <- as.data.frame(x)
    if (!"eid" %in% names(x) || anyDuplicated(names(x))) stop("Invalid columns: ",label)
    x$eid <- as.character(x$eid)
    if (anyNA(x$eid)||any(!nzchar(x$eid))||anyDuplicated(x$eid)) stop("Invalid/duplicate eid: ",label)
    x
}
numeric_clean <- function(x) {
    if (is.factor(x)) x <- as.character(x)
    z <- suppressWarnings(as.numeric(x)); z[!is.finite(z)] <- NA_real_; z
}
standardize <- function(x) {
    s <- sd(x); if(!is.finite(s)||s<=0) stop("Zero or invalid SD")
    as.numeric(scale(x, center=TRUE, scale=TRUE))
}
complete <- function(d) {
    ok <- complete.cases(d)
    for(n in names(d)) if(is.numeric(d[[n]])) ok <- ok & is.finite(d[[n]])
    ok
}
D <- read.csv(DICT,check.names=FALSE)
required <- c("field_id","family","category","label","side","atlas_key","transform","adjust_icv")
if(!all(required %in% names(D))) stop("Missing dictionary columns")
D$field_id <- as.character(D$field_id)
if(anyNA(D$adjust_icv) || any(!tolower(as.character(D$adjust_icv)) %in% c("true","false"))) stop("Invalid adjust_icv flag")
D$adjust_icv <- tolower(as.character(D$adjust_icv))=="true"
expected <- c(FA27=27L,MD27=27L,ThicknessDK=66L,Volume139=139L,WMH=1L)
stopifnot(!anyNA(D$field_id),!anyDuplicated(D$field_id),identical(table(factor(D$family,levels=names(expected))) |> as.integer(),unname(expected)),nrow(D)==260L)
# Latest agreed analysis: all outcomes directly standardized, including WMH.
# Preserve the supplied dictionary setting rather than silently editing its source.
D$dictionary_transform <- D$transform
D$transform <- "none"
write_tab(D,"00_dictionary_effective.csv")
cfg <- readRDS(file.path(PROJECT,"ukb/00data/analysis_config.rds"))
dat <- ids(readRDS(file.path(PROJECT,"ukb/00data/dat.base.imp.rds")),"cohort")
s0 <- readRDS(file.path(PROJECT,"ukb/04score/dose_score.rds"))
components <- attr(s0,"component_columns",exact=TRUE)
if (!identical(attr(s0,"score_version",exact=TRUE),"dose") || !length(components) ||
    anyDuplicated(components)) stop("Missing/invalid primary DOSE metadata")
score_metadata <- attributes(s0)[setdiff(names(attributes(s0)),c("names","row.names","class"))]
write_tab(data.frame(score_version=score_metadata$score_version,
 primary_scoring_method=if(is.null(score_metadata$primary_scoring_method)) "unspecified" else score_metadata$primary_scoring_method,
 contains_alcohol=score_metadata$score_contains_alcohol,
 components=paste(components,collapse=";")),"00_score_metadata.csv")
score <- ids(s0,"score"); rm(s0)
exposures <- c("dose",components)
if(!all(exposures %in% names(score))||!all(vapply(score[exposures],is.numeric,logical(1)))) stop("Invalid score columns")
if(!isTRUE(all.equal(score$dose,rowSums(score[components]),check.attributes=FALSE))) stop("DOSE total != component sum")
if(!setequal(dat$eid,score$eid)) stop("Score/cohort ID mismatch")
dat[exposures] <- score[match(dat$eid,score$eid),exposures,drop=FALSE]
covs0 <- switch(MODEL,M1=cfg$covs.m1,M2=unique(c(cfg$covs.m1,cfg$covs.m2.add)),stop("Invalid model"))
if(!length(covs0)||anyNA(covs0)) stop("Missing M2 covariate configuration")
# Only omit alcohol when the score metadata explicitly says it contains alcohol.
contains_alcohol <- attr(score,"score_contains_alcohol",exact=TRUE)
if(is.null(contains_alcohol)) contains_alcohol <- score_metadata$score_contains_alcohol
if(!is.logical(contains_alcohol)||length(contains_alcohol)!=1L||is.na(contains_alcohol)) stop("Explicit score_contains_alcohol metadata required")
remove <- REPLACE_BASELINE_COVS
if(contains_alcohol) remove <- unique(c(remove,cfg$covs.alcohol,"alcohol.c"))
covs <- setdiff(unique(covs0),remove)
if(!all(covs %in% names(dat))||any(covs %in% exposures)) stop("Missing/invalid cohort covariates")
write_tab(data.frame(variable=covs0,retained=covs0 %in% covs),"00_inherited_covariate_audit.csv")

# Read only explicitly selected instance-2 fields; never load all modalities for KNN.
if(!file.exists(RAW)) stop("Raw MRI table not found: ",RAW)
con <- if(grepl("\\.gz$",RAW)) gzfile(RAW,"rt") else file(RAW,"rt")
header <- strsplit(readLines(con,n=1L),"\t",fixed=TRUE)[[1]]; close(con)
D$column <- paste0("p",D$field_id,"_i2")
D$available <- D$column %in% header
write_tab(D,"00_MRI_field_mapping.csv")
if(VISIT_MODE=="imaging" && !file.exists(VISIT) && !all(c("p21003_i2","p54_i2") %in% header))
    stop("Missing MRI visit covariates. Run mechanism.py --mode visit first; expected: ",VISIT)
needed <- unique(c(D$column,"p26521_i2","p21003_i2","p54_i2","p53_i2"))
img <- ids(data.table::fread(RAW,select=c("eid",intersect(needed,header))),"raw MRI")
for(n in setdiff(names(img),c("eid","p53_i2"))) img[[n]] <- numeric_clean(img[[n]])
# Absent fields remain planned tests with failed status; none silently dropped.
for(n in setdiff(D$column,names(img))) img[[n]] <- NA_real_
if(!"p26521_i2" %in% names(img)) img$p26521_i2 <- NA_real_

# Optional verified CSV: eid, imaging_age, imaging_site; optional dates YYYY-MM-DD.
if(VISIT_MODE=="baseline") {
    if(!all(c("age","center") %in% names(dat))) stop("Baseline age/center missing")
    visit <- data.frame(eid=dat$eid,imaging_age=dat$age,imaging_site=dat$center)
    visit_source <- "baseline age/center (user requested); NOT imaging visit"
} else if(file.exists(VISIT)) {
    visit <- ids(data.table::fread(VISIT),"MRI visit covariates")
    if(!all(c("imaging_age","imaging_site") %in% names(visit))) stop("Visit CSV needs eid, imaging_age, imaging_site")
    visit_source <- VISIT
} else if(all(c("p21003_i2","p54_i2") %in% names(img))) {
    visit <- data.frame(eid=img$eid,imaging_age=img$p21003_i2,imaging_site=img$p54_i2)
    visit_source <- "raw imaging visit p21003_i2/p54_i2"
    if("p53_i2" %in% names(img)) visit$imaging_date <- as.character(img$p53_i2)
} else {
    stop("Need verified imaging visit age/site. Supply mri_visit_covariates.csv with eid, imaging_age (p21003_i2), imaging_site (p54_i2), or include these fields in the raw table. No baseline substitution.")
}
visit$imaging_age <- numeric_clean(visit$imaging_age)
visit$imaging_age[visit$imaging_age<=0] <- NA_real_
visit$imaging_site <- as.character(visit$imaging_site)
visit$imaging_site[!nzchar(trimws(visit$imaging_site)) | grepl("^-",visit$imaging_site)] <- NA_character_
joined <- dat[dat$eid %in% img$eid & dat$eid %in% visit$eid,c("eid",exposures,covs),drop=FALSE]
ii <- match(joined$eid,img$eid); iv <- match(joined$eid,visit$eid)
if(!nrow(joined)) stop("No cohort / MRI / visit overlap")
base <- data.frame(eid=joined$eid)
cn <- paste0("cov",seq_along(covs))
for(j in seq_along(covs)) base[[cn[j]]] <- joined[[covs[j]]]
base$imaging_age <- visit$imaging_age[iv]
base$imaging_site <- factor(visit$imaging_site[iv])
base$icv <- img$p26521_i2[ii]
base$icv[base$icv<=0] <- NA_real_
age_labels <- if(VISIT_MODE=="baseline") c("baseline_age","baseline_age_squared","baseline_center") else c("imaging_age","imaging_age_squared","imaging_site")
write_tab(data.frame(variable=c(covs,age_labels,"ICV (volume/WMH only)"),
                    source=c(rep("cohort M2",length(covs)),rep(visit_source,3),"p26521_i2")),"00_covariates_used.csv")
# Dates are not inferred from baseline age. If both supplied, enforce diet < MRI.
chronology <- "not verified: diet/MRI dates not both supplied"
if(all(c("diet_date","imaging_date") %in% names(visit))) {
    dd <- as.Date(visit$diet_date[iv]); md <- as.Date(visit$imaging_date[iv])
    ok <- !is.na(dd)&!is.na(md)&dd<md
    write_tab(data.frame(overlap=nrow(joined),eligible=sum(ok),excluded=sum(!ok)),"00_chronology.csv")
    joined <- joined[ok,,drop=FALSE];base <- base[ok,,drop=FALSE];ii <- ii[ok]
    chronology <- "enforced diet_date < imaging_date; missing dates excluded"
}
if(!nrow(joined)) stop("No participants after chronology restriction")
qc <- data.frame(field_id=D$field_id,family=D$family,available=D$available,
                 N_overlap=nrow(joined),N_observed=vapply(D$column,function(n)sum(is.finite(img[[n]][ii])),integer(1)))
write_tab(qc,"00_observed_MRI_counts.csv")
# Participants with no observed outcome cannot enter any of these models.
# This is NOT a requirement to have all modalities: one observed field suffices.
has_any <- Reduce(`|`,lapply(D$column,function(n)is.finite(img[[n]][ii])))
write_tab(data.frame(overlap=nrow(joined),any_observed=sum(has_any),none_observed=sum(!has_any)),"00_analysis_eligibility.csv")
joined <- joined[has_any,,drop=FALSE];base <- base[has_any,,drop=FALSE]
img <- img[ii[has_any],unique(c("eid",D$column)),drop=FALSE];ii <- seq_len(nrow(img))
rm(dat,score);invisible(gc())

run_pair <- function(e,k) {
    if(!D$available[k]) stop("Field absent from raw table")
    terms <- c(cn,"age_z","age_z2")
    bcols <- c(cn,"imaging_age","imaging_site",if(D$adjust_icv[k]) "icv")
    d <- base[bcols];d$x <- joined[[e]];d$y <- img[[D$column[k]]][ii]
    d$y[d$y<0] <- NA_real_
    if(D$family[k]=="ThicknessDK") d$y[d$y==0] <- NA_real_
    if(D$family[k]=="FA27") d$y[d$y>1] <- NA_real_
    # All MRI outcomes, including WMH: direct scale(), no log transformation.
    d <- droplevels(d[complete(d),,drop=FALSE])
    if(nrow(d)<MIN_N) stop("Insufficient observed complete cases: ",nrow(d))
    xsd <- sd(d$x);ysd <- sd(d$y)
    d$xz <- standardize(d$x);d$yz <- standardize(d$y)
    d$age_z <- standardize(d$imaging_age);d$age_z2 <- d$age_z^2
    site_fixed <- nlevels(d$imaging_site)==1L
    if(!site_fixed) terms <- c(terms,"imaging_site")
    constant <- cn[vapply(d[cn],function(x)length(unique(x))<2,logical(1))]
    if(length(constant)) stop("Constant inherited covariate: ",paste(covs[match(constant,cn)],collapse=", "))
    if(D$adjust_icv[k]) {d$icv <- standardize(d$icv);terms <- c(terms,"icv")}
    fit <- lm(reformulate(c("xz",terms),response="yz"),data=d,na.action=na.fail)
    if(fit$rank<length(coef(fit))) stop("Rank-deficient model")
    beta <- unname(coef(fit)["xz"])
    V <- if(ROBUST_HC3) sandwich::vcovHC(fit,type="HC3") else vcov(fit)
    se <- sqrt(V["xz","xz"]);df <- fit$df.residual
    if(!is.finite(beta)||!is.finite(se)||se<=0||df<=0) stop("Invalid coefficient or SE")
    ci <- beta+c(-1,1)*qt(.975,df)*se
    c(N=nrow(d),beta=beta,SE=se,LCI=ci[1],HCI=ci[2],P=2*pt(-abs(beta/se),df),
      exposure_SD=xsd,MRI_SD=ysd,beta_per_score_point=beta/xsd,site_constant=as.numeric(site_fixed))
}
mri_cache <- dose_result_cache(file.path(OUT, ".checkpoints/fits"),
    list(base, joined[exposures], img, ii, D, covs, cn, MODEL, VISIT_MODE, MIN_N, ROBUST_HC3),
    packages=c("stats", if (ROBUST_HC3) "sandwich"))
invisible(gc())
rows <- list()
for(e in exposures) {
    message("START MRI ",e,"; planned=",nrow(D))
    rows[[e]] <- mri_cache(e, {
    block <- vector("list", nrow(D))
    for(k in seq_len(nrow(D))) {
        err <- "";warn <- character()
        z <- tryCatch(withCallingHandlers(run_pair(e,k),warning=function(w){warn <<- c(warn,conditionMessage(w));invokeRestart("muffleWarning")}),error=function(x){err <<- conditionMessage(x);NULL})
        if(is.null(z)) z <- setNames(rep(NA_real_,10),c("N","beta","SE","LCI","HCI","P","exposure_SD","MRI_SD","beta_per_score_point","site_constant"))
        block[[k]] <- cbind(data.frame(exposure=e,exposure_type=if(e=="dose")"total" else "component",model=MODEL,age_center_source=VISIT_MODE),D[k,],as.data.frame(as.list(z)),data.frame(status=if(nzchar(err))"failed" else "ok",error=err,warnings=paste(unique(warn),collapse=" | ")))
    }
    do.call(rbind, block)
    })
    message("DONE MRI ",e)
    # Checkpoint only completed exposure blocks, with explicit partial status.
    dose_save_rds(list(completed_exposure=e,partial_results=do.call(rbind,rows)),file.path(OUT,"checkpoint_partial.rds"))
}
r <- do.call(rbind,rows);rownames(r)<-NULL
r$FDR_family <- r$FDR_exposure_all <- r$FDR_components_family <- NA_real_
for(e in exposures) {
    i <- which(r$exposure==e)
    r$FDR_exposure_all[i] <- p.adjust(r$P[i],"BH",n=nrow(D))
    for(f in unique(D$family)) {
        j <- i[r$family[i]==f]
        r$FDR_family[j] <- p.adjust(r$P[j],"BH",n=sum(D$family==f))
    }
}
for(f in unique(D$family)) {
    j <- which(r$family==f & r$exposure_type=="component")
    r$FDR_components_family[j] <- p.adjust(r$P[j],"BH",n=sum(D$family==f)*length(components))
}
r$FDR_all_tests <- p.adjust(r$P,"BH",n=nrow(D)*length(exposures))
write_tab(r,"01_all_linear_associations.csv")
write_tab(r[r$exposure=="dose",],"02_DOSE_MRI.csv")
write_tab(r[r$exposure_type=="component",],"03_components_MRI.csv")
write_tab(r[r$status!="ok" | nzchar(r$warnings),],"05_errors_warnings.csv")
counts <- do.call(rbind,lapply(split(r,interaction(r$exposure,r$family,drop=TRUE)),function(q)
    data.frame(exposure=q$exposure[1],family=q$family[1],planned=nrow(q),completed=sum(q$status=="ok"),P05=sum(q$P<.05,na.rm=TRUE),FDR05=sum(q$FDR_family<.05,na.rm=TRUE))))
write_tab(counts,"06_summary_counts.csv")
dose_save_rds(list(results=r,dictionary=D,score_metadata=score_metadata),file.path(OUT,"DOSE_MRI_revised.rds"))
writeLines(c("Observed MRI outcomes and ICV; no KNN and no blanket across-modality missingness filter.",
 "139 volume outcomes retained; 66 local DK thickness fields available in BIG40 (temporal poles absent); two global mean fields excluded.",
 "27 FA + 27 MD outcomes and WMH; all outcomes directly scale()-standardized, without log transformation. Do not confuse absent fields with nonsignificant results.",
 paste("M2 with age, age squared and centre from:",visit_source,"; ICV for volume/WMH only."),
 "A single observed centre is explicitly recorded as site_constant=1; its intercept absorbs this constant term.",
 "Primary FDR within exposure and prespecified modality family; broader corrections provided.",
 "Component predictors are scores, not intake amounts. Per-score-point beta also reported.",
 paste("Chronology:",chronology),paste("HC3:",ROBUST_HC3),
 "Associations only; no claim of MRI mediation of stroke. Input cohort may already have imputed covariates.",
 "No automated outlier trimming; inspect MRI QC and residual diagnostics before manuscript use."
),file.path(OUT,"README_analysis.txt"))
capture.output(sessionInfo(),file=file.path(OUT,"sessionInfo.txt"))
if(file.exists(file.path(OUT,"checkpoint_partial.rds"))) unlink(file.path(OUT,"checkpoint_partial.rds"))
print(counts)
if(all(r$status=="failed")) stop("All models failed; inspect error table")
message("DONE output: ",OUT)
}))

# ===========================================================================
# mechanism / mri_refit
# ===========================================================================
dose_register_step("mechanism", "mri_refit", expression({
# Independent refit of one real-data association using named covariates.
# Run after mechanism.R. Does not overwrite the primary association tables.
root <- dose_output_dir()
out <- file.path(root,"ukb/05mechanism/04mri")
cfg <- readRDS(file.path(root,"ukb/00data/analysis_config.rds"))
covs <- setdiff(unique(c(cfg$covs.m1,cfg$covs.m2.add)),c("age",cfg$covs.alcohol,"alcohol.c"))
d <- as.data.frame(readRDS(file.path(root,"ukb/00data/dat.base.imp.rds")))[c("eid","age","center",covs)]
s <- as.data.frame(readRDS(file.path(root,"ukb/04score/dose_score.rds")))[c("eid","dose")]
m <- as.data.frame(data.table::fread(file.path(out,"inputs/mri_selected.tsv"),select=c("eid","p25782_i2","p26521_i2")))
d <- merge(merge(d,s,by="eid"),m,by="eid")
d$p25782_i2[d$p25782_i2<0] <- NA_real_
d$p26521_i2[d$p26521_i2<=0] <- NA_real_
d <- d[complete.cases(d),]
d$center <- factor(d$center)
f <- reformulate(c("scale(dose)",covs,"scale(age)","I(as.numeric(scale(age))^2)","center","scale(p26521_i2)"),response="scale(p25782_i2)")
fit <- lm(f,data=d)
v <- coef(summary(fit))["scale(dose)",]
r <- read.csv(file.path(out,"02_DOSE_MRI.csv"))
r <- r[r$field_id==25782,]
stopifnot(nrow(r)==1L,r$N==nobs(fit),abs(r$beta-v[1])<1e-10,abs(r$SE-v[2])<1e-10,abs(r$P-v[4])<1e-10)
ci <- confint(fit,"scale(dose)")
stopifnot(abs(r$LCI-ci[1])<1e-10,abs(r$HCI-ci[2])<1e-10)
report <- c("PASS: independent raw-data refit of DOSE / field 25782 matches N, beta, SE, P and 95% CI.",
            paste("N =",nobs(fit),"residual df =",df.residual(fit)),paste("Formula:",paste(deparse(f),collapse=" ")))
writeLines(report,file.path(out,"verification_refit.txt"));cat(paste(report,collapse="\n"),"\n")
}))

# ===========================================================================
# mechanism / mri_wmh
# ===========================================================================
dose_register_step("mechanism", "mri_wmh", expression({
# Marginally adjusted mean prediction from the SAME linear WMH model as MRI analysis.
root <- dose_output_dir()
out <- file.path(root,"ukb/05mechanism/04mri")
dest <- commandArgs(trailingOnly=TRUE)
dest <- if(length(dest)) dest[1] else file.path(out,"figures/audit")
dir.create(dest,recursive=TRUE,showWarnings=FALSE)
cfg <- readRDS(file.path(root,"ukb/00data/analysis_config.rds"))
covs <- setdiff(unique(c(cfg$covs.m1,cfg$covs.m2.add)),c("age",cfg$covs.alcohol,"alcohol.c"))
d <- as.data.frame(readRDS(file.path(root,"ukb/00data/dat.base.imp.rds")))[c("eid","age","center",covs)]
s <- as.data.frame(readRDS(file.path(root,"ukb/04score/dose_score.rds")))[c("eid","dose")]
m <- as.data.frame(data.table::fread(file.path(out,"inputs/mri_selected.tsv"),select=c("eid","p25781_i2","p26521_i2")))
d <- merge(merge(d,s,by="eid"),m,by="eid")
d$p25781_i2[d$p25781_i2<0] <- NA_real_
d$p26521_i2[d$p26521_i2<=0] <- NA_real_
ok <- complete.cases(d)
for(n in names(d)) if(is.numeric(d[[n]])) ok <- ok & is.finite(d[[n]])
d <- droplevels(d[ok,]);d$center <- factor(d$center)
xmean <- mean(d$dose);xsd <- sd(d$dose);ymean <- mean(d$p25781_i2);ysd <- sd(d$p25781_i2)
d$xz <- as.numeric(scale(d$dose));d$yz <- as.numeric(scale(d$p25781_i2))
d$agez <- as.numeric(scale(d$age));d$agez2 <- d$agez^2;d$icvz <- as.numeric(scale(d$p26521_i2))
f <- reformulate(c("xz",covs,"agez","agez2",if(nlevels(d$center)>1) "center","icvz"),response="yz")
fit <- lm(f,data=d,na.action=na.fail)
r <- read.csv(file.path(out,"02_DOSE_MRI.csv"));r <- r[r$field_id==25781,]
v <- coef(summary(fit))["xz",];ci <- confint(fit,"xz")
stopifnot(nrow(r)==1,r$age_center_source=="baseline",r$N==nobs(fit),
 abs(r$beta-v[1])<1e-10,abs(r$SE-v[2])<1e-10,abs(r$P-v[4])<1e-10,
 abs(r$LCI-ci[1])<1e-10,abs(r$HCI-ci[2])<1e-10)
# Average model design over the empirical covariate distribution, replacing only exposure.
# t-based CI for the adjusted mean; not a prediction interval for an individual.
grid <- seq(quantile(d$dose,.01),quantile(d$dose,.99),length.out=121)
X <- model.matrix(fit);A <- matrix(rep(colMeans(X),each=length(grid)),nrow=length(grid))
colnames(A) <- colnames(X);A[,"xz"] <- (grid-xmean)/xsd
mu <- drop(A%*%coef(fit));se <- sqrt(rowSums((A%*%vcov(fit))*A));crit <- qt(.975,df.residual(fit))
pred <- data.frame(dose=grid,estimate_z=mu,lower_z=mu-crit*se,upper_z=mu+crit*se,
 estimate_mm3=ymean+ysd*mu,lower_mm3=ymean+ysd*(mu-crit*se),upper_mm3=ymean+ysd*(mu+crit*se))
write.csv(pred,file.path(dest,"WMH_adjusted_predictions.csv"),row.names=FALSE)
write.csv(data.frame(N=nobs(fit),beta=r$beta,SE=r$SE,P=r$P,FDR_family=r$FDR_family,
 dose_mean=xmean,dose_SD=xsd,WMH_mean=ymean,WMH_SD=ysd,grid_min=min(grid),grid_max=max(grid),
 method="Marginal mean; empirical covariate distribution; linear model; t-based 95% CI",
 interval="Mean confidence interval, not individual prediction interval",
 range="Observed DOSE 1st-99th percentiles",age_center_source="baseline"),file.path(dest,"WMH_prediction_metadata.csv"),row.names=FALSE)
cat("PASS: WMH refit matches existing N, beta, SE, P and CI. Marginal predictions exported.\n")
}))

# ===========================================================================
# mechanism / sem
# ===========================================================================
dose_register_step("mechanism", "sem", expression({
# SEM selected from the upstream Cox mediation results: all / top 10 / top 5 metabolites,
# each paired with all eligible inflammatory markers. No manuscript workbook input.
run_cox_sem <- function() {
options(stringsAsFactors = FALSE)

DIR_PROJECT <- dose_output_dir()
DIR_UKB <- dose_ukb_dir()
DIR_DATA <- file.path(DIR_PROJECT, "ukb/00data")
DIR_SCORE <- file.path(DIR_PROJECT, "ukb/04score")
DIR_MET <- file.path(DIR_PROJECT, "ukb/05mechanism/metabolite/M2")
DIR_INFLA <- file.path(DIR_PROJECT, "ukb/05mechanism/inflammation/M2")
DIR_OUT <- file.path(
    DIR_PROJECT, "ukb/05mechanism/sem/M2_Cox"
)

FILE_FUNCTIONS <- file.path(dose_script_dir(), "0.common.R")
FILE_COHORT <- file.path(DIR_DATA, "dat.base.imp.rds")
FILE_CONFIG <- file.path(DIR_DATA, "analysis_config.rds")
FILE_SCORE <- file.path(DIR_SCORE, "dose_score.rds")
FILE_MET <- file.path(DIR_UKB, "Rdata/met.rds")
FILE_INFLA <- file.path(DIR_INFLA, "inflammation_9.rds")
FILE_MET_ASSOC <- file.path(DIR_MET, "01_all_metabolite_associations.csv")
FILE_INFLA_ASSOC <- file.path(DIR_INFLA,
                              "01_all_inflammation_associations.csv")
FILE_MET_COX <- file.path(DIR_MET, "07_mediation_Cox_product.csv")
FILE_INFLA_COX <- file.path(DIR_INFLA, "07_mediation_Cox_product.csv")

OUTCOME <- "cvd_stroke_i"
EVENT_VAR <- paste0(OUTCOME, ".Yt2e")
# Rank eligible metabolites by P_proportion ascending, with marker ID as the tie-breaker.
# Eligibility: upstream candidate, positive proportion, unadjusted P_proportion < 0.05.
# No additional FDR filter. Alternative ranking: "P_ACME".
MET_RANK_BY <- "P_proportion"
SELECTION_RULE <- paste0("Cox_positive_proportion_P05_", MET_RANK_BY)
MEASUREMENT_SPEC <- "cox_selected"
STOP_ON_IMPROPER_SOLUTION <- FALSE
HEYWOOD_TOLERANCE <- 1e-6
MIN_LATENT_INDICATORS <- 3L
DRAW_SEM_PLOT <- FALSE
MAX_PLOT_INDICATORS <- 40L
SEED <- 12345L

input_files <- c(
    functions = FILE_FUNCTIONS, cohort = FILE_COHORT,
    config = FILE_CONFIG, dose = FILE_SCORE,
    metabolite = FILE_MET, inflammation = FILE_INFLA,
    metabolite_associations = FILE_MET_ASSOC,
    inflammation_associations = FILE_INFLA_ASSOC,
    metabolite_Cox = FILE_MET_COX,
    inflammation_Cox = FILE_INFLA_COX
)
check_files <- if ("--selection-only" %in% commandArgs(trailingOnly = TRUE)) {
    input_files[c("functions", "metabolite_associations", "inflammation_associations", "metabolite_Cox", "inflammation_Cox")]
} else input_files
if (any(!file.exists(check_files))) {
    stop(
        "Missing input file(s). Run mechanism.R and mechanism.R first: ",
        paste(check_files[!file.exists(check_files)], collapse = "; ")
    )
}

required_packages <- if ("--selection-only" %in% commandArgs(trailingOnly = TRUE)) character() else c("dplyr", "lavaan", "openxlsx")
missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
    stop("Install required packages first: ",
         paste(missing_packages, collapse = ", "))
}
if (!"--selection-only" %in% commandArgs(trailingOnly = TRUE)) suppressPackageStartupMessages(library(dplyr))
dir.create(DIR_OUT, recursive = TRUE, showWarnings = FALSE)
set.seed(SEED)

make_alias_map <- function(selection, prefix) {
    selection$alias <- sprintf("%s%03d", prefix, seq_len(nrow(selection)))
    selection
}

make_covariates <- function(data, covariates) {
    form <- stats::reformulate(covariates)
    mm <- stats::model.matrix(form, data = data)
    if ("(Intercept)" %in% colnames(mm)) {
        mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
    }
    if (!ncol(mm)) stop("M2 design matrix contains no usable covariate column.")
    keep <- vapply(seq_len(ncol(mm)), function(j) {
        x <- mm[, j]
        all(is.finite(x)) && stats::sd(x) > 0
    }, logical(1))
    mm <- mm[, keep, drop = FALSE]
    if (!ncol(mm)) stop("All M2 design-matrix columns have zero variance.")

    # Keep 0/1 indicators unchanged and Z-standardize other numeric columns.
    for (j in seq_len(ncol(mm))) {
        u <- unique(mm[, j])
        if (!all(u %in% c(0, 1))) mm[, j] <- as.numeric(scale(mm[, j]))
    }
    original <- colnames(mm)
    alias <- sprintf("Cov%03d", seq_len(ncol(mm)))
    colnames(mm) <- alias
    list(
        data = as.data.frame(mm, check.names = FALSE),
        map = data.frame(alias = alias, original = original,
                         stringsAsFactors = FALSE)
    )
}

make_domain <- function(data, aliases, latent_name, score_name,
                        min_indicators = 3L) {
    if (length(aliases) >= min_indicators) {
        return(list(
            data = data,
            variable = latent_name,
            measurement = paste0(
                latent_name, " =~ ", paste(aliases, collapse = " + ")
            ),
            model_type = "latent factor",
            n_indicators = length(aliases)
        ))
    }
    # A one- or two-indicator latent factor is generally not identified without
    # extra constraints. Use an equal-weight standardized observed score.
    data[[score_name]] <- if (length(aliases) == 1L) {
        data[[aliases]]
    } else {
        as.numeric(scale(rowMeans(data[aliases])))
    }
    list(
        data = data, variable = score_name, measurement = character(),
        model_type = paste0(
            "exploratory observed equal-weight score (",
            length(aliases), " selected indicator(s))"
        ),
        n_indicators = length(aliases)
    )
}

rhs <- function(terms) paste(terms, collapse = " + ")

met_eligible <- s4_read_mediation_table(DIR_MET, "metabolite")
inf_eligible <- s4_read_mediation_table(DIR_INFLA, "inflammation")
if (!nrow(met_eligible) || !nrow(inf_eligible))
    stop("Joint SEM requires at least one eligible metabolite and inflammatory marker (positive proportion, P_proportion < 0.05).")
if (!MET_RANK_BY %in% c("P_proportion", "P_ACME")) stop("Unknown MET_RANK_BY")
if (any(!is.finite(met_eligible[[MET_RANK_BY]]))) stop("Nonfinite ranking P values")
met_eligible <- met_eligible[order(met_eligible[[MET_RANK_BY]], met_eligible$marker), ]
inf_eligible <- inf_eligible[order(inf_eligible$P_proportion, inf_eligible$marker), ]
MODEL_SIZES <- c(all_met_all_inf = nrow(met_eligible),
                 top10_met_all_inf = min(10L, nrow(met_eligible)),
                 top5_met_all_inf = min(5L, nrow(met_eligible)))
met_eligible$selection_rank <- seq_len(nrow(met_eligible))
inf_eligible$selection_rank <- seq_len(nrow(inf_eligible))
met_sel <- make_alias_map(met_eligible, "Met")
inf_sel <- make_alias_map(inf_eligible, "Inf")
selection_audit_all <- rbind(met_sel, inf_sel)
utils::write.csv(selection_audit_all, file.path(DIR_OUT, "Cox_selection_audit.csv"), row.names = FALSE)
for (spec in names(MODEL_SIZES)) {
    sel <- rbind(head(met_sel, MODEL_SIZES[[spec]]), inf_sel)
    write.csv(sel, file.path(DIR_OUT, paste0(spec, "_planned_markers.csv")), row.names = FALSE)
    message(spec, ": ", paste(sel$marker, collapse = ", "))
}
if ("--selection-only" %in% commandArgs(trailingOnly = TRUE)) return(invisible(NULL))

cfg <- readRDS(FILE_CONFIG)
if (!OUTCOME %in% cfg$Ys ||
    !identical(cfg$primary_development_outcome, OUTCOME) ||
    !identical(cfg$primary_model, "M2")) {
    stop("analysis_config.rds does not define cvd_stroke_i/M2 as primary.")
}
m2_covariates <- s4_get_m2_covariates(cfg)
dat <- s4_validate_ids(as.data.frame(readRDS(FILE_COHORT)), "dat")
dose <- s4_validate_dose(readRDS(FILE_SCORE))$data
met <- s4_validate_ids(as.data.frame(readRDS(FILE_MET)), "met")
inflam <- s4_validate_ids(as.data.frame(readRDS(FILE_INFLA)), "inflam")

if (!setequal(dat$eid, dose$eid)) stop("Score/cohort ID mismatch.")
dat$dose <- dose$dose[match(dat$eid, dose$eid)]
if (!all(c(EVENT_VAR, m2_covariates) %in% names(dat))) {
    stop("Cohort lacks: ",
         paste(setdiff(c(EVENT_VAR, m2_covariates), names(dat)), collapse = ", "))
}
if (!all(met_sel$marker %in% names(met))) {
    stop("met.rds lacks selected metabolites: ",
         paste(setdiff(met_sel$marker, names(met)), collapse = ", "))
}
if (!all(inf_sel$marker %in% names(inflam))) {
    stop("inflammation_9.rds lacks selected markers: ",
         paste(setdiff(inf_sel$marker, names(inflam)), collapse = ", "))
}

sem_raw <- dat %>%
    dplyr::select(eid, dose, dplyr::all_of(EVENT_VAR),
                  dplyr::all_of(m2_covariates)) %>%
    dplyr::inner_join(
        met %>% dplyr::select(eid, dplyr::all_of(met_sel$marker)),
        by = "eid"
    ) %>%
    dplyr::inner_join(
        inflam %>% dplyr::select(eid, dplyr::all_of(inf_sel$marker)),
        by = "eid"
    )
rm(dat, met, inflam, dose)
invisible(gc())

analysis_vars <- unique(c(
    "dose", EVENT_VAR, m2_covariates, met_sel$marker, inf_sel$marker
))
complete <- stats::complete.cases(sem_raw[analysis_vars])
for (v in analysis_vars) {
    if (is.numeric(sem_raw[[v]]) || is.integer(sem_raw[[v]])) {
        complete <- complete & is.finite(sem_raw[[v]])
    }
}
sem_raw <- droplevels(sem_raw[complete, , drop = FALSE])
if (!nrow(sem_raw)) stop("No complete cases for the joint SEM.")
if (!all(sem_raw[[EVENT_VAR]] %in% 0:1)) stop(EVENT_VAR, " must be coded 0/1.")
if (sum(sem_raw[[EVENT_VAR]] == 1L) < 20L) stop("Too few IS events for SEM.")

sem_data <- data.frame(
    DOSE = s4_transform_direct_z(sem_raw$dose, "dose"),
    IS = as.integer(sem_raw[[EVENT_VAR]]),
    check.names = FALSE
)
for (i in seq_len(nrow(met_sel))) {
    sem_data[[met_sel$alias[i]]] <- s4_transform_log1p_z(
        sem_raw[[met_sel$marker[i]]], met_sel$marker[i]
    )
}
for (i in seq_len(nrow(inf_sel))) {
    sem_data[[inf_sel$alias[i]]] <- s4_transform_direct_z(
        sem_raw[[inf_sel$marker[i]]], inf_sel$marker[i]
    )
}

cov_design <- make_covariates(sem_raw, m2_covariates)
sem_data <- cbind(sem_data, cov_design$data)

common_sem_data <- sem_data
full_inf_sel <- inf_sel
full_met_sel <- met_sel
output_base <- DIR_OUT
fit_and_export <- function(spec, n_met, output_dir) {
MEASUREMENT_SPEC <- spec
DIR_OUT <- output_dir
dir.create(DIR_OUT, recursive = TRUE, showWarnings = FALSE)
if (file.exists(file.path(DIR_OUT, "RUN_ERROR.txt"))) unlink(file.path(DIR_OUT, "RUN_ERROR.txt"))
inf_sel <- full_inf_sel
met_sel <- head(full_met_sel, n_met)
# Both models use the exact same complete-case cohort and transformations.
sem_data <- common_sem_data[, c("DOSE", "IS", met_sel$alias, inf_sel$alias, cov_design$map$alias)]
met_domain <- make_domain(
    sem_data, met_sel$alias, "Metabolites", "Metabolite_score",
    MIN_LATENT_INDICATORS
)
sem_data <- met_domain$data
inf_domain <- make_domain(
    sem_data, inf_sel$alias, "Inflammation", "Inflammation_score",
    MIN_LATENT_INDICATORS
)
sem_data <- inf_domain$data

met_var <- met_domain$variable
inf_var <- inf_domain$variable
cov_alias <- cov_design$map$alias
measurement_lines <- c(met_domain$measurement, inf_domain$measurement)

# MODERN/SSB-style serial-parallel structural model:
# DOSE -> Metabolites; DOSE + Metabolites -> Inflammation;
# DOSE + Metabolites + Inflammation -> binary IS.
model_lines <- c(
    "# Measurement model (automatically generated)",
    measurement_lines,
    "",
    "# Structural model (M2 covariates enter every endogenous regression)",
    paste0(met_var, " ~ a*DOSE + ", rhs(cov_alias)),
    paste0(inf_var, " ~ d*", met_var, " + e*DOSE + ", rhs(cov_alias)),
    paste0("IS ~ cprime*DOSE + b_met*", met_var,
           " + b_inf*", inf_var, " + ", rhs(cov_alias)),
    "",
    "# Defined effects on the lavaan/probit coefficient scale",
    "ind_met_direct := a*b_met",
    "ind_inf_direct := e*b_inf",
    "ind_serial_met_inf := a*d*b_inf",
    "total_indirect := ind_met_direct + ind_inf_direct + ind_serial_met_inf",
    "total_effect := cprime + total_indirect"
)
model_sem <- paste(model_lines, collapse = "\n")
writeLines(model_sem, file.path(DIR_OUT, "01_SEM_model_syntax.txt"))

sem_cache <- dose_result_cache(file.path(DIR_OUT, ".checkpoints/fit"),
    list(sem_data, model_sem), packages="lavaan")
fitted <- sem_cache("sem", {
warnings_sem <- character()
fit_sem <- withCallingHandlers(
    lavaan::sem(
        model = model_sem,
        data = sem_data,
        ordered = "IS",
        estimator = "WLSMV",
        optim.method = "BFGS",
        parameterization = "theta",
        std.lv = TRUE,
        fixed.x = TRUE,
        meanstructure = TRUE,
        control = list(iter.max = 20000L)
    ),
    warning = function(w) {
        warnings_sem <<- c(warnings_sem, conditionMessage(w))
        invokeRestart("muffleWarning")
    }
)
list(fit=fit_sem, warnings=warnings_sem)
})
fit_sem <- fitted$fit
warnings_sem <- fitted$warnings
rm(fitted)

converged <- isTRUE(lavaan::lavInspect(fit_sem, "converged"))
if (!converged) {
    dose_save_rds(fit_sem, file.path(DIR_OUT, "SEM_fit.rds"))
    writeLines(c(warnings_sem, "SEM did not converge"), file.path(DIR_OUT, "09_SEM_warnings.txt"))
    stop("SEM did not converge; diagnostic object saved.")
}

parameters <- lavaan::parameterEstimates(
    fit_sem, standardized = TRUE, ci = TRUE
)
parameters$FDR_within_operator <- NA_real_
for (operator in intersect(c("=~", "~", ":="), unique(parameters$op))) {
    idx <- which(parameters$op == operator & is.finite(parameters$pvalue))
    parameters$FDR_within_operator[idx] <- stats::p.adjust(
        parameters$pvalue[idx], method = "BH"
    )
}
structural <- parameters %>%
    dplyr::filter(op == "~", lhs %in% c(met_var, inf_var, "IS"))
indirect <- parameters %>% dplyr::filter(op == ":=")
loadings <- parameters %>% dplyr::filter(op == "=~")

# Explicit measurement-model QC. A standardized loading above 1 or a negative
# indicator residual variance is a Heywood/improper solution and invalidates a
# confirmatory interpretation even when the optimizer reports convergence.
indicator_residuals <- parameters %>%
    dplyr::filter(
        op == "~~", lhs == rhs,
        lhs %in% c(met_sel$alias, inf_sel$alias)
    ) %>%
    dplyr::transmute(
        alias = lhs, residual_variance = est,
        residual_variance_std = std.all
    )
indicator_qc <- loadings %>%
    dplyr::transmute(
        latent = lhs, alias = rhs,
        loading = est, standardized_loading = std.all
    ) %>%
    dplyr::left_join(indicator_residuals, by = "alias") %>%
    dplyr::left_join(
        dplyr::bind_rows(
            met_sel %>% dplyr::select(alias, marker, domain),
            inf_sel %>% dplyr::select(alias, marker, domain)
        ),
        by = "alias"
    ) %>%
    dplyr::mutate(
        loading_gt_one = abs(standardized_loading) >
            1 + HEYWOOD_TOLERANCE,
        negative_residual_variance = residual_variance <
            -HEYWOOD_TOLERANCE,
        improper = loading_gt_one | negative_residual_variance
    ) %>%
    dplyr::arrange(domain, dplyr::desc(abs(standardized_loading)))
post_check <- isTRUE(lavaan::lavInspect(fit_sem, "post.check"))
improper_solution <- !post_check || any(indicator_qc$improper, na.rm = TRUE)
if (improper_solution) warnings_sem <- c(warnings_sem, "Improper solution: see indicator QC and lavaan post.check.")
if (improper_solution) {
    bad <- indicator_qc$marker[indicator_qc$improper]
    warning(
        "Improper SEM solution for indicator(s): ",
        paste(bad, collapse = ", "),
        ". Do not interpret the SEM as a valid confirmatory model."
    )
    if (STOP_ON_IMPROPER_SOLUTION) {
        stop("STOP_ON_IMPROPER_SOLUTION=TRUE.")
    }
}

fit_names <- c("chisq.scaled", "df.scaled", "pvalue.scaled", "cfi.scaled",
               "tli.scaled", "rmsea.scaled", "rmsea.ci.lower.scaled",
               "rmsea.ci.upper.scaled", "srmr")
fit_measures <- tryCatch(
    data.frame(
        measure = fit_names,
        value = as.numeric(lavaan::fitMeasures(fit_sem, fit_names))
    ),
    error = function(e) data.frame(
        measure = "fitMeasures_error", value = NA_real_,
        message = conditionMessage(e)
    )
)
r2 <- tryCatch(
    data.frame(
        variable = names(lavaan::inspect(fit_sem, "r2")),
        R2 = as.numeric(lavaan::inspect(fit_sem, "r2"))
    ),
    error = function(e) data.frame(variable = character(), R2 = numeric())
)
gradient <- tryCatch(
    max(abs(lavaan::lavInspect(fit_sem, "optim.gradient")), na.rm = TRUE),
    error = function(e) NA_real_
)

model_qc <- data.frame(
    N = nrow(sem_data), IS_events = sum(sem_data$IS),
    IS_event_percent = 100 * mean(sem_data$IS),
    N_metabolites_eligible = nrow(met_eligible),
    N_inflammation_eligible = nrow(inf_eligible),
    N_metabolites_selected = nrow(met_sel),
    N_inflammation_selected = nrow(inf_sel),
    metabolite_model = met_domain$model_type,
    inflammation_model = inf_domain$model_type,
    selection_source = paste(FILE_MET_COX, FILE_INFLA_COX, sep="; "),
    metabolite_ranking = MET_RANK_BY,
    selection_rule = SELECTION_RULE,
    measurement_spec = MEASUREMENT_SPEC,
    estimator = "WLSMV", optimizer = "BFGS", ordered_outcome = "IS",
    converged = converged, improper_solution = improper_solution,
    N_improper_indicators = sum(indicator_qc$improper, na.rm = TRUE),
    max_abs_standardized_loading = max(
        abs(indicator_qc$standardized_loading), na.rm = TRUE
    ),
    min_indicator_residual_variance = min(
        indicator_qc$residual_variance, na.rm = TRUE
    ),
    max_abs_gradient = gradient,
    stringsAsFactors = FALSE
)

marker_map <- dplyr::bind_rows(met_sel, inf_sel)
selection_audit <- selection_audit_all
selection_audit$retained_in_measurement_model <- selection_audit$marker %in% marker_map$marker
selection_audit$measurement_spec <- MEASUREMENT_SPEC
utils::write.csv(marker_map,
                 file.path(DIR_OUT, "00_selected_marker_alias_map.csv"),
                 row.names = FALSE)
utils::write.csv(
    selection_audit,
    file.path(DIR_OUT, "00a_Cox_selection_measurement_audit.csv"),
    row.names = FALSE
)
utils::write.csv(cov_design$map,
                 file.path(DIR_OUT, "00_M2_covariate_design_map.csv"),
                 row.names = FALSE)
utils::write.csv(model_qc, file.path(DIR_OUT, "02_SEM_sample_model_QC.csv"),
                 row.names = FALSE)
utils::write.csv(fit_measures, file.path(DIR_OUT, "03_SEM_fit_measures.csv"),
                 row.names = FALSE)
utils::write.csv(parameters, file.path(DIR_OUT, "04_SEM_all_parameters.csv"),
                 row.names = FALSE)
utils::write.csv(structural, file.path(DIR_OUT, "05_SEM_structural_paths.csv"),
                 row.names = FALSE)
utils::write.csv(indirect, file.path(DIR_OUT, "06_SEM_indirect_effects.csv"),
                 row.names = FALSE)
utils::write.csv(loadings, file.path(DIR_OUT, "07_SEM_factor_loadings.csv"),
                 row.names = FALSE)
utils::write.csv(r2, file.path(DIR_OUT, "08_SEM_R2.csv"), row.names = FALSE)
utils::write.csv(
    indicator_qc, file.path(DIR_OUT, "08a_SEM_indicator_QC.csv"),
    row.names = FALSE
)
writeLines(unique(warnings_sem), file.path(DIR_OUT, "09_SEM_warnings.txt"))
dose_save_rds(fit_sem, file.path(DIR_OUT, "SEM_fit.rds"))

notes <- c(
    paste0("Membership from upstream Cox mediation CSVs: ", FILE_MET_COX, "; ", FILE_INFLA_COX),
    paste0("Eligible markers: ", nrow(met_eligible), " metabolites and ", nrow(inf_eligible), " inflammatory markers."),
    paste0("Metabolite ranking: ", MET_RANK_BY, " ascending; ties by marker ID."),
    "Eligibility: upstream candidate, positive proportion, unadjusted P_proportion < 0.05. No additional FDR filter.",
    "Three nested models: all / top 10 / top 5 eligible metabolites, each with all eligible inflammatory markers. All use the same complete-case cohort based on the union of eligible indicators.",
    paste0("Measurement specification: ", MEASUREMENT_SPEC),
    "Metabolites log1p then Z; inflammation directly Z; DOSE Z; M2 covariates in every structural regression.",
    "WLSMV with ordered binary IS: effects are probit coefficients, not Cox hazard ratios; follow-up time is not modeled.",
    "Same-data selection is exploratory, not independent or causal validation. Inspect fit and improper solutions."
)
writeLines(notes, file.path(DIR_OUT, "README_SEM.txt"))

workbook <- list(
    "Model QC" = model_qc,
    "Fit measures" = fit_measures,
    "Selected markers" = marker_map,
    "Cox selection audit" = selection_audit,
    "M2 covariate map" = cov_design$map,
    "Structural paths" = structural,
    "Indirect effects" = indirect,
    "Factor loadings" = loadings,
    "R2" = r2,
    "Indicator QC" = indicator_qc,
    "All parameters" = parameters,
    "Notes" = data.frame(note = notes)
)
openxlsx::write.xlsx(
    workbook,
    file.path(DIR_OUT, "DOSE_metabolite_inflammation_IS_SEM.xlsx"),
    overwrite = TRUE, asTable = TRUE
)

if (DRAW_SEM_PLOT && converged &&
    nrow(marker_map) <= MAX_PLOT_INDICATORS &&
    requireNamespace("semPlot", quietly = TRUE)) {
    grDevices::pdf(file.path(DIR_OUT, "Fig_SEM_paths.pdf"),
                   width = 14, height = 9, onefile = TRUE)
    semPlot::semPaths(
        fit_sem, what = "std", whatLabels = "std", layout = "tree2",
        residuals = FALSE, intercepts = FALSE, thresholds = FALSE,
        nCharNodes = 0, edge.label.cex = 0.7, sizeMan = 4, sizeLat = 7
    )
    bounded_dev_off()
}

s4_write_manifest(input_files, DIR_OUT)
capture.output(lavaan::summary(fit_sem, standardized = TRUE, fit.measures = TRUE,
                       rsquare = TRUE),
               file = file.path(DIR_OUT, "SEM_summary.txt"))

print(model_qc)
print(fit_measures)
print(indirect)
message("Completed joint SEM. Results: ", DIR_OUT)

invisible(list(qc = model_qc, fit = fit_measures))
}
comparison <- list()
for (spec in names(MODEL_SIZES)) {
    comparison[[spec]] <- tryCatch({
        res <- fit_and_export(spec, MODEL_SIZES[[spec]], file.path(output_base, spec))
        row <- cbind(model = spec, status = "completed", res$qc)
        for (j in seq_len(nrow(res$fit))) row[[res$fit$measure[j]]] <- res$fit$value[j]
        row
    }, error = function(e) {
        dir.create(file.path(output_base, spec), recursive = TRUE, showWarnings = FALSE)
        writeLines(conditionMessage(e), file.path(output_base, spec, "RUN_ERROR.txt"))
        warning(spec, ": ", conditionMessage(e))
        data.frame(model = spec, status = "failed", error = conditionMessage(e))
    })
    write.csv(dplyr::bind_rows(comparison), file.path(output_base, "three_model_comparison.csv"), row.names = FALSE)
}
}
run_cox_sem()
}))

# ===========================================================================
# mechanism / sem_plot
# ===========================================================================
dose_register_step("mechanism", "sem_plot", expression({
# Display-only selection from the existing SEM. No refitting or raw data needed.
# All arrows and loadings are standardized (std.all); P values are the
# original SEM parameter tests. Selection is NOT mediation ranking.
args <- commandArgs(trailingOnly = TRUE)
ROOT <- dose_output_dir()
DEFAULT_ROOT <- file.path(ROOT, 'ukb/05mechanism/sem/M2_Cox')
INPUT_DIRS <- if (length(args)) args else file.path(DEFAULT_ROOT,
    c('all_met_all_inf', 'top10_met_all_inf', 'top5_met_all_inf'))
plot_one <- function(DIR_IN) {
DIR_OUT <- file.path(DIR_IN, 'plots')
TOP_MET <- 3L
TOP_INF <- 3L
dir.create(DIR_OUT, recursive = TRUE, showWarnings = FALSE)
read_result <- function(name) read.csv(file.path(DIR_IN, name), check.names = FALSE)
mapping <- read_result('00_selected_marker_alias_map.csv')
loading <- read_result('07_SEM_factor_loadings.csv')
paths <- read_result('05_SEM_structural_paths.csv')
effects <- read_result('06_SEM_indirect_effects.csv')
parameters <- read_result('04_SEM_all_parameters.csv')
qc <- read_result('02_SEM_sample_model_QC.csv')
fm <- read_result('03_SEM_fit_measures.csv')
stopifnot(isTRUE(qc$converged[1]))
loading$marker <- mapping$marker[match(loading$rhs, mapping$alias)]
stopifnot(!anyNA(loading$marker))
rv <- parameters[parameters$op == '~~' & parameters$lhs == parameters$rhs, ]
loading$residual_variance <- rv$est[match(loading$rhs, rv$lhs)]
loading$improper <- abs(loading$std.all) > 1 |
    (!is.na(loading$residual_variance) & loading$residual_variance < 0)
loading$abs_loading_rank <- ave(abs(loading$std.all), loading$lhs,
    FUN = function(x) rank(-x, ties.method = 'first'))
# Compact reference-style figure shows three strongest absolute loadings per domain.
loading$selected_for_display <- (loading$lhs == 'Metabolites' & loading$abs_loading_rank <= TOP_MET) |
    (loading$lhs == 'Inflammation' & loading$abs_loading_rank <= TOP_INF)
loading$source_display_name <- mapping$display_name[match(loading$rhs, mapping$alias)]
labels <- c(bc_NEUT = 'Neutrophils', bb_CRP = 'CRP', SII = 'SII',
    NLR = 'NLR', bc_PLT = 'Platelets', bc_MONO = 'Monocytes', bc_LYMPH = 'Lymphocytes')
pretty_name <- function(x) {
    out <- gsub('_', '-', sub('^met_', '', x))
    out <- sub('\\.pct$', ' %', out)
    hit <- x %in% names(labels)
    out[hit] <- labels[x[hit]]
    out
}
loading$display_name <- pretty_name(loading$marker)
loading <- loading[order(loading$lhs, loading$abs_loading_rank), ]
selected <- loading[loading$selected_for_display, ]
write.csv(loading, file.path(DIR_OUT, 'SEM_indicator_selection_audit.csv'), row.names = FALSE)
write.csv(selected, file.path(DIR_OUT, 'SEM_main_indicators.csv'), row.names = FALSE)
write.csv(paths[nzchar(paths$label), ], file.path(DIR_OUT, 'SEM_plotted_paths.csv'), row.names = FALSE)
write.csv(effects, file.path(DIR_OUT, 'SEM_indirect_effects.csv'), row.names = FALSE)

library(grid)
ink <- '#202426'
cols <- c(dose = '#C6DBE7', met = '#CBC9E2', inf = '#F6EBD5', is = '#F0CDCE')
borders <- c(dose = '#6F9FB6', met = '#9490C4', inf = '#E7C88F', is = '#C18D90')
txt <- function(label, x, y, size = 11, col = ink, face = 'plain', just = 'centre', rot = 0)
    grid.text(label, x, y, just = just, rot = rot,
        gp = gpar(fontsize = size, col = col, fontface = face, fontfamily = 'sans'))
node <- function(label, x, y, key, width = .17) {
    grid.roundrect(x, y, width = width, height = .10, r = unit(3.5, 'mm'),
        gp = gpar(fill = cols[key], col = borders[key], lwd = 1.5))
    txt(label, x, y, 15)
}
# A segmented Bezier path provides a portable vector gradient (no raster assets).
# Coordinates below use a 12 x 4.8 inch device.
grad_arrow <- function(p0, p1, p2, p3, color, width = 5, dashed = FALSE) {
    t <- seq(0, 1, length.out = 100)
    xy <- outer((1-t)^3,p0) + outer(3*(1-t)^2*t,p1) +
          outer(3*(1-t)*t^2,p2) + outer(t^3,p3)
    shades <- colorRampPalette(c('#F3F6FA', color))(99)
    head_length <- if (width >= 9) .21 else .14
    head_halfwidth <- if (width >= 9) .13 else .075
    distance_to_tip <- sqrt(rowSums(sweep(xy,2,p3)^2 *
        matrix(rep(c(12,4.8)^2,each=nrow(xy)),ncol=2)))
    for (k in 1:99) {
        if (distance_to_tip[k+1] >= head_length*.85 && (!dashed || (k %% 10) < 6))
            grid.lines(xy[k:(k+1),1],xy[k:(k+1),2],
                gp = gpar(col = shades[k], lwd = width, lineend = 'round'))
    }
    # Arrow head computed in physical coordinates to maintain its shape.
    tip <- p3*c(12,4.8)
    delta <- (p3-xy[97,])*c(12,4.8)
    u <- delta/sqrt(sum(delta^2)); n <- c(-u[2],u[1])
    vertices <- rbind(tip, tip-head_length*u+head_halfwidth*n, tip-head_length*u-head_halfwidth*n)
    grid.polygon(vertices[,1]/12,vertices[,2]/4.8,gp=gpar(fill=color,col=NA))
}
path_label <- function(key,x,y,rot=0) {
    r <- paths[paths$label == key, ]
    stopifnot(nrow(r)==1,is.finite(r$std.all),is.finite(r$pvalue))
    sig <- if(r$pvalue<.001) '***' else if(r$pvalue<.01) '**' else if(r$pvalue<.05) '*' else 'ns'
    label <- bquote(.(sprintf('%.3f',r$std.all))^.(sig))
    txt(label,x,y,14,rot=rot)
}
table_panel <- function(domain,title,top,fill) {
    d <- selected[selected$lhs==domain, ]
    grid.rect(.795,top-.049,.35,.098,gp=gpar(fill=fill,col=NA))
    txt(title,.795,top-.049,14,col='white',face='bold')
    for(i in seq_len(nrow(d))) {
        yy <- top-.098-(i-.5)*.089
        grid.rect(.795,yy,.35,.089,gp=gpar(fill=if(i%%2) '#F3F3F3' else '#E0E0E0',col=NA))
        txt(d$display_name[i],.747,yy,12.5)
        txt(sprintf('%.3f%s',d$std.all[i],if(d$improper[i]) ' !' else ''),.903,yy,12.5,
            col=if(d$improper[i]) '#A53B39' else ink)
    }
}
measure <- function(key) fm$value[match(key,fm$measure)]
draw <- function() {
    grid.newpage()
    txt('D.',.028,.935,23,just='left')
    txt('E.',.574,.935,23,just='left')
    grad_arrow(c(.205,.85),c(.145,.85),c(.105,.72),c(.105,.565),'#E7C88F',4)
    grad_arrow(c(.375,.85),c(.425,.85),c(.475,.72),c(.475,.565),'#AAA5D0',4)
    grad_arrow(c(.105,.455),c(.105,.30),c(.150,.19),c(.200,.19),'#D8AB9E',4)
    grad_arrow(c(.475,.455),c(.475,.30),c(.425,.19),c(.380,.19),'#AA96B7',4,
        dashed=paths$pvalue[paths$label=='b_met']>=.05)
    # Draw horizontal path first so central vertical path passes above it.
    grad_arrow(c(.378,.51),c(.325,.51),c(.245,.51),c(.202,.51),'#71A8AF',9)
    grad_arrow(c(.290,.785),c(.290,.63),c(.290,.41),c(.290,.24),'#D39A9E',10)
    path_label('e',.104,.735,55)
    path_label('a',.476,.735,-55)
    path_label('b_inf',.110,.315,-55)
    path_label('b_met',.471,.315,55)
    path_label('d',.243,.574)
    path_label('cprime',.311,.387,-90)
    node('DOSE',.290,.85,'dose')
    node('Inflammation',.105,.51,'inf',.185)
    node('Metabolites',.475,.51,'met',.185)
    node('IS',.290,.19,'is',.18)
    table_panel('Metabolites','Metabolites',.90,'#9994C7')
    table_panel('Inflammation','Inflammation',.485,'#D49698')
    txt('Standardized coefficients and loadings.  * P < 0.05; ** P < 0.01; *** P < 0.001; ns, not significant.',
        .035,.043,8.5,col='#586069',just='left')
}
stem <- file.path(DIR_OUT, 'SEM_main_indicators_path')
bounded_png(paste0(stem, '.png'), width = 3600, height = 1440, res = 300, type = 'cairo')
draw(); bounded_dev_off()
pdf(paste0(stem, '.pdf'), width = 12, height = 4.8, useDingbats = FALSE)
draw(); bounded_dev_off()
svg(paste0(stem, '.svg'), width = 12, height = 4.8)
draw(); bounded_dev_off()
# Complete loadings as a separate multi-page vector figure, avoiding 23-row crowding.
pdf(file.path(DIR_OUT, 'SEM_all_indicator_loadings.pdf'), width = 9, height = 10, useDingbats = FALSE)
for (domain in unique(loading$lhs)) {
    dd <- loading[loading$lhs == domain, ]
    for (start in seq(1L, nrow(dd), by = 18L)) {
        page <- dd[start:min(start+17L, nrow(dd)), ]
        grid.newpage()
        txt(paste(basename(DIR_IN), domain), .06, .95, 16, face = 'bold', just = 'left')
        txt('Indicator', .07, .885, 12, face = 'bold', just = 'left')
        txt('Standardized loading', .92, .885, 12, face = 'bold', just = 'right')
        for (j in seq_len(nrow(page))) {
            yy <- .845 - (j-1)*.040
            txt(page$display_name[j], .07, yy, 11, just = 'left')
            txt(sprintf('%.3f%s', page$std.all[j], if (page$improper[j]) ' *' else ''),
                .92, yy, 11, col = if (page$improper[j]) '#B14235' else ink, just = 'right')
        }
        txt('* Negative residual variance or absolute standardized loading > 1.', .07, .07, 9, just = 'left')
    }
}
bounded_dev_off()
writeLines(c(
    'Figure shows top three absolute standardized loadings per domain, preserving their signs; all fitted indicators are retained in the model.',
    sprintf('Full model includes %d metabolic and %d inflammatory indicators.', sum(loading$lhs == 'Metabolites'), sum(loading$lhs == 'Inflammation')),
    'This is not a ranking of individual mediation effects or biological importance.',
    'Negative residual variances / absolute standardized loadings >1 are flagged, not silently excluded.',
    'Arrows and table values use std.all. P values come from original unstandardized parameter tests.',
    'Inspect the fit indices and QC flags; coefficients alone cannot establish mediation.',
    'No model was refitted; all inputs are existing result CSV files.'
), file.path(DIR_OUT, 'README_plot.txt'))
writeLines(c(
    sprintf('Exploratory SEM: %d metabolic and %d inflammatory indicators; N=%d, IS events=%d.',
        qc$N_metabolites_selected,qc$N_inflammation_selected,qc$N,qc$IS_events),
    'Panel D: standardized structural path coefficients. Stars refer to original parameter P values, not FDR-adjusted values.',
    'Panel E: three indicators with the largest absolute standardized factor loadings in each domain; signed loadings shown.',
    'This display selection does not refit or reduce the model. The complete loadings are in SEM_all_indicator_loadings.pdf.',
    sprintf('CFI=%.3f; TLI=%.3f; RMSEA=%.3f; SRMR=%.3f.',
        measure('cfi.scaled'),measure('tli.scaled'),measure('rmsea.scaled'),measure('srmr')),
    'Overall fit remains limited; these exploratory paths do not establish causal mediation.',
    if(any(loading$improper)) '! indicates an improper indicator estimate; see the indicator selection audit.' else 'No indicator Heywood flag detected.'
),file.path(DIR_OUT,'SEM_figure_caption.txt'))
print(selected[, c('marker', 'std.all', 'residual_variance', 'improper')], row.names = FALSE)
message('Plots and selection audit saved to: ', normalizePath(DIR_OUT))

}
for (DIR_IN in INPUT_DIRS) {
    if (!file.exists(file.path(DIR_IN, '02_SEM_sample_model_QC.csv')) ||
        file.exists(file.path(DIR_IN, 'RUN_ERROR.txt'))) {
        warning('No completed model to plot: ', DIR_IN)
        next
    }
    plot_one(DIR_IN)
}
}))

# ===========================================================================
# mechanism / sem_sensitivity
# ===========================================================================
dose_register_step("mechanism", "sem_sensitivity", expression({
# Exploratory sensitivity analysis; do not select a model by mediation P value.
# Twelve models: four selection rules x three measurement specifications.
# No existing SEM outputs are overwritten. Each invocation has its own directory.
ROOT <- dose_output_dir()
UKB <- dose_ukb_dir()
args <- commandArgs(trailingOnly = TRUE)
PLAN_ONLY <- '--plan-only' %in% args
COR_CUTOFF <- 0.90
SEED <- 12345L
set.seed(SEED)
BASE <- file.path(ROOT, 'ukb/05mechanism')
OUT <- tempfile(paste0(format(Sys.time(), '%Y%m%d_%H%M%S'), '_'),
                tmpdir = file.path(BASE, 'sem/sensitivity'))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
read_domain <- function(domain, filename) {
    z <- read.csv(file.path(BASE, domain, 'M2', filename), check.names = FALSE)
    z$marker <- z[[domain]]
    stopifnot(!anyDuplicated(z$marker))
    z
}
select_markers <- function(domain, rule) {
    if (rule == 'xm_top10') {
        z <- read_domain(domain, paste0('01a_DOSE_', domain, '_XM_all.csv'))
        z <- z[is.finite(z$FDR_XM) & z$FDR_XM < .05, ]
        z <- z[order(z$FDR_XM, z$P_XM, z$marker), ]
        if (domain == 'metabolite') z <- head(z, 10)
    } else {
        fn <- if (rule == 'aft_positive') '05_mediation_LHR_AFT.csv' else '07_mediation_Cox_product.csv'
        z <- read_domain(domain, fn)
        keep <- if (rule == 'aft_positive') {
            is.finite(z$FDR_ACME) & z$FDR_ACME < .05 &
                is.finite(z$FDR_proportion) & z$FDR_proportion < .05
        } else is.finite(z$P_proportion) & z$P_proportion < .05
        if (rule != 'cox_all') keep <- keep & is.finite(z$Proportion_mediated) & z$Proportion_mediated > 0
        z <- z[keep, ]
    }
    # Alphabetical order fixes pruning priority independently of effect/P values.
    sort(z$marker)
}
rules <- c('cox_all', 'cox_positive', 'aft_positive', 'xm_top10')
sets <- setNames(lapply(rules, function(r) list(
    met = select_markers('metabolite', r), inf = select_markers('inflammation', r))), rules)
plan <- expand.grid(rule = rules, measurement = c('latent', 'pruned', 'pca'), stringsAsFactors = FALSE)
plan$id <- paste(plan$rule, plan$measurement, sep = '_')
plan$n_met_selected <- vapply(plan$rule, function(r) length(sets[[r]]$met), integer(1))
plan$n_inf_selected <- vapply(plan$rule, function(r) length(sets[[r]]$inf), integer(1))
write.csv(plan, file.path(OUT, 'analysis_plan.csv'), row.names = FALSE)
for (r in rules) write.csv(data.frame(marker = c(sets[[r]]$met, sets[[r]]$inf),
    domain = c(rep('metabolite', length(sets[[r]]$met)), rep('inflammation', length(sets[[r]]$inf)))),
    file.path(OUT, paste0(r, '_selected.csv')), row.names = FALSE)
writeLines(c(
    'All models are exploratory; report all attempted models, including failures.',
    'cox_all: Cox P_proportion<0.05, either sign; cox_positive also requires proportion>0.',
    'aft_positive: ACME and proportion BH-FDR<0.05, proportion>0.',
    'xm_top10: top 10 metabolites by X-M FDR/P and all X-M FDR<0.05 inflammatory markers.',
    'Cox/AFT candidate files retain upstream X-M FDR and M-Y P selection.',
    'latent: full selected indicators, one factor per domain.',
    'pruned: inflammation restricted to selected raw NEUT/CRP/PLT/MONO/LYMPH; excludes derived ratios.',
    'pruned: sequentially discard indicators with absolute Pearson correlation>=0.90 to a retained indicator.',
    'Pruning priority is alphabetical, not based on mediation significance. All exclusions are recorded.',
    'pca: one PC1 per domain; variance explained and signed weights exported. No latent measurement model.',
    'A single complete-case cohort for the union of ALL planned indicators is used for every model.',
    'Metabolites log1p then Z; inflammation directly Z; DOSE Z; M2 covariates in all regressions.',
    'IS is binary WLSMV/probit, not a Cox or AFT survival SEM. Follow-up time is not modeled.',
    'Quality flags require convergence, small gradient, finite SE, post.check, no negative variance, no |loading|>1.',
    'CFI>=0.90,TLI>=0.90,RMSEA<=0.08,SRMR<=0.08 are descriptive screening heuristics, not universal proof of fit.',
    'Do not compare AIC/BIC or declare a winner across different observed-variable sets.',
    'PCA path models are saturated/just identified here: fit indices cannot demonstrate superiority.',
    'Same-data mediator selection makes subsequent inference exploratory; this is not independent validation.',
    'PCA scores are treated as fixed observed variables; SE do not include score-estimation or selection uncertainty.',
    'Effect stability across specifications is descriptive because latent factors/scores have different meanings.',
    'No automated ranking by P values or winner selection is performed.'
), file.path(OUT, 'README.txt'))
print(plan, row.names = FALSE)
message('Output: ', OUT)

main <- function() {
    if (PLAN_ONLY) return(invisible(NULL))
    if (!requireNamespace('lavaan', quietly = TRUE)) stop('Install lavaan before running models.')
    only <- sub('^--only=', '', args[grepl('^--only=', args)])
    run_plan <- plan
    if (length(only)) {
        stopifnot(length(only) == 1, only %in% plan$id)
        run_plan <- plan[plan$id == only, ]
    }
    cfg <- readRDS(file.path(ROOT, 'ukb/00data/analysis_config.rds'))
    covs <- s4_get_m2_covariates(cfg)
    ev <- 'cvd_stroke_i.Yt2e'
    dat <- s4_validate_ids(as.data.frame(readRDS(file.path(ROOT, 'ukb/00data/dat.base.imp.rds'))), 'cohort')
    dose <- s4_validate_dose(readRDS(file.path(ROOT, 'ukb/04score/dose_score.rds')))$data
    stopifnot(setequal(dat$eid, dose$eid))
    dat$dose <- dose$dose[match(dat$eid, dose$eid)]
    met_names <- sort(unique(unlist(lapply(sets, `[[`, 'met'))))
    inf_names <- sort(unique(unlist(lapply(sets, `[[`, 'inf'))))
    met <- s4_validate_ids(as.data.frame(readRDS(file.path(UKB, 'Rdata/met.rds'))), 'met')
    inf <- s4_validate_ids(as.data.frame(readRDS(file.path(BASE, 'inflammation/M2/inflammation_9.rds'))), 'inf')
    dat <- merge(dat[, c('eid', 'dose', ev, covs)], met[, c('eid', met_names)], by = 'eid')
    dat <- merge(dat, inf[, c('eid', inf_names)], by = 'eid')
    rm(met, inf, dose); gc()
    keep <- complete.cases(dat)
    for (nm in names(dat)) if (is.numeric(dat[[nm]])) keep <- keep & is.finite(dat[[nm]])
    dat <- droplevels(dat[keep, ])
    stopifnot(nrow(dat) > 0, all(dat[[ev]] %in% 0:1), sum(dat[[ev]]) >= 20)
    z <- data.frame(DOSE = s4_transform_direct_z(dat$dose, 'DOSE'), IS = as.integer(dat[[ev]]))
    mp <- data.frame(marker = c(met_names, inf_names), alias = c(sprintf('Met%03d', seq_along(met_names)), sprintf('Inf%03d', seq_along(inf_names))))
    for (j in seq_len(nrow(mp))) z[[mp$alias[j]]] <- if (mp$marker[j] %in% met_names)
        s4_transform_log1p_z(dat[[mp$marker[j]]], mp$marker[j]) else s4_transform_direct_z(dat[[mp$marker[j]]], mp$marker[j])
    mm <- model.matrix(reformulate(covs), dat)
    mm <- mm[, colnames(mm) != '(Intercept)', drop = FALSE]
    mm <- mm[, apply(mm, 2, sd) > 0, drop = FALSE]
    for (j in seq_len(ncol(mm))) if (!all(mm[, j] %in% 0:1)) mm[, j] <- as.numeric(scale(mm[, j]))
    cv <- sprintf('Cov%03d', seq_len(ncol(mm)))
    write.csv(data.frame(alias = cv, original = colnames(mm)), file.path(OUT, 'covariate_map.csv'), row.names = FALSE)
    colnames(mm) <- cv; z <- cbind(z, mm)
    write.csv(mp, file.path(OUT, 'marker_alias_map.csv'), row.names = FALSE)
    write.csv(data.frame(N = nrow(z), events = sum(z$IS)), file.path(OUT, 'common_cohort.csv'), row.names = FALSE)
    rm(dat, mm); gc()
    aliases <- function(nms) mp$alias[match(nms, mp$marker)]
    prune <- function(nms) {
        retained <- character(); audit <- data.frame(marker = character(), excluded_due_to = character(), abs_r = numeric())
        for (nm in nms) {
            rr <- if (length(retained)) abs(cor(z[[aliases(nm)]], z[, aliases(retained), drop = FALSE])) else numeric()
            if (length(rr) && max(rr) >= COR_CUTOFF) audit <- rbind(audit, data.frame(marker = nm,
                excluded_due_to = retained[which.max(rr)], abs_r = max(rr))) else retained <- c(retained, nm)
        }
        list(keep = retained, audit = audit)
    }
    rows <- list(); all_effects <- list()
    for (i in seq_len(nrow(run_plan))) {
        spec <- run_plan[i, ]; folder <- file.path(OUT, spec$id); dir.create(folder)
        message('Fitting ', i, '/', nrow(run_plan), ': ', spec$id)
        warnings <- character()
        result <- tryCatch(withCallingHandlers({
            m <- sets[[spec$rule]]$met; f <- sets[[spec$rule]]$inf
            if (spec$measurement == 'pruned') {
                raw <- intersect(f, c('bc_NEUT', 'bb_CRP', 'bc_PLT', 'bc_MONO', 'bc_LYMPH'))
                writeLines(setdiff(f, raw), file.path(folder, 'excluded_derived_inflammation.txt'))
                pm <- prune(m); pf <- prune(raw); m <- pm$keep; f <- pf$keep
                write.csv(rbind(pm$audit, pf$audit), file.path(folder, 'correlation_exclusions.csv'), row.names = FALSE)
            }
            write.csv(mp[mp$marker %in% c(m, f), ], file.path(folder, 'used_markers.csv'), row.names = FALSE)
            if (!length(m) || !length(f)) stop('No selected indicator in one domain.')
            data_fit <- z[, c('DOSE', 'IS', aliases(c(m, f)), cv)]
            measurement <- character()
            if (spec$measurement == 'pca') {
                for (domain in c('Metabolites', 'Inflammation')) {
                    nm <- if (domain == 'Metabolites') m else f
                    pc <- prcomp(z[, aliases(nm), drop = FALSE], center = TRUE, scale. = FALSE)
                    # Fix sign by largest absolute weight, not outcome association.
                    sign_pc <- sign(pc$rotation[which.max(abs(pc$rotation[,1])),1])
                    data_fit[[domain]] <- as.numeric(scale(pc$x[,1] * sign_pc))
                    write.csv(data.frame(marker = nm, PC1_weight = pc$rotation[,1]*sign_pc,
                        PC1_variance_fraction = pc$sdev[1]^2/sum(pc$sdev^2)),
                        file.path(folder, paste0(domain, '_PC1.csv')), row.names = FALSE)
                }
                data_fit <- data_fit[, c('DOSE','IS','Metabolites','Inflammation',cv)]
            } else {
                if (length(m) < 3 || length(f) < 3) stop('Fewer than three indicators: latent model skipped, no automatic fallback.')
                measurement <- c(paste('Metabolites =~', paste(aliases(m), collapse=' + ')),
                    paste('Inflammation =~', paste(aliases(f), collapse=' + ')))
            }
            syntax <- paste(c(measurement,
                paste('Metabolites ~ a*DOSE +', paste(cv, collapse=' + ')),
                paste('Inflammation ~ d*Metabolites + e*DOSE +', paste(cv, collapse=' + ')),
                paste('IS ~ cprime*DOSE + b_met*Metabolites + b_inf*Inflammation +', paste(cv, collapse=' + ')),
                'ind_met := a*b_met', 'ind_inf := e*b_inf', 'ind_serial := a*d*b_inf',
                'total_indirect := a*b_met + e*b_inf + a*d*b_inf',
                'total_effect := cprime + a*b_met + e*b_inf + a*d*b_inf'), collapse='\n')
            writeLines(syntax, file.path(folder, 'model.txt'))
            fit <- lavaan::sem(syntax, data = data_fit, ordered = 'IS', estimator = 'WLSMV',
                parameterization = 'theta', std.lv = TRUE, fixed.x = TRUE,
                optim.method = 'BFGS', control = list(iter.max = 20000L))
            dose_save_rds(fit, file.path(folder, 'fit.rds'))
            converged <- isTRUE(lavaan::lavInspect(fit, 'converged'))
            if (!converged) stop('Model did not converge; fit object saved.')
            pe <- lavaan::parameterEstimates(fit, standardized = TRUE, ci = TRUE)
            write.csv(pe, file.path(folder, 'parameters.csv'), row.names = FALSE)
            ie <- pe[pe$op == ':=', ]; ie$model <- spec$id; all_effects[[spec$id]] <- ie
            write.csv(ie, file.path(folder, 'indirect_effects.csv'), row.names = FALSE)
            measures <- lavaan::fitMeasures(fit, c('df.scaled','cfi.scaled','tli.scaled','rmsea.scaled','srmr'))
            write.csv(data.frame(measure=names(measures), value=as.numeric(measures)), file.path(folder,'fit_measures.csv'),row.names=FALSE)
            neg <- sum(pe$op == '~~' & pe$lhs == pe$rhs & pe$est < -1e-6, na.rm=TRUE)
            over <- sum(pe$op == '=~' & abs(pe$std.all) > 1+1e-6, na.rm=TRUE)
            grad <- max(abs(lavaan::lavInspect(fit,'optim.gradient')))
            post <- isTRUE(lavaan::lavInspect(fit,'post.check'))
            finite <- all(is.finite(pe$se[pe$op %in% c('~','=~',':=')]))
            admissible <- post && neg == 0 && over == 0 && finite && is.finite(grad) && grad < 1e-3
            assessable <- spec$measurement != 'pca' && is.finite(measures['df.scaled']) && measures['df.scaled'] > 0
            screen <- if (assessable) isTRUE(all(is.finite(measures)) && measures['cfi.scaled'] >= .90 &&
                measures['tli.scaled'] >= .90 && measures['rmsea.scaled'] <= .08 && measures['srmr'] <= .08) else NA
            data.frame(model=spec$id, status='completed', N=nrow(z), n_met=length(m), n_inf=length(f),
                converged=converged, negative_variances=neg, loadings_over_one=over, max_gradient=grad,
                admissible=admissible, fit_assessable=assessable, descriptive_fit_screen=screen,
                CFI=unname(measures['cfi.scaled']), TLI=unname(measures['tli.scaled']),
                RMSEA=unname(measures['rmsea.scaled']), SRMR=unname(measures['srmr']),
                interpretation=if (!admissible) 'Improper/unreliable; do not interpret mediation' else if (!assessable)
                    'Observed-score sensitivity; fit cannot rank this model' else if (!screen)
                    'Inadequate global fit' else 'Eligible for substantive review; not independent validation')
        }, warning=function(w) {warnings <<- c(warnings, conditionMessage(w)); invokeRestart('muffleWarning')}),
        error=function(e) data.frame(model=spec$id, status='failed_or_skipped', interpretation=conditionMessage(e)))
        writeLines(unique(warnings), file.path(folder,'warnings.txt'))
        rows[[spec$id]] <- result
        cols <- unique(unlist(lapply(rows,names)))
        complete_rows <- lapply(rows,function(r) {for (nm in setdiff(cols,names(r))) r[[nm]] <- NA; r[,cols]})
        write.csv(do.call(rbind,complete_rows),file.path(OUT,'model_comparison.csv'),row.names=FALSE)
        if (length(all_effects)) write.csv(do.call(rbind,all_effects),file.path(OUT,'all_indirect_effects.csv'),row.names=FALSE)
        gc()
    }
    capture.output(sessionInfo(), file=file.path(OUT,'sessionInfo.txt'))
    message('Complete. Review model_comparison.csv and README.txt; no P-value-based winner selected.')
}
main()
}))

# ===========================================================================
# mechanism / summary
# ===========================================================================
dose_register_step("mechanism", "summary", expression({
# Mechanism overview: a-j, vector PDFs, no statistical models are refitted.
# Optional environment variables: DOSE_PROJECT, S4_PLOT_OUT, S4_SEM_MODEL.
# Equal axes means equal physical panel lengths, NOT equal numerical units.
ROOT <- dose_output_dir()
OUT <- Sys.getenv('S4_PLOT_OUT', file.path(ROOT, 'ukb/05mechanism/summary_figures'))
SEM_MODEL <- Sys.getenv('S4_SEM_MODEL', 'top10_met_all_inf')
if (dir.exists(file.path(ROOT, '.r-lib'))) .libPaths(c(file.path(ROOT, '.r-lib'), .libPaths()))
pkgs <- c('ggplot2','ggrepel','patchwork','jsonlite','dplyr','stringr')
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]
if(length(missing)) stop('Install packages: ', paste(missing,collapse=', '))
suppressPackageStartupMessages({library(ggplot2); library(patchwork)})
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
TEXT_SCALE <- 1.95 # Uniform font-size multiplier, independent of data-point sizes.
POINT_SCALE <- 1.6875 # Common diameter multiplier for every data-point panel.
POS <- '#FF684F'; NEG <- '#86CFF1'; NS <- '#CFCFCF'
COL <- c(Negative=NEG, Positive=POS, 'Not significant'=NS)
NODE <- c(DOSE='#C6DBE7', IS='#F0CDCE', Metabolites='#CBC9E2', Inflammation='#F6EBD5')
mechanism <- file.path(ROOT,'ukb/05mechanism')
inputs <- character()
read_data <- function(path) {
    inputs <<- unique(c(inputs,path))
    read.csv(path,check.names=FALSE)
}
star <- function(p) ifelse(is.na(p),'',ifelse(p<.001,'***',ifelse(p<.01,'**',ifelse(p<.05,'*','ns'))))
base_theme <- theme_classic(base_size=10,base_family='Arial') +
    theme(axis.text=element_text(color='#202020'), axis.title=element_text(size=10),
          legend.title=element_text(size=9),legend.text=element_text(size=8),
          plot.margin=margin(10,12,10,10), plot.title=element_text(size=11,face='bold'))

# Reuse the original data preparation and plot definitions, stopping before exports.
# Existing analysis scripts and their output PDFs are not changed.
load_definition <- function(module, step, stop_assignment) {
    path <- file.path(dose_script_dir(),paste0(module,".R"))
    inputs <<- unique(c(inputs,path))
    env <- new.env(parent=globalenv())
    for(expr in dose_step_code(module,step)) {
        head <- if(is.call(expr)) as.character(expr[[1]]) else ''
        if(head %in% c('<-','=') && identical(as.character(expr[[2]]),stop_assignment)) break
        if(head %in% c('dir.create','write.csv','writeLines','ggsave','bounded_ggsave')) next
        eval(expr,env)
    }
    env
}
old_project <- Sys.getenv('DOSE_PROJECT',unset=NA_character_)
Sys.setenv(DOSE_PROJECT=ROOT)
sig <- load_definition('signature','metabolite_plot','save_plot')
if(is.na(old_project)) Sys.unsetenv('DOSE_PROJECT') else Sys.setenv(DOSE_PROJECT=old_project)
# a: editable vector illustration, using actual screened/selected counts.
# The ball-and-stick motifs are schematic metabolites, not molecular structures.
screened <- read_data(file.path(mechanism,'metabolite/M2/01a_DOSE_metabolite_XM_all.csv'))
p_a <- ggplot()+coord_fixed(xlim=c(0,8),ylim=c(0,4.5),expand=FALSE,clip='off')+
    theme_void(base_family='Arial')+theme(plot.margin=margin(10,8,8,8))
poly_a <- function(x,y,fill,col=NA,lw=.4) annotate('polygon',x=x,y=y,fill=fill,colour=col,linewidth=lw)
ellipse_a <- function(x,y,rx,ry,fill,col=NA) {
    t<-seq(0,2*pi,length.out=80);poly_a(x+rx*cos(t),y+ry*sin(t),fill,col)
}
round_a <- function(x0,x1,y0,y1,fill,col=NA,r=.16) {
    theta<-seq(0,pi/2,length.out=15)
    xx<-c(x1-r+r*cos(theta),x0+r+r*cos(theta+pi/2),x0+r+r*cos(theta+pi),x1-r+r*cos(theta+3*pi/2))
    yy<-c(y1-r+r*sin(theta),y1-r+r*sin(theta+pi/2),y0+r+r*sin(theta+pi),y0+r+r*sin(theta+3*pi/2))
    poly_a(xx,yy,fill,col)
}
text_a <- function(x,y,label,size=3.1,face='plain',color='#263840')
    annotate('text',x=x,y=y,label=label,size=size,fontface=face,colour=color,lineheight=1.15)
p_a<-p_a+text_a(1.95,4.27,paste0('UK Biobank | ',format(sig$summary$n_participants,big.mark=',',trim=TRUE),' participants'),3.1,'bold')+
    round_a(1.65,3.9,.42,3.95,'#F8F2DE','#C9B895')+
    text_a(2.775,3.68,paste0('Metabolites: ',nrow(screened)),3.4,'bold')
# Centrifuged blood tube, including plasma, buffy coat, erythrocytes and cap.
p_a<-p_a+poly_a(c(.22,1.22,1.22,1.12,.83,.61,.32,.22),c(3.54,3.54,1.02,.76,.39,.39,.76,1.02),'#E7F0F3','#78858A',.65)+
    poly_a(c(.29,1.15,1.15,1.04,.8,.64,.4,.29),c(1.32,1.32,1.02,.79,.48,.48,.79,1.02),'#84212A')+
    annotate('rect',xmin=.29,xmax=1.15,ymin=1.32,ymax=1.95,fill='#EEDCC8')+
    annotate('rect',xmin=.29,xmax=1.15,ymin=1.95,ymax=2.98,fill='#E2BD42')+
    ellipse_a(.72,2.98,.43,.09,'#F2D169','#C3A546')+
    annotate('segment',x=.38,xend=.38,y=1.12,yend=3.38,colour='white',alpha=.65,linewidth=2.2)+
    round_a(.12,1.32,3.45,4.0,'#D6D8D5','#7B8181',.1)+
    ellipse_a(.72,4.0,.6,.13,'#E6E7E4','#7B8181')+
    ellipse_a(.72,4.02,.43,.07,'#CED1CE','#A1A5A1')
for(x in seq(.23,1.21,length.out=9)) p_a<-p_a+annotate('segment',x=x,xend=x,y=3.51,yend=3.91,colour='#A9AEAA',linewidth=.5)
for(y in seq(1.45,3.3,length.out=7)) p_a<-p_a+annotate('segment',x=.54,xend=.96,y=y,yend=y,colour='#86918D',linewidth=.4)
p_a<-p_a+annotate('path',x=c(.91,1.4,1.65),y=c(2.5,3.13,3.13),linewidth=.4,colour='#53636B')+
    text_a(.72,.16,'Plasma',2.8)
# Repeated small-molecule symbols retain a scientific illustration style in PDF.
mc<-c('#8AB875','#6CA8BC','#CB9BB8','#E5AD6B','#9C93C7','#74B5A3')
for(i in 1:12) {
    cx<-c(2.1,2.78,3.46)[(i-1)%%3+1];cy<-3.18-floor((i-1)/3)*.72
    t<-seq(0,2*pi,length.out=6)[1:5]+i*.4
    xx<-cx+c(0,.19*cos(t)); yy<-cy+c(0,.19*sin(t))
    for(k in 2:6) p_a<-p_a+annotate('segment',x=xx[1],y=yy[1],xend=xx[k],yend=yy[k],colour=mc[(i-1)%%6+1],linewidth=.9)
    for(k in 1:6) p_a<-p_a+ellipse_a(xx[k],yy[k],if(k==1) .10 else .065,if(k==1) .10 else .065,mc[(i-1)%%6+1])
}
# Right-hand statistical filtering, matching the reference's broad arrows/boxes.
p_a<-p_a+text_a(5.95,3.83,'Multivariable linear\nregression (M2)',3.8,'bold')+
    annotate('segment',x=4.25,xend=6.0,y=3.16,yend=3.16,colour='#D9CCE4',linewidth=3.2)+
    annotate('segment',x=6.0,xend=6.0,y=3.16,yend=2.82,colour='#D9CCE4',linewidth=3.2,
        arrow=grid::arrow(type='closed',length=grid::unit(3,'mm')))+
    round_a(4.3,7.85,1.99,2.74,'#D9CCE4')+
    text_a(6.075,2.48,paste0('DOSE-associated: ',sig$summary$n_fdr_features),3.5,'bold')+
    text_a(6.075,2.18,paste0('BH-FDR < 0.05 | QC retained: ',sig$summary$n_qc_features),2.7)+
    annotate('segment',x=4.8,xend=4.8,y=1.87,yend=1.16,colour='#AED2B9',linewidth=3.1,
        arrow=grid::arrow(type='closed',length=grid::unit(3,'mm')))+
    text_a(6.3,1.67,'LASSO',4,'bold')+
    text_a(6.3,1.35,paste0(sig$summary$folds,'-fold cross-validation'),2.8)+
    round_a(4.3,7.85,.37,1.06,'#AED2B9')+
    text_a(6.075,.80,paste0('Metabolomic signature: ',sig$summary$n_selected),3.5,'bold')+
    text_a(6.075,.53,'selected metabolites',2.9)+
    text_a(5.95,.13,'Median imputation | log1p | standardization',2.5)
# b: table-style forest plot; the numerical axis spans only the forest column.
forest_b <- sig$forest
forest_b$estimate[forest_b$group=='Q1'] <- 'Ref.'
forest_range <- range(c(.6,1,forest_b$LCI,forest_b$HCI),finite=TRUE)
forest_ticks <- pretty(forest_range,n=4)
forest_ticks <- forest_ticks[forest_ticks > 0]
forest_range <- range(c(forest_range,forest_ticks))
forest_x <- function(hr) .29+(hr-forest_range[1])/diff(forest_range)*.34
p_b <- ggplot(forest_b)+
    geom_segment(aes(x=forest_x(LCI),xend=forest_x(HCI),y=y,yend=y),colour='#78C9AD',linewidth=1)+
    geom_point(aes(x=forest_x(HR),y=y),colour='#379AAF',size=3.5*POINT_SCALE)+
    geom_text(aes(x=forest_x(HR),y=y+.40,label=stars),size=4,fontface='bold')+
    geom_text(aes(x=.22,y=y,label=group),hjust=1,size=3.4)+
    geom_text(aes(x=.74,y=y,label=estimate),hjust=0,size=3.3)+
    annotate('segment',x=forest_x(1),xend=forest_x(1),y=.28,yend=5.92,colour='#CE8778',linetype='dashed',linewidth=.65)+
    annotate('segment',x=.01,xend=1.04,y=5.34,yend=5.34,linewidth=.6)+
    annotate('text',x=.015,y=5.67,label='DOSE diet metabolomic signature',hjust=0,size=3.5,fontface='bold')+
    annotate('text',x=.74,y=5.67,label='HR (95% CI)',hjust=0,size=3.5,fontface='bold')+
    annotate('segment',x=.27,xend=.685,y=.25,yend=.25,linewidth=.6)+
    annotate('segment',x=forest_x(forest_ticks),xend=forest_x(forest_ticks),y=.25,yend=.17,linewidth=.5)+
    annotate('text',x=forest_x(forest_ticks),y=-.03,label=format(forest_ticks,trim=TRUE),size=3.2)+
    annotate('text',x=.475,y=-.48,label='Hazard ratio',size=3.5)+
    coord_cartesian(xlim=c(0,1.07),ylim=c(-.65,6),expand=FALSE,clip='off')+
    theme_void(base_family='Arial')+theme(plot.margin=margin(14,12,8,12))
# c: original Cox-product results on the coefficient scale (log HR), not exp(beta).
# a is the linear DOSE -> signature coefficient; b adjusts for DOSE.
med_c <- sig$med
stopifnot(isTRUE(all.equal(med_c$ACME,med_c$a_estimate*med_c$b_adjusted_estimate,tolerance=1e-10)))
beta_label <- function(effect,estimate,p) sprintf('\u03b2 %s = %.4f%s',effect,estimate,star(p))
c_box <- function(x,y,label,fill,border) annotate('label',x=x,y=y,label=label,
    fill=fill,colour=border,size=3.3,fontface='bold',linewidth=.45,
    label.padding=grid::unit(.65,'lines'),label.r=grid::unit(.35,'lines'))
c_arrow <- function(x,y,xe,ye) annotate('segment',x=x,y=y,xend=xe,yend=ye,
    colour='#293C48',linewidth=.7,arrow=grid::arrow(length=grid::unit(3,'mm'),type='closed'))
p_c <- ggplot()+
    annotate('text',x=.5,y=.96,label=beta_label('indirect',med_c$ACME,med_c$P_ACME),size=3.3)+
    c_box(.5,.79,'DOSE-related\nmetabolomic signature','#E3ECD9','#617F45')+
    c_box(.12,.18,'DOSE',NODE['DOSE'],'#476B80')+
    c_box(.88,.18,'IS',NODE['IS'],'#955E63')+
    c_arrow(.19,.29,.34,.66)+c_arrow(.66,.66,.81,.29)+c_arrow(.27,.18,.73,.18)+
    annotate('text',x=.20,y=.48,label=beta_label('a',med_c$a_estimate,med_c$a_P),angle=58,size=3.1)+
    annotate('text',x=.80,y=.48,label=beta_label('b',med_c$b_adjusted_estimate,med_c$b_adjusted_P),angle=-58,size=3.1)+
    annotate('text',x=.5,y=.48,size=3.3,lineheight=1.3,
        label=paste0(beta_label('total',med_c$c_total_model_estimate,med_c$c_total_model_P),
            sprintf('\nMediation proportion: %.2f%%%s\n95%% CI %.2f%% - %.2f%%',
                100*med_c$Proportion_mediated,star(med_c$P_proportion),
                100*med_c$Proportion_LCI,100*med_c$Proportion_HCI)))+
    annotate('text',x=.5,y=.09,label=beta_label('direct',med_c$c_direct_model_estimate,med_c$c_direct_model_P),size=3.3)+
    annotate('text',x=.5,y=-.04,label='Cox effects on the log-hazard scale',size=2.5,colour='#586069')+
    coord_cartesian(xlim=c(0,1),ylim=c(-.09,1.04),expand=FALSE,clip='off')+
    theme_void(base_family='Arial')+theme(plot.margin=margin(15,20,15,20))

# d: derive the plotting table from this run and the supplied metabolite dictionary.
met <- read_data(file.path(mechanism,'metabolite/M2/01a_DOSE_metabolite_XM_all.csv'))
met_dict_path <- file.path(dose_ukb_dir(), 'common/met.lst')
met_dict <- read.delim(met_dict_path, check.names=FALSE, stringsAsFactors=FALSE)
inputs <- unique(c(inputs, met_dict_path))
if (all(c('met_name','group','full_name') %in% names(met_dict))) {
    ix <- match(sub('^met_', '', met$metabolite), sub('^met_', '', met_dict$met_name))
    met$category <- gsub('_',' ',met_dict$group[ix],fixed=TRUE)
    met$label <- gsub('_',' ',met_dict$full_name[ix],fixed=TRUE)
} else {
    met$category <- NA_character_; met$label <- NA_character_
}
met$category[is.na(met$category) | !nzchar(met$category)] <- 'Unclassified'
met$label[is.na(met$label) | !nzchar(met$label)] <- sub('^met_', '', met$metabolite[is.na(met$label) | !nzchar(met$label)])
met$abs_beta <- abs(met$beta_XM)
category_out <- file.path(mechanism,'metabolite/M2/category_plot')
dir.create(category_out, recursive=TRUE, showWarnings=FALSE)
write.csv(met, file.path(category_out,'DOSE_metabolite_category_plot_data.csv'),row.names=FALSE)
met <- met[order(met$category,-met$beta_XM,met$metabolite),]
met$x <- seq_len(nrow(met)); met$direction <- ifelse(met$FDR_XM<.05,ifelse(met$beta_XM>0,'Positive','Negative'),'Not significant')
cats <- unique(met$category)
# Add visual spacing between categories; x is an ordering index, not a quantity.
met$x <- met$x + 7*(match(met$category,cats)-1)
band <- do.call(rbind,lapply(cats,function(k) {
    z <- met$x[met$category==k]; data.frame(category=k,xmin=min(z)-4,xmax=max(z)+4,mid=mean(range(z)))
}))
band$color <- c('#9DCADD','#247AA7','#B1DD8A','#31A343','#F9ADC5','#E92325','#FDBF6F','#FF8000',
                '#CBB2D6','#6A3D9A','#FFFF99','#B15928','#76C7B7','#F2EDA0','#BBB8CB','#FA7774')[seq_len(nrow(band))]
yr <- range(met$beta_XM); span <- diff(yr); band$ymin <- yr[1]-.075*span; band$ymax <- yr[1]-.055*span
illustrative <- c('met_Omega_3.pct','met_DHA.pct','met_DHA','met_Unsaturation','met_Omega_3','met_Albumin','met_GlycA','met_Omega_6_by_Omega_3','met_S_HDL_TG.pct','met_L_HDL_C','met_XL_HDL_C','met_LDL_size')
label_met <- met[met$metabolite %in% unique(c(head(met$metabolite[order(-abs(met$beta_XM))],6),illustrative)),]
p_d <- ggplot(met,aes(x,beta_XM)) +
    geom_rect(data=band,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),inherit.aes=FALSE,fill=band$color)+
    geom_hline(yintercept=0,linetype='dotted',color='#777777')+
    geom_point(aes(fill=direction,size=abs_beta),shape=21,stroke=.18,color='#333333',alpha=.95)+
    ggrepel::geom_text_repel(data=label_met,aes(label=label),size=2.4,seed=42,max.overlaps=Inf,
        min.segment.length=0,box.padding=.25,max.time=2)+
    scale_fill_manual(values=COL,name='Association')+
    scale_size_area(max_size=5.2*POINT_SCALE,breaks=c(.02,.04,.06,.08),name='Absolute beta')+
    scale_x_continuous(breaks=band$mid,labels=band$category,expand=expansion(mult=.01))+
    scale_y_continuous(expand=expansion(mult=c(.025,.14)))+
    labs(x=NULL,y='Beta per 1-SD higher DOSE')+base_theme+
    theme(aspect.ratio=.5,axis.text.x=element_text(angle=55,hjust=1,vjust=1,size=7),axis.ticks.x=element_blank(),
          legend.position='right',legend.direction='vertical',legend.justification='top',legend.key.size=grid::unit(.3,'cm'))+
    guides(fill=guide_legend(override.aes=list(size=3*POINT_SCALE)),size=guide_legend(ncol=1))

# Pool both domains so identical mediation percentages have identical areas.
mediation_data <- lapply(c('metabolite','inflammation'),function(domain)
    read_data(file.path(mechanism,domain,'M2',paste0('Fig_DOSE_',domain,'_mediation_beta_bubble_plot_data.csv'))))
names(mediation_data) <- c('metabolite','inflammation')
MEDIATION_MAX <- max(4.5,unlist(lapply(mediation_data,function(d) d$proportion_percent)))
bubble <- function(domain) {
    d <- mediation_data[[domain]]
    stopifnot(all(is.finite(d$beta_a)),all(is.finite(d$beta_b)),all(d$proportion_percent>0))
    d$direction <- ifelse(d$beta_a<0,'Negative','Positive')
    d$label <- d[[domain]]
    labels <- if(domain=='metabolite') head(d[order(-d$proportion_percent),],8) else d
    ggplot(d,aes(abs(beta_a),abs(beta_b)))+
        geom_point(aes(size=proportion_percent,fill=direction),shape=21,color='white',stroke=.3,alpha=.9)+
        ggrepel::geom_text_repel(data=labels,aes(label=label),size=2.7,seed=42,box.padding=.4,
            point.padding=.5,min.segment.length=0,max.overlaps=Inf,max.time=2)+
        scale_fill_manual(values=COL[c('Negative','Positive')],drop=FALSE,name='DOSE association')+
        scale_size_area(max_size=10*POINT_SCALE,breaks=c(.5,1,2,4),limits=c(0,MEDIATION_MAX),name='Proportion mediated (%)')+
        scale_x_continuous(limits=c(0,NA),expand=expansion(mult=c(.02,.15)))+
        scale_y_continuous(limits=c(0,NA),expand=expansion(mult=c(.02,.15)))+
        labs(x=expression('Absolute value of '*beta[a]),y=expression('Absolute value of '*beta[b]))+
        base_theme+theme(aspect.ratio=1,legend.position='right',legend.direction='vertical',legend.justification='top',legend.box='vertical',legend.margin=margin(0,0,0,0))+
        guides(fill=guide_legend(order=1,override.aes=list(size=3*POINT_SCALE)),size=guide_legend(order=2,ncol=1,override.aes=list(fill='#444444')))
}
p_e <- bubble('metabolite'); p_g <- bubble('inflammation')
volcano <- function(domain) {
    d <- read_data(file.path(mechanism,domain,'M2',paste0('01a_DOSE_',domain,'_XM_all.csv')))
    stopifnot(all(is.finite(d$beta_XM)),all(is.finite(d$P_XM)),all(is.finite(d$FDR_XM)))
    d$direction <- ifelse(d$FDR_XM<.05,ifelse(d$beta_XM>0,'Positive','Negative'),'Not significant')
    d$label <- sub('^(bc_|bb_)','',d[[domain]])
    labels <- if(domain=='inflammation') d else do.call(rbind,lapply(c('Negative','Positive'),function(s) {
        z<-d[d$direction==s,]; head(z[order(z$P_XM,-abs(z$beta_XM)),],8)
    }))
    p <- ggplot(d,aes(beta_XM,-log10(pmax(P_XM,.Machine$double.xmin))))+
        geom_vline(xintercept=0,color='#B58AD0',linetype='dashed',linewidth=.4)+
        geom_point(aes(color=direction),size=POINT_SCALE*if(domain=='protein') 1.6 else 2.5,alpha=.8)+
        ggrepel::geom_text_repel(data=labels,aes(label=label,color=direction),size=2.5,seed=42,
            box.padding=.3,min.segment.length=0,max.overlaps=Inf,max.time=2,show.legend=FALSE)+
        scale_color_manual(values=COL,breaks=c('Negative','Not significant','Positive'),name=NULL,drop=FALSE)+
        scale_y_continuous(limits=c(0,NA),expand=expansion(mult=c(.025,.14)))+
        scale_x_continuous(expand=expansion(mult=c(.12,.15)))+
        labs(x=paste0('Beta for DOSE-',domain,' association'),y=expression(-log[10](P)))+
        base_theme+theme(aspect.ratio=1,legend.position='right',legend.direction='vertical',legend.justification='top')
    if(any(d$FDR_XM<.05)) p<-p+geom_hline(yintercept=-log10(max(d$P_XM[d$FDR_XM<.05])),color='#B58AD0',linetype='dashed',linewidth=.4)
    p
}
p_f <- volcano('inflammation'); p_i <- volcano('protein')

# h: reuse the original SEM drawing functions and geometry without refitting.
# Evaluate only the preparation/definition portion of plot_one; skip all writes.
sem_dir <- file.path(mechanism,'sem/M2_Cox',SEM_MODEL)
sem_source <- file.path(dose_script_dir(), 'mechanism.R')
inputs <- unique(c(inputs,sem_source,list.files(sem_dir,pattern='csv$',full.names=TRUE)))
sem <- new.env(parent=globalenv())
sem$DIR_IN <- sem_dir
sem_expr <- dose_step_code('mechanism','sem_plot')
plot_one_expr <- Filter(function(e) is.call(e) && identical(e[[1]],as.name('<-')) &&
    identical(e[[2]],as.name('plot_one')),as.list(sem_expr))[[1]]
for (expr in as.list(plot_one_expr[[3]][[3]])[-1]) {
    head <- if(is.call(expr)) as.character(expr[[1]]) else ''
    if(head %in% c('<-','=') && identical(as.character(expr[[2]]),'stem')) break
    if(head %in% c('dir.create','write.csv','writeLines')) next
    eval(expr,sem)
}
# Retain the curved gradient arrows, original node positions and right-hand tables.
# Only replace the original internal D./E. labels with this figure's h label.
body(sem$draw) <- as.call(c(list(as.name('{')),Filter(function(e) {
    !(is.call(e) && (identical(e[[1]],as.name('grid.newpage')) ||
        (identical(e[[1]],as.name('txt')) && is.character(e[[2]]) && e[[2]] %in% c('D.','E.'))))
},as.list(body(sem$draw))[-1])))
sem$cols[c('dose','is')] <- NODE[c('DOSE','IS')]
paths <- sem$paths; selected <- sem$selected
p_h <- NULL
# Match the original 12:4.8 geometry and scale text/strokes with the viewport.
draw_h <- function(width,height) {
    scale <- min((width-.15)/12,(height-.35)/4.8)
    grid::grid.text('h',x=grid::unit(10,'pt'),y=grid::unit(1,'npc')-grid::unit(10,'pt'),
        just=c('left','top'),gp=grid::gpar(fontsize=16*TEXT_SCALE,fontface='bold',fontfamily='Arial'))
    grid::pushViewport(grid::viewport(width=grid::unit(12*scale,'in'),
        height=grid::unit(4.8*scale,'in'),gp=grid::gpar(cex=scale,lex=scale)))
    sem$draw()
    grid::popViewport()
}
# j: exactly the default-background enrichment specified by the user.
enrichment <- load_definition('mechanism','pathway_plot','name')
p_j <- enrichment$p +
    scale_size_area(max_size=5.8*POINT_SCALE,limits=c(0,max(enrichment$d$Count)),
        breaks=c(5,10,20,30,39),name='Overlap') +
    theme(legend.position='right',legend.direction='vertical')

# Scale theme text and annotation/repel text together; never scale point geoms.
enlarge_text <- function(p) {
    for(k in seq_along(p$layers)) {
        lyr <- p$layers[[k]]
        if(any(grepl('Text|Label',class(lyr$geom)))) {
            old <- lyr$aes_params$size
            if(is.null(old)) old <- 3.88
            p$layers[[k]]$aes_params$size <- old*TEXT_SCALE
        }
    }
    th <- theme_get()+p$theme
    for(nm in names(th)) {
        if(inherits(th[[nm]],'element_text') && is.numeric(th[[nm]]$size) && !inherits(th[[nm]]$size,'rel'))
            th[[nm]]$size <- th[[nm]]$size*TEXT_SCALE
    }
    p+th
}
# Preserve SEM geometry while enlarging its labels and loading-table text.
eval(quote(txt <- function(label,x,y,size=11,col=ink,face='plain',just='centre',rot=0)
    grid.text(label,x,y,just=just,rot=rot,
        gp=gpar(fontsize=size*TEXT_SCALE,col=col,fontface=face,fontfamily='sans'))),sem)
plots <- list(a=p_a,b=p_b,c=p_c,d=p_d,e=p_e,f=p_f,g=p_g,h=p_h,i=p_i,j=p_j)
plots[names(plots)!='h'] <- lapply(plots[names(plots)!='h'],enlarge_text)
names_long <- c('signature_selection','signature_cox','signature_mediation','metabolite_association',
    'metabolite_mediation','inflammation_association','inflammation_mediation','sem','protein_association','pathway_enrichment')
# Each independent cell is a gtable; this preserves square panels in the composite.
# Use the same physical cell dimensions for individual exports and overview.
cells <- data.frame(letter=letters[1:10],x=c(0,10.5,21,0,21,0,10,20,0,10),
    y=c(0,0,0,6,6,16.5,16.5,16.5,24.5,24.5),
    w=c(10.5,10.5,11,21,11,10,10,12,10,22),h=c(6,6,6,10.5,10.5,8,8,8,8,8))
PAGE_WIDTH <- 32; PAGE_HEIGHT <- 32.5
# Fix the panel tracks in inches, rather than merely setting aspect.ratio.
# The same gtable is drawn at its native size in both individual and composite PDFs.
PANEL_SIDE <- 5.3
fixed_letters <- c('d','e','f','g','i')
fixed_grobs <- list()
fix_panel <- function(p,letter) {
    g <- ggplotGrob(p)
    z <- g$layout[g$layout$name=='panel',]
    stopifnot(nrow(z)==1)
    g$widths[z$l:z$r] <- grid::unit(if(letter=='d') 3.125*PANEL_SIDE else PANEL_SIDE,'in')
    g$heights[z$t:z$b] <- grid::unit(if(letter=='d') 6.6 else PANEL_SIDE,'in')
    g$respect <- FALSE
    g
}
draw_fixed <- function(g,width,height) {
    gw <- grid::convertWidth(sum(g$widths),'in',valueOnly=TRUE)
    gh <- grid::convertHeight(sum(g$heights),'in',valueOnly=TRUE)
    stopifnot(gw <= width,gh <= height)
    grid::pushViewport(grid::viewport(width=grid::unit(gw,'in'),height=grid::unit(gh,'in')))
    grid::grid.draw(g)
    grid::popViewport()
}
tag_plot <- function(p,letter) p+labs(tag=letter)+theme(plot.tag=element_text(face='bold',size=16*TEXT_SCALE),plot.tag.position='topleft')
for(i in seq_along(plots)) {
    if(names(plots)[i]=='h') {
        cairo_pdf(file.path(OUT,'h_sem.pdf'),width=12,height=5.1)
        grid::grid.newpage(); draw_h(12,5.1); bounded_dev_off()
        next
    }
    plots[[i]]<-tag_plot(plots[[i]],names(plots)[i])
    if(names(plots)[i] %in% fixed_letters) {
        cairo_pdf(file.path(OUT,paste0(names(plots)[i],'_',names_long[i],'.pdf')),width=cells$w[i],height=cells$h[i])
        grid::grid.newpage()
        fixed_grobs[[names(plots)[i]]] <- fix_panel(plots[[i]],names(plots)[i])
        draw_fixed(fixed_grobs[[names(plots)[i]]],cells$w[i],cells$h[i]); bounded_dev_off()
        next
    }
    bounded_ggsave(file.path(OUT,paste0(names(plots)[i],'_',names_long[i],'.pdf')),plots[[i]],
        width=cells$w[i],height=cells$h[i],device=cairo_pdf,bg='white',limitsize=FALSE)
}
draw_overview <- function() {
    grid::grid.newpage()
    for(i in seq_along(plots)) {
        z<-cells[i,]
        vp <- grid::viewport(x=z$x/PAGE_WIDTH,y=1-z$y/PAGE_HEIGHT,
            width=z$w/PAGE_WIDTH,height=z$h/PAGE_HEIGHT,just=c('left','top'))
        if(names(plots)[i]=='h') {
            grid::pushViewport(vp); draw_h(z$w,z$h); grid::popViewport()
        } else if(names(plots)[i] %in% fixed_letters) {
            grid::pushViewport(vp); draw_fixed(fixed_grobs[[names(plots)[i]]],z$w,z$h); grid::popViewport()
        } else print(plots[[i]],newpage=FALSE,vp=vp)
    }
}
cairo_pdf(file.path(OUT,'Fig_s4_mechanism_summary.pdf'),width=PAGE_WIDTH,height=PAGE_HEIGHT)
draw_overview(); bounded_dev_off()
bounded_png(file.path(OUT,'Fig_s4_mechanism_summary.png'),width=round(PAGE_WIDTH*150),height=round(PAGE_HEIGHT*150),res=150,type='cairo')
draw_overview(); bounded_dev_off()
write.csv(cells,file.path(OUT,'panel_layout_inches.csv'),row.names=FALSE)
write.csv(paths[nzchar(paths$label)&!is.na(paths$label),],file.path(OUT,'h_SEM_paths.csv'),row.names=FALSE)
write.csv(selected,file.path(OUT,'h_SEM_displayed_loadings.csv'),row.names=FALSE)
result_inputs <- c(list.files(file.path(ROOT,'ukb/06signature/metabolite/association_mediation'),full.names=TRUE,pattern='csv$'),
    file.path(ROOT,'ukb/06signature/metabolite',c('metabolite_LASSO_summary.json','metabolite_LASSO_selected.csv')),
    file.path(enrichment$root,'01_all_five_databases.csv'))
inputs<-unique(c(inputs,result_inputs))
write.csv(data.frame(path=inputs,md5=unname(tools::md5sum(inputs))),file.path(OUT,'input_manifest.csv'),row.names=FALSE)
writeLines(c('Mechanism overview; all estimates are existing results; no models refitted.',
    'a: signature construction; b: M2 Cox; c: Cox-product signature mediation; a is linear beta, b/direct/total are Cox log-HR coefficients; indirect=a*b; DOSE-adjusted b. No exponentiation.',
    'd/f/i: all screened associations; significant means source BH-FDR < 0.05; nonsignificant points gray.',
    'e/g: saved upstream Cox mediation plot data; |beta_a| versus |beta_b|; bubble AREA represents calculated mediation percentage.',
    'd/e/f/g/i: positive DOSE association red #FF684F, negative blue #86CFF1.',
    'c/h: DOSE light blue #C6DBE7, IS light red #F0CDCE.',
    'd: y-axis length 6.6 inches; e/f/g/i: y-axis length exactly 5.3 inches. e/f/g/i: x-axis length exactly 5.3 inches; d: 16.5625 inches. Fixed gtable panel tracks in all exports.',
    paste('e/g shared bubble scale: area proportional to mediation percent; limits 0 to',MEDIATION_MAX,'percent; max diameter 16.875 mm; shared breaks 0.5, 1, 2, 4.'),
    paste('h SEM model:',SEM_MODEL),
    'h: std.all structural coefficients and three strongest absolute standardized loadings per domain; all ten metabolites remain in the fitted model.',
    'Stars: * P<0.05, ** P<0.01, *** P<0.001; ns: P>=0.05. SEM P values are original parameter tests.',
    'f/i: y is raw -log10(P); dashed horizontal line is largest observed raw P passing BH-FDR<0.05.',
    'j: user-specified Enrichr default background; three KEGG, six Reactome and six GO BP terms; not the measured-background analysis.',
    'Overview is 32 x 32.5 inches, with a-c / d-e / f-h / i-j rows. h preserves the original SEM 12:4.8 drawing geometry, gradient curves and right-hand loading tables.',
    'Individual PDFs and overview contain vector text and shapes.'),file.path(OUT,'README.txt'))
message('Saved 10 individual vector PDFs and overview PDF/PNG to: ',OUT)
}))

if (sys.nframe() == 0L) dose_run_module("mechanism")
