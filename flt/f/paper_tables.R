source(file.path(Sys.getenv("FLT_CODE"), "runtime.R"))
library(data.table)
root <- flt_root
paper <- Sys.getenv("FLT_PAPER")
reports <- Sys.getenv("FLT_PAPER_REPORTS")
tables <- list()
coverage <- list()

add <- function(name, file, previous = "", eligible = TRUE) {
    exists <- file.exists(file) && eligible
    coverage[[name]] <<- data.table(table = name, status = if (exists) "available" else "unavailable",
                               source = file, previous_manuscript = previous)
    if (!exists) return(invisible(NULL))
    d <- fread(file, data.table = FALSE)
    if (any(tolower(names(d)) %in% c("eid", "seqn", "iid"))) stop("Participant-level table excluded: ", file)
    tables[[name]] <<- d
}
paths <- c(Participant_characteristics = "01cox/demographic_characteristics.csv",
           Cox_all_contrasts = "01cox/cox_results_table.csv",
           Cox_sensitivity = "01cox/cox-sensitivity.csv",
           General_liking_items = "01cox/general-liking-items.csv",
           Intake = "01intake/spearman_selected_actual_liking.csv",
           Intake_items = "01intake/spearman_selected_intake_items.csv",
           Intake_profiles = "01intake/radar_project_intake.csv",
           Intake_samples = "01intake/sample_counts.csv",
           Behavioral = "01lm/lm_results.csv", MRI = "01lm/mri/mri_lm_results.csv",
           Protein_limma = "01protein/limma/all_limma_results.csv",
           Concordant_food_protein = "01protein/limma/shared_same_direction_proteins.csv",
           Protein_linear = "01protein/linear/FLT_protein_linear_all.csv",
           Colocalization = "03placo/coloc-results.csv",
           PheWAS = "03phewas/lead_phewas_results.tsv", PheWAS_sources = "03phewas/source_coverage.tsv",
           PheWAS_traits = "03phewas/trait_counts.tsv", PheWAS_categories = "03phewas/trait_class_locus_counts.tsv",
           SMR_matched = "04smr/matched.tsv", SMR_probes = "04smr/probe-tests.tsv")
for (name in names(paths)) add(name, file.path(root, paths[[name]]))
add("PLACO_summary", flt_placo_file("_summary.tsv"))
add("PLACO_hits", flt_placo_file("_hits.tsv"))
add("Cohort_flow", file.path(Sys.getenv("FLT_REPORTS"), "cohort-flow.csv"))
mr <- c(MR_methods = "uvmr_main_table.csv", MR_raw = "uvmr_results.csv",
        MR_instruments = "uvmr_snp_counts.csv", MR_sensitivity = "uvmr_sensitivity_table.csv",
        MR_failures = "uvmr_failures.csv", MR_trait_metadata = "steiger_trait_info_used.csv",
        MR_bidirectional = "uvmr_bidirectional.csv", MR_leave_one_out = "uvmr_leave_one_out.csv")
for (name in names(mr)) add(name, file.path(flt_run("uvmr"), mr[[name]]))
mv <- c(MVMR_estimates = "mvmr_all_results.csv", MVMR_conditional_F = "mvmr_conditional_F.csv",
        MVMR_pleiotropy = "mvmr_conditional_Q.csv", MVMR_failures = "failed_models.csv",
        MVMR_pairs = "mvmr_pairs.csv", MVMR_input = "gwas_manifest.tsv")
for (name in names(mv)) add(name, file.path(flt_run("mvmr"), mv[[name]]))
add("Genetic_correlation", file.path(flt_run("ldsc"), "ldsc_hdl_results.csv"))
add("GWAS_panel", file.path(Sys.getenv("FLT_CONFIG"), "panels.tsv"))
annotations <- c(FUMA_loci_genes = "gene_mapping/G732-mdd_2025_Fig6a_locus_genes.tsv",
                 FUMA_genes = "gene_mapping/G732-mdd_2025_Fig6b_gene_evidence.tsv",
                 MAGMA_pathways = "MAGMA_pathways_all.csv", MAGMA_GTEx54 = "MAGMA_GTEx54_all.csv",
                 MAGMA_general30 = "MAGMA_GTEx_general_all.csv",
                 deTS_GTEx = "deTS/deTS_GTEx_all.csv", deTS_ENCODE = "deTS/deTS_ENCODE_all.csv",
                 Metascape = "metascape/Metascape_significant_cluster_representatives.csv")
fuma <- file.path(paper, "03fuma")
fuma_current <- file.exists(file.path(paper, "Fig5.png"))
for (name in names(annotations))
    add(name, file.path(fuma, annotations[[name]]),
        eligible = fuma_current && (name != "Metascape" || file.exists(Sys.getenv("FLT_METASCAPE_RESULTS"))))
for (type in c("cox", "rcs")) {
    name <- if (type == "cox") "protein_cox_all.csv" else "protein_rcs_tests_all.csv"
    for (food in c("G717", "G824", "G732")) {
        add(paste("Protein", type, food, sep = "_"), file.path(root, "01protein", type, food, name))
    }
}
cox <- tables$Cox_all_contrasts
pcx <- rbindlist(tables[grepl("^Protein_cox_", names(tables))], fill = TRUE)
count_if <- function(x, predicate) if (is.null(x)) NA_integer_ else sum(predicate, na.rm = TRUE)
checks <- data.table(
    Metric = c("Dictionary traits", "Primary prospective outcomes", "Planned continuous Cox pairs",
               "Concordant limma food-protein pairs", "COL18A1 among depression Cox hits",
               "Candidate genes meeting at least four criteria"),
    Previous = c(182, 5, 910, 15, 1, 23),
    Recomputed = c(if (is.null(cox)) NA_integer_ else uniqueN(cox$name),
                   if (is.null(cox)) NA_integer_ else uniqueN(cox$Disease),
                   count_if(cox, cox$Analysis == "continuous"),
                   if (is.null(tables$Concordant_food_protein)) NA_integer_ else nrow(tables$Concordant_food_protein),
                   if (!nrow(pcx)) NA_integer_ else as.integer(any(toupper(pcx$protein) == "COL18A1" &
                       pcx$outcome == "depression" & pcx$significant %in% TRUE)),
                   count_if(tables$FUMA_genes, tables$FUMA_genes$Evidence_count >= 4)))
checks[, comparison := fifelse(is.na(Recomputed), "unavailable",
                               fifelse(Previous == Recomputed, "matches", "differs"))]
tables$Coverage <- rbindlist(coverage)
tables$Draft_comparison <- checks
fwrite(tables$Coverage, file.path(reports, "Table-coverage.csv"))
fwrite(checks, file.path(reports, "Draft-comparison.csv"))
flt_workbooks(tables, file.path(paper, "FLT.source-tables.xlsx"))
