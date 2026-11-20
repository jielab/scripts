source(file.path(Sys.getenv("FLT_CODE"), "runtime.R"))
root <- Sys.getenv("FLT_OUT")
out <- file.path(root,"main_plot")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
read_table <- function(f) {
    p <- file.path(root, f)
    if (!file.exists(p)) stop("Missing computed table: ", p)
    d <- data.table::fread(p, data.table=FALSE)
    if (any(tolower(names(d)) %in% c("eid","seqn","iid"))) stop("Participant-level file excluded: ",p)
    d
}
tables <- list(
    Participant_characteristics=read_table("01cox/demographic_characteristics.csv"),
    Cox_continuous=flt_cox(), Cox_all_contrasts=read_table("01cox/cox_results_table.csv"),
    Forward_MR=flt_uv("forward"), Reverse_MR=flt_uv("reverse"), MVMR=flt_mv(),
    MR_sensitivity=as.data.frame(data.table::fread(file.path(flt_run("uvmr"),"uvmr_sensitivity_table.csv"))),
    Intake=read_table("01intake/spearman_selected_actual_liking.csv"),
    Intake_profiles=read_table("01intake/radar_project_intake.csv"),
    Behavioral=read_table("01lm/lm_results.csv"), MRI=read_table("01lm/mri/mri_lm_results.csv"),
    Concordant_food_protein=read_table("01protein/limma/shared_same_direction_proteins.csv"),
    Protein_linear=read_table("01protein/linear/FLT_protein_linear_all.csv"),
    Protein_limma=read_table("01protein/limma/all_limma_results.csv"),
    Genetic_correlation=as.data.frame(data.table::fread(file.path(flt_run("ldsc_hdl"),"ldsc_hdl_results.csv"))),
    PLACO_summary=read_table("03placo/G732__mdd_2025_PLACO_summary.tsv"),
    SMR_matched=read_table("04smr/matched.tsv"), FUMA_genes=read_table("03fuma/gene_mapping/G732__mdd_2025_Fig6b_gene_evidence.tsv"),
    MAGMA_pathways=read_table("03fuma/MAGMA_pathways_all.csv"), MAGMA_tissues=read_table("03fuma/MAGMA_GTEx54_all.csv"),
    Metascape=read_table("03fuma/metascape/Metascape_significant_cluster_representatives.csv"),
    PheWAS=read_table("03phewas/lead_phewas_results.tsv")
)
for (type in c("cox","rcs")) {
    name <- if (type=="cox") "protein_cox_all.csv" else "protein_rcs_tests_all.csv"
    paths <- list.files(file.path(root,"01protein",type),pattern=paste0("^",name,"$"),recursive=TRUE)
    if (!length(paths)) stop("Missing computed protein ",type," results")
    tables[[paste0("Protein_",type)]] <- data.table::rbindlist(lapply(paths,function(p) {
        d <- read_table(file.path("01protein",type,p)); d$pair_directory <- dirname(p); d
    }),fill=TRUE)
}
summary_checks <- data.frame(
    metric=c("Dictionary traits","Primary Cox outcomes","Cox trait-outcome pairs",
      "Genes supported by at least four evidence criteria"),
    manuscript=c(182,5,910,23),
    recomputed=c(length(unique(tables$Cox_continuous$name)),length(unique(tables$Cox_continuous$disease)),
      nrow(tables$Cox_continuous),sum(tables$FUMA_genes$Evidence_count >= 4, na.rm=TRUE)))
cox <- tables$Protein_cox
dep_hit <- unique(toupper(as.character(cox$protein[cox$outcome=="depression" & cox$significant %in% TRUE])))
summary_checks <- rbind(summary_checks, data.frame(metric=c("Concordant food-protein pairs","COL18A1 among depression Cox hits"),
    manuscript=c(15,1),recomputed=c(nrow(tables$Concordant_food_protein),as.integer("COL18A1" %in% dep_hit))))
summary_checks$comparison <- ifelse(summary_checks$manuscript==summary_checks$recomputed,"matches","differs")
tables$Manuscript_checks <- summary_checks
file <- tempfile(tmpdir=out, fileext=".xlsx")
openxlsx::write.xlsx(tables,file,asTable=FALSE,overwrite=TRUE)
if (!file.rename(file,file.path(out,"FLT.source-tables.xlsx"))) stop("Unable to publish workbook")
