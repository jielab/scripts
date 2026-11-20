# GWAS-export rules shared by the pipeline and small synthetic regressions.
ukb_gwas_traits <- c("cvd_cad", "t2dm", "stroke", "cvd_stroke_i", "stroke_o",
    "ad", "amr", "asthma", "copd", "cvd_afib", "cvd_htn", "epilepsy",
    "masld", "mi", "ra", "sle2")

ukb_gwas_batch <- function(gen_batch) {
    value <- trimws(as.character(gen_batch))
    value[is.na(gen_batch) | value %in% c("", "NA")] <- NA_character_
    if (any(!is.na(value) & !grepl("^-?[0-9]+$", value)))
        stop("gen_batch must contain integer UKB genotyping batch codes", call. = FALSE)
    # Negative UK BiLEVE batch codes are valid. Missing values remain NA.
    factor(ifelse(is.na(value), NA_character_, paste0("b", value)))
}

ukb_gwas_events <- function(dat, end_date) {
    required <- c("gen_batch", paste0("PC", 1:10),
        paste0("fod_icd10_", ukb_gwas_traits), paste0("cnt_icd10_", ukb_gwas_traits))
    missing <- setdiff(required, names(dat))
    if (length(missing)) stop("Missing GWAS inputs: ", paste(missing, collapse = ", "),
        ". Rebuild the corresponding vip_phe/hla_pca/icd/merge stages first.", call. = FALSE)
    dat$batch <- ukb_gwas_batch(dat$gen_batch)
    for (trait in ukb_gwas_traits) {
        # Preserve the existing CVD definitions. New traits use their own
        # evidence counts, exactly as t2dm does, including undated diagnoses.
        domain <- if (trait %in% c("cvd_cad", "stroke", "cvd_stroke_i", "stroke_o"))
            "cvd" else paste0("cnt_icd10_", trait)
        dat <- t2e(dat, domain, paste0("fod_icd10_", trait), "birth_date",
            "date_attend", "date_lost", "date_death", end_date, trait, "year")
    }
    dat
}

ukb_gwas_columns <- function(dat) {
    dat %>% mutate(center = factor(center)) %>%
        dplyr::select(FID, IID, ethnic.c, center, batch, tdi, edu.sco, age, sex,
            bmi, height, ldl, drug.lipid, drug.dm, drug.htn, bald,
            matches("^bald1|^hap"), all_of(paste0("PC", 1:10)),
            matches("^cvd_cad|^t2dm|^stroke_|t2e$"), ends_with(".adu"),
            -starts_with("happy_"), -matches("(r2e|bi2e)$"))
}
