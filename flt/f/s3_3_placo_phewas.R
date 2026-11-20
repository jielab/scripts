source(file.path(Sys.getenv("FLT_CODE"), "runtime.R"))
# PLACO lead-SNP PheWAS：先运行本脚本，再运行 s3_3_placo_phewas_plot.R。
# 1. 准备 analysis/03placo/FUMA_job787285/leadSNPs.txt（使用 rsID 列）。
# 2. GWAS Atlas：https://atlas.ctglab.nl/PheWAS
#    在 Search SNP 中逐个输入 lead rsID，下载 CSV，并命名为 rsID.csv。
#    保存到 analysis/03phewas/GWAS_Atlas，例如 rs123.csv；保留原始列名。
#    无年份限制，脚本统一筛选 P < 1e-5；缺少文件会记录为未查询。
# 3. 直接运行本脚本：自动查询 OpenGWAS、下载 Catalog、读取 Atlas，再分类汇总。
#    首次需安装 data.table、ieugwasr；OpenGWAS 需有效 JWT（见 credential_file）。
#    已有本机令牌配置会自动读取，不需要设置 PHEWAS_MODE 等环境变量。
# 4. 直接运行绘图脚本。修改 Atlas 文件或 phenotype_dictionary.tsv 后，重跑本脚本即可。
#    已完成查询自动使用缓存；只离线汇总时将 query_online 改为 FALSE。
# 不合并无依据的跨数据库记录；表型/类别统计按不同 FUMA 区域去重。

library(data.table)

# 一、路径与参数
project <- Sys.getenv("FLT_OUT")
infile <- file.path(Sys.getenv("FLT_FUMA_DIR"), "leadSNPs.txt")
outdir <- file.path(project, "03phewas")
atlas_dir <- file.path(outdir, "GWAS_Atlas")
query_online <- tolower(Sys.getenv("FLT_PHEWAS_ONLINE", "true")) == "true"    # 默认直接查询；FALSE 只读取已有缓存和 Atlas 文件
refresh <- FALSE        # TRUE 重新查询 SNP；Catalog 仍使用本地已有下载文件
run_opengwas <- TRUE
run_catalog <- TRUE
p_limit <- 1e-5
catalog_file <- file.path(outdir, "cache/catalog.zip")
catalog_url <- "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations_ontology-annotated-full.zip"
# 令牌文件内容：OPENGWAS_JWT=你的令牌；不要把令牌写入共享分析脚本。
credential_file <- file.path(Sys.getenv("USERPROFILE", unset = path.expand("~")), ".config/opengwas/.Renviron")
if (query_online && !nzchar(Sys.getenv("OPENGWAS_JWT")) && file.exists(credential_file)) readRenviron(credential_file)
if (!file.exists(infile)) stop("找不到 leadSNPs 文件：", infile)
dir.create(file.path(outdir, "cache"), recursive = TRUE, showWarnings = FALSE)
dir.create(atlas_dir, recursive = TRUE, showWarnings = FALSE)

# 表型标准化、分类及汇总：仅在分析阶段执行。
# Match the same phenotype across sources, and count each FUMA region once per term/category.
# Source with data.table loaded. Keep source-specific phenotype labels and study IDs.
classify_phewas <- function(d, dictionary) {
  text_or_na <- function(x) {
    x <- trimws(as.character(x))
    x[!nzchar(x) | is.na(x)] <- NA_character_
    x
  }
  key <- function(x) {
    x <- tolower(gsub("[[:space:]]+", " ", trimws(x)))
    x <- gsub("[\u2010-\u2015]", "-", x)
    trimws(x)
  }
  column <- function(x, col) {
    if (col %in% names(x)) as.character(x[[col]]) else rep(NA_character_, nrow(x))
  }
  raw <- copy(d)
  if (!"legacy_trait_class" %in% names(raw)) raw[, legacy_trait_class := column(d, "trait_class")]
  raw[, atlas_domain := column(d, "atlas_domain")]
  raw[, trait := trimws(as.character(trait))]
  if (anyNA(raw$trait) || any(!nzchar(raw$trait))) stop("Missing phenotype names.")
  one <- function(x) {
    x <- unique(na.omit(text_or_na(x)))
    if (length(x) == 1L) x else NA_character_
  }
  if (!nrow(raw)) {
    raw[, (c("canonical_trait", "trait_class", "classification_basis", "analysis_scope")) := character()]
    raw[, input_trait := logical()]
    return(raw)
  }
  traits <- raw[, .(legacy = one(legacy_trait_class), domain = one(atlas_domain)), by = trait]
  traits[, trait_key := key(trait)]
  # Strip only clear metadata wrappers. Preserve clinical definitions, contrasts,
  # side of the body, timepoints and measurement types in the phenotype name.
  suffix <- " [(](ukb data field [0-9_]+|firth correction|spa correction|saige|univariate|ma gwama|two-way meta|three-way meta|mtag)[)]$"
  previous <- NULL
  while (!identical(previous, traits$trait_key)) {
    previous <- traits$trait_key
    traits[, trait_key := sub(suffix, "", trait_key, perl = TRUE)]
  }
  traits[, trait_key := sub("^impedance measures *- *", "", trait_key)]
  traits[grepl("^non-cancer illness code[, ]*self-reported:", trait_key),
         trait_key := sub("^non-cancer illness code[, ]*self-reported: *", "self-reported ", trait_key)]
  traits[, canonical_key := trait_key]

  alias <- c(
    "Body mass index" = "^(bmi|body mass index( [(]bmi[)])?|adult bmi (sex-combined|females|males)|body mass index [(](female|male|sex-combined)[)])$",
    "Body weight" = "^(weight|body weight)$",
    "Standing height" = "^(height|standing height)$",
    "Body fat percentage" = "^((whole )?body fat percentage|body fat percent)$",
    "Waist-to-hip ratio" = "^(waist[- ]hip ratio|waist-to-hip ratio|waist to hip ratio|whr)$",
    "Major depressive disorder" = "^(major depression|major depressive disorder|mdd)$",
    "Broad depression" = "^(broad depression|depression [(]broad[)])$",
    "Depressive symptoms" = "^depressive symptoms$",
    "Neuroticism" = "^(neuroticism( (sum )?score)?|neurociticism)$",
    "Mood swings" = "^(experiencing mood swings|mood swings( [(]mood[)])?)$",
    "Years of education" = "^(years of schooling|years of education|educational attainment [(]years of education[)])$",
    "College or university completion" = "^(educational attainment [(]college completion[)]|qualifications: college or university degree)$",
    "Chicken liking" = "^(chicken liking|liking (for |of )?chicken)$",
    "Glycated haemoglobin (HbA1c)" = "^(hba1c( levels?)?|glycated ha?emoglobin( hba1c)?( levels?)?)$",
    "HDL cholesterol" = "^(hdl cholesterol|high[- ]density[- ]lipoprotein cholesterol)( levels?)?$",
    "LDL cholesterol" = "^(ldl cholesterol|low[- ]density[- ]lipoprotein cholesterol)( levels?)?$",
    "Forced vital capacity (FVC)" = "^(fvc|forced vital capacity|lung function [(]fvc[)])$",
    "Age at first sexual intercourse" = "^age (at first|first had) sexual intercourse$",
    "Ever had prolonged loss of interest" = "^(depression - )?ever had prolonged loss of interest in normal activities$"
  )
  for (label in names(alias)) {
    ix <- grepl(alias[[label]], traits$trait_key, perl = TRUE)
    traits[ix, canonical_key := key(label)]
  }

  # Named rules use a shared vocabulary. Specific phenotypes outrank source labels.
  patterns <- c(
    "Composite / multivariate" = "pleiotropy|multivariate analysis|^multi-trait|^first pc of|\\bbipolar and major depression status\\b|^(adhd|bipolar disorder|schizophrenia|major depressive disorder) vs |\\basthma and major depressive disorder\\b|\\bbipolar disorder or major depressive disorder\\b",
    "Administrative / participation" = "^time to (complete round|answer)$|^duration (to entering value|of questionnaire|screen displayed)|^duration to first press|questionnaires? completed|participation in .*questionnaire|acceptance of an invitation|^main speciality|^methods of (admission|discharge)|^destinations on discharge|^operation code|^operative procedures|^number of operations|^interpolated year|^intended management of patient|^home location|^any icdmain",
    "Gene expression" = "^ensg[0-9]+([.][0-9]+)?$",
    "Brain imaging" = "^(idp_|aparc-|aseg_|amygnuclei_|hippsubfield_|ba-exvivo_|brainstem_|wg_[lr]h_|a2009s |rfmri)|^volume of (left|right) (medial|central) nucleus",
    "Sleep / circadian rhythm" = "\\binsomnia\\b|\\bchronotype\\b|\\bsleep (duration|length|timing|quality|apnoea|apnea)\\b|^morning person$|\\bnap(ping|s)?\\b",
    "Psychiatric / mental health" = "depress|neurotic|neurocit|well-being|mood swings|miserableness|feeling miserable|life satisfaction|positive affect|\\badhd\\b|\\bautism\\b|schizophren|\\bbipolar\\b|post-traumatic stress|prolonged loss of interest|\\bmanic\\b|substance use disorder|\\banxiety\\b|externalizing behavio",
    "Cognition" = "\\bintelligence\\b|cognitive (ability|performance|function)|^fi[0-9]+ *:|symbol digit|digits remembered|puzzles correctly|numeric path|correct(ly)? (identify )?matches|self-reported math ability",
    "Education / socioeconomic" = "education|schooling|school performance|math class|^qualifications:|college or university|household income|deprivation|year immigrated|length of time at current address",
    "Diet / food preference" = "\\bliking\\b|^vegetarianism$|^spread type:|^(coffee|tea|beef|lamb|pork|cereal|bread|vegetable|fruit|salt|cheese|fish|poultry|chicken|cooked vegetable|added salt|dried fruit|fresh fruit)s? (consumption|intake|type:)",
    "Body size / composition" = "body mass index|\\bbmi\\b|^(body )?weight$|^standing height$|^sitting height$|^height$|\\b(fat mass|fat percentage|fat-free mass|lean mass|predicted mass|body water mass|obesity|adiposity)\\b|\\bfat ratio\\b|\\b(waist|hip) circumference\\b|waist.to.hip ratio|^impedance of |^comparative body size",
    "Glycaemia / lipids / metabolites" = "diabet|hba1c|glycated ha?emoglobin|\\bcholesterol\\b|triglyceride|lipoprotein|hyperlipid|\\bglucose\\b|\\binsulin\\b|\\blinoleic acid\\b|\\btyrosine\\b|\\bleucine\\b|\\bmannose\\b|hydroxyvitamin d|hydroxybutyrate|apolipo|\\burate\\b|phosphatidylcholin|phosphoglycerid|\\bcholines\\b|basal metabolic rate",
    "Cardiovascular" = "hypertension|antihypertensive|arterial pressure|blood pressure|\\bangina\\b|coronary|myocardial infarction|heart failure|\\bvaricose veins\\b|\\bstroke\\b",
    "Blood cells / haematology" = "white blood cell|\\bleukocyte|lymphocyte|reticulocyte|erythrocyte|platelet|red (blood )?cell|neutrophil|eosinophil|granulocyte|basophil|ha?ematocrit|ha?emoglobin|ha?ematopoiesis|\\bferritin\\b|transferrin saturation",
    "Immune / inflammation" = "baff-r|cd16-cd56|\\bcd8 on|natural killer|sj(.|<u[+]00f6>)gren|\\bsicca syndrome|\\bsle\\b|lupus|sarcoidosis|connective tissue disorder|c-reactive protein|immunological or systemic disorder|epstein-barr virus.*antibody",
    "Respiratory" = "asthma|\\bcopd\\b|forced vital capacity|\\bfvc\\b|\\bfev1\\b|cough on most days|impaired spr[i]?ometry",
    "Reproduction" = "hysterectomy|contraceptive pill|menarche|menopause|first sexual intercourse|first (live )?birth|number of children|age of first facial hair",
    "Activity / lifestyle / environment" = "physical activity|sedentary|screen time|time spent watching|plays computer games|smoking|\\balcohol\\b|\\bair pollution\\b|driv(e|ing) faster|duration of walks|log acceleration|time spent using computer|mobile phone use|solarium|sunlamp|inhaled tobacco|workplace very cold",
    "Endocrine" = "thyrotoxic|thyroid|graves' disease|\\bigf[ -]?1\\b|sex hormone.binding globulin|\\bshbg\\b",
    "Renal / urinary" = "nephropathy|estimated glomerular filtration|\\bcreatinine\\b|urinary tract",
    "Neurological" = "myasthenia gravis|\\bmigraine\\b|\\bheadache syndromes\\b|parkinson|alzheimer|dementia|epilepsy",
    "Cancer" = "\\bcancer\\b|\\bcarcinoma\\b|\\bneoplasm\\b|\\badenocarcinoma\\b",
    "Musculoskeletal" = "fibromyalgia|dorsalgia|soft tissue disorders|bone mineral density|osteoarthritis|heel broadband ultrasound|disc herniation|injuries to the knee",
    "Other clinical phenotypes" = "corneal|\\bkeratitis\\b|short-sightedness|intra-ocular pressure|bullous disorder|\\bgout\\b|\\bhernia\\b|heel broadband ultrasound|fibromyalgia|dorsalgia|\\bfrailty index\\b|haemorrhoidal disease|lichen planus",
    "Proteins / molecular traits" = "measurement$|\\bneurofilament light chain|tensin-2|trypsin-3|beta-defensin|\\bcomplement c[24]\\b|beta-crystallin|kallikrein-8|plexin-b2|serotransferrin|neurexin|^otu[0-9_]+|abundance in stool"
  )
  traits[, `:=`(trait_class = NA_character_, classification_basis = NA_character_)]
  for (name in names(patterns)) {
    ix <- is.na(traits$trait_class) & grepl(patterns[[name]], traits$canonical_key, perl = TRUE)
    traits[ix, `:=`(trait_class = name, classification_basis = "phenotype_rule")]
  }

  # Existing labels and Atlas domains are fallback hints, never separate categories.
  legacy_map <- c(
    "Anthropometry / body composition" = "Body size / composition",
    "Psychiatric / psychological" = "Psychiatric / mental health",
    "Lipids / metabolism" = "Glycaemia / lipids / metabolites",
    "Blood cells" = "Blood cells / haematology", "Education / cognition" = "Cognition",
    "Food preference / intake" = "Diet / food preference", "Sleep" = "Sleep / circadian rhythm",
    "General health / healthcare" = "Other clinical phenotypes", "Neurological" = "Neurological",
    "Cancer" = "Cancer", "Endocrine" = "Endocrine",
    "Gastrointestinal / liver" = "Gastrointestinal / liver", "Kidney / urinary" = "Renal / urinary",
    "Musculoskeletal" = "Musculoskeletal", "Skin / pigmentation" = "Other clinical phenotypes",
    "Eye / hearing" = "Other clinical phenotypes", "Ageing / mortality" = "Other clinical phenotypes",
    "Cardiovascular" = "Cardiovascular", "Respiratory" = "Respiratory", "Reproduction" = "Reproduction",
    "Gene expression" = "Gene expression", "Molecular traits / proteins" = "Proteins / molecular traits",
    "Immune / inflammatory" = "Immune / inflammation", "Socioeconomic / social" = "Education / socioeconomic",
    "Smoking" = "Activity / lifestyle / environment", "Alcohol" = "Activity / lifestyle / environment",
    "Physical activity / sedentary behaviour" = "Activity / lifestyle / environment"
  )
  atlas_map <- c(
    Psychiatric = "Psychiatric / mental health", Cognitive = "Cognition",
    Activities = "Activity / lifestyle / environment", Environment = "Activity / lifestyle / environment",
    Environmental = "Activity / lifestyle / environment", Metabolic = "Glycaemia / lipids / metabolites",
    Nutritional = "Diet / food preference", Immunological = "Immune / inflammation",
    Respiratory = "Respiratory", Reproduction = "Reproduction", Cardiovascular = "Cardiovascular",
    `Social Interactions` = "Education / socioeconomic", `Body Structures` = "Body size / composition",
    Skeletal = "Musculoskeletal", Muscular = "Musculoskeletal",
    `Connective Tissue` = "Other clinical phenotypes", Gastrointestinal = "Gastrointestinal / liver",
    Dermatological = "Other clinical phenotypes", Neurological = "Neurological",
    Endocrine = "Endocrine", Mortality = "Other clinical phenotypes"
  )
  traits[, domain := fcoalesce(domain, fifelse(startsWith(legacy, "Atlas: "), sub("^Atlas: ", "", legacy), NA_character_))]
  traits[is.na(trait_class) & !is.na(legacy_map[legacy]),
         `:=`(trait_class = unname(legacy_map[legacy]), classification_basis = "legacy_class_crosswalk")]
  traits[is.na(trait_class) & !is.na(atlas_map[domain]),
         `:=`(trait_class = unname(atlas_map[domain]), classification_basis = "atlas_domain_crosswalk")]

  # Reuse the best resolved category for equivalent names; flag disagreements.
  traits[, priority := match(classification_basis, c("phenotype_rule", "legacy_class_crosswalk", "atlas_domain_crosswalk"))]
  best <- traits[!is.na(priority)][order(priority), {
    same_rank <- unique(trait_class[priority == priority[1L]])
    .(resolved = if (length(same_rank) == 1L) same_rank else NA_character_,
      basis = if (length(same_rank) == 1L) classification_basis[1L] else "category_conflict")
  }, by = canonical_key]
  traits[best, on = "canonical_key", `:=`(trait_class = i.resolved, classification_basis = i.basis)]
  display <- setNames(names(alias), key(names(alias)))
  traits[, canonical_trait := fcoalesce(unname(display[canonical_key]), canonical_key)]
  traits[, manual_class := FALSE]
  if (file.exists(dictionary)) {
    dict <- fread(dictionary, colClasses = "character", na.strings = c("", "NA"), encoding = "UTF-8")
    stopifnot(all(c("trait", "trait_class") %in% names(dict)))
    dict[, `:=`(trait = trimws(trait), trait_class = text_or_na(trait_class))]
    if (!"canonical_trait" %in% names(dict)) dict[, canonical_trait := NA_character_]
    dict[, canonical_trait := text_or_na(canonical_trait)]
    if (anyDuplicated(dict$trait)) stop("phenotype_dictionary.tsv contains duplicate trait names.")
    message("Manual dictionary: ", dictionary, "; rows: ", nrow(dict),
            "; matched raw terms: ", sum(dict$trait %chin% traits$trait))
    traits[dict, on = "trait", `:=`(
      canonical_trait = fcoalesce(i.canonical_trait, canonical_trait),
      trait_class = fcoalesce(i.trait_class, trait_class), manual_class = !is.na(i.trait_class),
      classification_basis = fifelse(!is.na(i.trait_class) | !is.na(i.canonical_trait),
                                     "manual_dictionary", classification_basis))]
  }
  traits[, canonical_key := key(canonical_trait)]
  # Manual decisions apply to all aliases of the same final concept. Different
  # classes assigned manually to one concept must be resolved in the dictionary.
  decisions <- traits[!is.na(trait_class), {
    manual <- unique(trait_class[manual_class])
    all_classes <- unique(trait_class)
    .(resolved = if (length(manual) == 1L) manual else if (length(all_classes) == 1L) all_classes else NA_character_,
      manual_conflict = length(manual) > 1L)
  }, by = canonical_key]
  if (any(decisions$manual_conflict)) stop("Conflicting manual classes assigned to one canonical trait.")
  traits[decisions, on = "canonical_key", `:=`(
    classification_basis = fifelse(is.na(i.resolved), "category_conflict",
      fifelse(is.na(trait_class) | trait_class != i.resolved, "canonical_trait_consensus", classification_basis)),
    trait_class = i.resolved)]
  # Use consistent spelling for all aliases of a concept.
  labels <- traits[order(classification_basis != "manual_dictionary", canonical_trait),
                   .(label = canonical_trait[1L]), by = canonical_key]
  traits[labels, on = "canonical_key", canonical_trait := i.label]
  traits[is.na(trait_class), trait_class := "Unclassified"]
  traits[is.na(classification_basis), classification_basis := "unclassified"]
  traits[, analysis_scope := fcase(
    trait_class == "Unclassified", "unclassified",
    trait_class == "Administrative / participation", "administrative",
    trait_class == "Composite / multivariate", "composite",
    trait_class %chin% c("Gene expression", "Proteins / molecular traits"), "molecular",
    default = "phenotype")]
  for (col in intersect(names(raw), c("canonical_trait", "canonical_key", "trait_class", "classification_basis", "analysis_scope", "input_trait")))
    raw[, (col) := NULL]
  raw <- merge(raw, traits[, .(trait, canonical_trait, trait_class, classification_basis, analysis_scope)],
               by = "trait", all.x = TRUE, sort = FALSE)
  raw[, input_trait := pair == "G732__mdd_2025" &
        key(canonical_trait) %chin% c("chicken liking", "major depressive disorder")]
  raw
}

summarise_phewas <- function(d, write_tsv) {
  trait_counts <- d[, .(n_records = .N, n_leads = uniqueN(lead_snp),
                        n_loci = uniqueN(locus), n_sources = uniqueN(source),
                        n_raw_terms = uniqueN(trait), min_p = min(p)),
                    by = .(pair, canonical_trait, trait_class, analysis_scope, input_trait)]
  class_counts <- d[, .(n_records = .N, n_leads = uniqueN(lead_snp),
                        n_loci = uniqueN(locus), n_traits = uniqueN(canonical_trait),
                        n_sources = uniqueN(source)),
                    by = .(pair, trait_class, analysis_scope)]
  audit <- d[, .(n_records = .N, n_raw_terms = uniqueN(trait),
                 n_canonical_traits = uniqueN(canonical_trait), n_loci = uniqueN(locus)),
             by = .(pair, source, analysis_scope)]
  write_tsv(trait_counts[order(pair, -n_loci, canonical_trait)], "trait_counts.tsv")
  write_tsv(class_counts[order(pair, -n_loci, trait_class)], "trait_class_locus_counts.tsv")
  write_tsv(audit, "classification_summary.tsv")
  write_tsv(d[analysis_scope == "unclassified",
             .(n_loci = uniqueN(locus), n_records = .N, sources = paste(sort(unique(source)), collapse = ";")),
             by = .(trait, canonical_trait, classification_basis)][order(-n_loci, -n_records, trait)],
            "classification_needs_review.tsv")
  write_tsv(unique(d[, .(trait, canonical_trait, trait_class, classification_basis, analysis_scope)]),
            "phenotype_classification_auto.tsv")
  invisible(list(traits = trait_counts, classes = class_counts))
}


# 2. 辅助函数：三种来源统一成相同列，缓存只承认成功完成的查询。
num <- function(x) suppressWarnings(as.numeric(x))
# 主结果独立保存，中间查询/导入核查表集中保存，减少零散文件。
run_audit <- list()
main_tables <- c("query_snps.tsv", "lead_phewas_results.tsv", "source_coverage.tsv", "phenotype_classification_auto.tsv",
                 "classification_needs_review.tsv", "classification_summary.tsv", "trait_counts.tsv", "trait_class_locus_counts.tsv")
write_tsv <- function(x, name) {
  if (name %in% main_tables || name == "phenotype_dictionary.tsv")
    fwrite(x, file.path(outdir, name), sep="\t", na="", bom=TRUE)
  else run_audit[[name]] <<- copy(x)
}
init_tsv <- function(x, name) if (!file.exists(file.path(outdir, name))) write_tsv(x, name)
field <- function(d, ...) {
  nm <- intersect(c(...), names(d))
  if (length(nm)) d[[nm[1]]] else rep(NA_character_, nrow(d))
}
empty_hits <- function() data.table(query_snp = character(), source = character(),
  study_id = character(), trait = character(), p = double(), beta = double(),
  se = double(), n = character(), ea = character(), nea = character(), record_id = character())
normalise <- function(d, source, expected_snp = NULL) {
  d <- as.data.table(d)
  if (!nrow(d)) return(empty_hits())
  h <- data.table(query_snp = as.character(field(d, "query_snp", "rsid", "SNP", "variant")),
    source = source, study_id = as.character(field(d, "study_id", "id", "STUDY ACCESSION")),
    trait = as.character(field(d, "trait", "DISEASE/TRAIT")), p = num(field(d, "p", "P-VALUE")),
    beta = num(field(d, "beta")), se = num(field(d, "se")),
    n = as.character(field(d, "n", "sample_size", "INITIAL SAMPLE SIZE")),
    ea = as.character(field(d, "ea", "effect_allele")),
    nea = as.character(field(d, "nea", "other_allele")),
    record_id = as.character(field(d, "record_id", "study_id", "id")))
  if (anyNA(h$query_snp) || any(!grepl("^rs[0-9]+$", h$query_snp)) ||
      any(!is.finite(h$p) | h$p < 0 | h$p > 1) || anyNA(h$trait) || any(!nzchar(h$trait)))
    stop(source, "：返回数据缺少有效 SNP、表型或 P 值。")
  if (!is.null(expected_snp) && any(h$query_snp != expected_snp)) stop("查询 SNP 与返回 SNP 不一致。")
  unique(h[p < p_limit])
}
cache_path <- function(source, snp) file.path(outdir, "cache", paste0(source, "_", snp, ".rds"))
read_cache <- function(source, snp) {
  f <- cache_path(source, snp)
  if (!file.exists(f)) return(NULL)
  tryCatch({
    z <- suppressWarnings(readRDS(f))
    if (!isTRUE(z$complete) || !identical(z$p_limit, p_limit) || !is.data.frame(z$data) ||
        !all(names(empty_hits()) %in% names(z$data))) return(NULL)
    normalise(z$data, source, snp)
  }, error = function(e) NULL)
}
save_cache <- function(source, snp, data, provenance) {
  f <- cache_path(source, snp); tmp <- paste0(f, ".tmp")
  saveRDS(list(data = data, complete = TRUE, p_limit = p_limit,
    queried_at = as.character(Sys.time()), provenance = provenance), tmp)
  if (file.exists(f)) unlink(f)
  if (!file.rename(tmp, f)) stop("缓存保存失败：", f)
}
pending <- function(source) query_snps[refresh | vapply(query_snps,
  function(s) is.null(read_cache(source, s)), logical(1))]

# 3. 只查询本次 FUMA lead SNP；不自动用 IndSigSNP 替代。
raw_lead <- fread(infile)
stopifnot(all(c("rsID", "GenomicLocus", "chr", "pos", "p", "uniqID") %in% names(raw_lead)))
lead <- raw_lead[, .(pair = "G732__mdd_2025", locus = as.character(GenomicLocus), snp = rsID, chromosome = chr, pos = pos, uniqID = uniqID, p_placo = num(p))]
stopifnot(nrow(lead) > 0, !anyDuplicated(lead$snp), all(grepl("^rs[0-9]+$", lead$snp)), all(is.finite(lead$p_placo)))
query_snps <- sort(unique(lead$snp))
init_tsv(data.table(trait = character(), canonical_trait = character(), trait_class = character()), "phenotype_dictionary.tsv")
write_tsv(lead, "lead.tsv")
write_tsv(data.table(query_snp = query_snps), "query_snps.tsv")
message("本次 lead SNP 数量：", length(query_snps))

# 4. 查询 OpenGWAS。仅查询未完成/损坏的缓存；认证和额度错误立即停止。
if (query_online && run_opengwas && length(todo <- pending("OpenGWAS"))) {
  if (!requireNamespace("ieugwasr", quietly = TRUE)) stop("请先安装 ieugwasr。")
  if (!nzchar(ieugwasr::get_opengwas_jwt())) stop("缺少 OpenGWAS 令牌，请在 credential_file 指定的文件中配置 OPENGWAS_JWT。")
  for (s in todo) tryCatch({
    message("OpenGWAS：", s)
    d <- ieugwasr::phewas(variants = s, pval = p_limit)
    if (!is.data.frame(d) || nrow(d) >= 100000L) stop("返回格式异常或达到 API 行数上限。")
    if (nrow(d) && !all(c("p", "trait", "id") %in% names(d))) stop("OpenGWAS 返回列不完整。")
    # PheWAS may return other rsIDs at the queried position. Keep exact IDs only.
    if (nrow(d) && "rsid" %in% names(d)) {
      other <- d[is.na(d$rsid) | d$rsid != s, , drop = FALSE]
      if (nrow(other)) write_tsv(as.data.table(other), paste0("OpenGWAS_nonmatching_", s, ".tsv"))
      d <- d[!is.na(d$rsid) & d$rsid == s, , drop = FALSE]
    }
    h <- normalise(d, "OpenGWAS", s)
    save_cache("OpenGWAS", s, h, paste("ieugwasr", packageVersion("ieugwasr")))
    Sys.sleep(1)
  }, error = function(e) {
    message(s, " 查询失败，不记作零结果：", conditionMessage(e))
    if (grepl("429|allowance|401|403|token|JWT", conditionMessage(e), ignore.case = TRUE)) stop(e)
  })
}

# 5. Catalog：先筛出目标 SNP 的记录，再拆分 rsID，避免逐行遍历全量数据库。
if (query_online && run_catalog && length(todo <- pending("GWASCatalog"))) {
  if (!file.exists(catalog_file)) {
    dir.create(dirname(catalog_file), recursive = TRUE, showWarnings = FALSE)
    old <- getOption("timeout"); options(timeout = max(1800, old))
    tryCatch(download.file(catalog_url, paste0(catalog_file, ".part"), mode = "wb", method = "libcurl"), finally = options(timeout = old))
    if (!file.rename(paste0(catalog_file, ".part"), catalog_file)) stop("Catalog 下载保存失败。")
  }
  tab <- catalog_file
  if (grepl("\\.zip$", tab, ignore.case = TRUE)) {
    members <- unzip(tab, list = TRUE)$Name
    members <- members[grepl("\\.(tsv|txt)$", members, ignore.case = TRUE)]
    if (length(members) != 1L) stop("请将 catalog_file 指向解压后的 tsv。")
    tab <- unzip(tab, files = members, exdir = file.path(outdir, "cache", "catalog"))
  }
  hdr <- names(fread(tab, nrows = 0)); required <- c("SNPS", "P-VALUE", "DISEASE/TRAIT", "STUDY ACCESSION")
  stopifnot(all(required %in% hdr))
  d <- fread(tab, select = intersect(c(required, "INITIAL SAMPLE SIZE", "SNP_ID_CURRENT"), hdr), colClasses = "character")
  d[, `:=`(p = num(`P-VALUE`), record_id = paste0("catalog_row_", seq_len(.N)))]
  d <- d[is.finite(p) & p >= 0 & p < p_limit & !grepl("\\s[xX]\\s", SNPS)]
  current_id <- as.character(field(d, "SNP_ID_CURRENT"))
  current_id <- ifelse(grepl("^[0-9]+$", current_id), paste0("rs", current_id), current_id)
  keep <- grepl(paste0("\\b(", paste(todo, collapse = "|"), ")\\b"), d$SNPS) | current_id %chin% todo
  d <- d[keep]; current_id <- current_id[keep]
  ids <- Map(function(a, b) intersect(unique(c(a, b)), todo), regmatches(d$SNPS, gregexpr("rs[0-9]+", d$SNPS)), current_id)
  ix <- data.table(query_snp = as.character(unlist(ids)), row_id = rep(seq_len(nrow(d)), lengths(ids)))
  h <- if (nrow(ix)) normalise(cbind(ix[, .(query_snp)], d[ix$row_id]), "GWASCatalog") else empty_hits()
  provenance <- paste(basename(catalog_file), "MD5", unname(tools::md5sum(catalog_file)))
  for (s in todo) save_cache("GWASCatalog", s, h[query_snp == s], provenance)
}

# 6. 汇总 lead SNP；查询失败/未查询不记为零关联。
# 6. 读取 Atlas、合并三来源并保存分类结果
pieces <- list(empty_hits()); cov <- list()
for (source in c("OpenGWAS", "GWASCatalog")) for (s in query_snps) {
  h <- read_cache(source, s)
  if (!is.null(h)) pieces[[length(pieces) + 1L]] <- h
  cov[[length(cov) + 1L]] <- data.table(query_snp = s, source = source, status = if (is.null(h)) "not_complete" else "complete")
}
# 文件名必须为 rsID.csv；导出文件本身无 SNP 列，身份依据人工命名。
files <- list.files(atlas_dir, pattern = "[.]csv$", full.names = TRUE, ignore.case = TRUE)
ids <- sub("[.]csv$", "", basename(files), ignore.case = TRUE)
if (any(!grepl("^rs[0-9]+$", ids)) || anyDuplicated(ids) || any(!ids %chin% query_snps)) stop("Atlas 文件需按 lead rsID.csv 唯一命名，且必须属于查询清单。")
atlas_parts <- list(empty_hits()); metadata <- list(); audits <- list()
ac <- data.table(query_snp=query_snps, status="not_queried")
for (i in seq_along(files)) {
  z <- fread(files[i], encoding="UTF-8", colClasses="character")
  required_atlas <- c("atlas ID", "PMID", "Year", "Domain", "Trait", "P-value", "N", "EA", "NEA")
  if (!all(required_atlas %in% names(z))) stop("Atlas CSV 列不完整：", files[i])
  z[, p_num := num(`P-value`)]
  if (any(!is.finite(z$p_num) | z$p_num < 0 | z$p_num > 1) || anyNA(z$Trait) || any(!nzchar(z$Trait))) stop("Atlas CSV 数据无效：", files[i])
  h <- z[, .(query_snp=ids[i], study_id=paste0("atlas:", `atlas ID`), trait=Trait, p=p_num, n=N, ea=EA, nea=NEA, record_id=paste0("atlas:", `atlas ID`))]
  atlas_parts[[length(atlas_parts)+1L]] <- normalise(h, "GWASAtlas", ids[i])
  metadata[[length(metadata)+1L]] <- z[, .(query_snp=ids[i], study_id=paste0("atlas:", `atlas ID`), trait=Trait, atlas_domain=Domain, pmid=PMID, year=Year)]
  audits[[length(audits)+1L]] <- data.table(query_snp=ids[i], file=basename(files[i]), md5=unname(tools::md5sum(files[i])), n_raw=nrow(z), n_pass=sum(z$p_num < p_limit), year_filter="none", identity_basis="filename supplied by user")
  ac[query_snp == ids[i], status := "complete"]
}
atlas <- unique(rbindlist(atlas_parts, fill=TRUE))
atlas_meta <- if (length(metadata)) unique(rbindlist(metadata)) else
  data.table(query_snp=character(), study_id=character(), trait=character(), atlas_domain=character(), pmid=character(), year=character())
write_tsv(atlas, "atlas_import.tsv")
write_tsv(ac, "atlas_coverage.tsv")
write_tsv(atlas_meta, "atlas_metadata.tsv")
if (length(audits)) write_tsv(rbindlist(audits), "atlas_file_audit.tsv")
pieces[[length(pieces) + 1L]] <- normalise(atlas, "GWASAtlas")
cov[[length(cov) + 1L]] <- merge(data.table(query_snp = query_snps), ac, all.x = TRUE, by = "query_snp")[, source := "GWASAtlas"]
coverage <- rbindlist(cov, use.names = TRUE); coverage[is.na(status), status := "not_complete"]
complete <- coverage[, .(ok = .N == 3L && all(status == "complete")), by = query_snp][ok == TRUE, query_snp]
hits <- unique(rbindlist(pieces)); observed <- unique(hits$query_snp)
decisions <- lead[, .(pair, locus, lead_snp = snp, query_snp = snp)]
decisions[, decision := fifelse(query_snp %chin% complete, fifelse(query_snp %chin% observed, "complete_with_hits", "complete_no_hit"), "partial_sources")]
result <- merge(lead, hits, by.x = "snp", by.y = "query_snp", allow.cartesian = TRUE)
setnames(result, "snp", "query_snp")
result[, `:=`(lead_snp = query_snp, is_fallback = FALSE, coordinate_build = "GRCh37_FUMA")]
atlas_meta[, source := "GWASAtlas"]
result <- merge(result, atlas_meta, by = c("query_snp", "study_id", "trait", "source"), all.x = TRUE)
# Reuse available legacy category hints when rebuilding a previously classified report.
old_report <- file.path(outdir, "cache", "classification_reference.tsv")
if (!file.exists(old_report)) old_report <- file.path(outdir, "lead_phewas_results.tsv")
if (file.exists(old_report)) {
  old_header <- names(fread(old_report, nrows = 0L))
  old_class <- intersect(c("legacy_trait_class", "trait_class"), old_header)
  if (length(old_class) && "trait" %in% old_header) {
    prior <- fread(old_report, select = c("trait", old_class[1L]), colClasses = "character", na.strings = c("", "NA"))
    setnames(prior, old_class[1L], "legacy_trait_class")
    prior[, trait := trimws(trait)]
    prior <- unique(prior[!is.na(legacy_trait_class) & nzchar(trimws(legacy_trait_class))])
    prior <- prior[!duplicated(trait) & !duplicated(trait, fromLast = TRUE)]
    result[, trait := trimws(trait)]
    result[prior, on = "trait", legacy_trait_class := i.legacy_trait_class]
  }
}
result <- classify_phewas(result, file.path(outdir, "phenotype_dictionary.tsv"))
summarise_phewas(result, write_tsv)
write_tsv(coverage, "source_coverage.tsv")
write_tsv(decisions, "representative_decisions.tsv")
write_tsv(hits, "current_all_sources.tsv")
write_tsv(result, "lead_phewas_results.tsv")
writeLines(c(paste("Report:", Sys.time()), paste("All-source complete SNPs:", length(complete), "/", length(query_snps)),
  "Lead SNPs only. Coordinates are FUMA GRCh37; queries use rsID. Effects remain source-reported.",
  "Association records are retained by source. Phenotype and category summaries count distinct FUMA loci."), file.path(outdir, "run_notes.txt"))
message("当前关联记录：", nrow(result), "；三来源全部完成：", length(complete), "/", length(query_snps))

saveRDS(run_audit, file.path(outdir, "cache", "run_audit.rds"))
message("输出目录：", normalizePath(outdir, winslash = "/"))
