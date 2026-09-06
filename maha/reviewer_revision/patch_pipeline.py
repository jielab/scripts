from pathlib import Path
root = Path(__file__).resolve().parents[1]
p = root / 'f/MAHA_ukb.R'
s = p.read_text(encoding='utf-8')
def sub(a,b):
    global s
    assert a in s, a[:120]
    s = s.replace(a,b)
sub('source(file.path(.this_dir, "comm.f.R"))', 'source(file.path(.this_dir, "comm.f.R"))\nsource(file.path(.this_dir, "ukb_followup.R"))')
sub('covs <- c("age",', 'covs <- c("age.landmark",')
sub('diet.inc <- c(maha_ref,', 'vars.risk <- c("age.landmark", "sex", "tdi", "PC1", "PC2")\n\ndiet.inc <- c(maha_ref,')
a=s.index('for (Y in names(dx.lst)) {', s.index('maha_run_step("data_qc"'))
b=s.index('\ndat <- dat %>% mutate',a)
s=s[:a]+'''# End of the entire recorded WebQ window also covers recalls whose food
# values may have survived occasion-level QC in the existing score builder.
webq_raw <- readRDS(ukb_required_inputs[["diet_raw_rds"]])
webq_window <- maha_webq_window(webq_raw)
rm(webq_raw); invisible(gc())
dat$eid <- as.character(dat$eid)
dat <- maha_set_landmark(dat, webq_window, diet.inc.s100)
ukb_time_audit <- list()
for (Y in names(dx.lst)) {
  if (Y != "death" && !paste0("fod_icd10_", Y) %in% names(dat)) next
  dat[grep(paste0("^", Y, "\\\\.([Y]?(t2e|r2e)|b2e|bi2e)$"), names(dat))] <- NULL
  dat <- maha_outcome_followup(dat, Y, date_follow_end)
  ukb_time_audit[[Y]] <- attr(dat, "maha_outcome_audit")
}
write_xlsx(bind_rows(ukb_time_audit), "Reviewer.followup_audit.xlsx")
writeLines(c(maha_followup_version,
  "Survival time origin: max(enrolment, last of all recorded WebQ dates).",
  "Outcome-specific events on/before landmark excluded; positive follow-up required.",
  "Cox age: enrolment age plus elapsed time to landmark. Other covariates remain enrolment measurements.",
  "Existing diet scores unchanged; scoring is not optimized using outcomes.",
  paste("Administrative censoring:", date_follow_end)), cohort_file("Reviewer.methods.txt"))
''' +s[b:]
sub('comp <- comp %>% dplyr::select(eid, all_of(ukb_component_cols))','comp <- comp %>% mutate(eid = as.character(eid)) %>% dplyr::select(eid, all_of(ukb_component_cols))')
sub('ukb_validation <- ukb_validate_maha(dat)', '''dat$maha_n_components <- rowSums(!is.na(dat[ukb_component_cols]))
component_audit <- bind_rows(lapply(Y.inc, function(y) {
  cc <- complete.cases(dat[, c(paste0(y, c(".t2e", ".Yt2e")), covs, diet.inc.pts)])
  data.frame(Outcome = y, N = sum(cc), N_complete9 = sum(cc & dat$maha_n_components == 9),
    proportion_complete9 = if (sum(cc)) mean(dat$maha_n_components[cc] == 9) else NA_real_)
}))
write_xlsx(list(complete_case_survival = component_audit,
  score_sample = as.data.frame(table(dat$maha_n_components[!is.na(dat[[maha_sum]])]))),
  "Reviewer.component_completeness.xlsx")
ukb_validation <- ukb_validate_maha(dat)''')
a=s.index('phewas.rds <- cohort_file('); b=s.index('\nlab <- setNames',a)
s=s[:a]+'''phewas.res <- maha_raw_phewas(dat1, Xs, vars.basic,
  file.path(indir, "Rdata", "PheWAS.rds"), cohort_file("phewas_raw_v1"),
  replace = any(c("all", "fig1") %in% strsplit(Sys.getenv("MAHA_REPLACE_STEPS"), ",", fixed = TRUE)[[1]]))
writeLines(c("Raw P preserved; BH FDR computed separately for each score.",
  "Main and supplementary PheWAS plots cap -log10(P) at 50 for display only.",
  "P=0 is retained in p_raw and flagged; logp uses .Machine$double.xmin for underflow.",
  paste("Input signature:", phewas.res$signature)), cohort_file("Reviewer.phewas_methods.txt"))
''' +s[b:]
sub('\t\tlogp = -log10(p),\n\t\tFDR = p.adjust(p, "BH"),','\t\tlogp = logp,')
sub('filter(!is.na(Diet), is.finite(p), p > 0)', 'filter(!is.na(Diet), is.finite(p), p >= 0)')
a=s.index('phewas_xmap <-'); b=s.index('\nphewas_centers <-',a)
s=s[:a]+'''phewas_xmap <- res %>%
  distinct(category, phenotype) %>%
  arrange(category, suppressWarnings(as.numeric(phenotype)), phenotype) %>%
  mutate(phe_x = row_number())
''' +s[b:]
sub('aes(phe_x, logp, color = category)', 'aes(phe_x, logp_display, color = category)')
sub('y = expression(-log[10](italic(p))))', 'y = "-log10(P); display cap 50")')
sub('logp.MAHA.cap = pmin(logp.MAHA, 15), logp.DASH.cap = pmin(logp.DASH, 15)', 'logp.MAHA.cap = pmin(logp.MAHA, 50), logp.DASH.cap = pmin(logp.DASH, 50)')
sub('xlim = c(0, 16), ylim = c(0, 16)', 'xlim = c(0, 52), ylim = c(0, 52)')
sub('seq(0, 15, 3)', 'seq(0, 50, 10)')
sub('labs(x = "DASH: -log10(P)", y = "MAHA: -log10(P)", title = "a. Association strength")', 'labs(x = "DASH: -log10(P); display cap 50", y = "MAHA: -log10(P); display cap 50", title = "a. Association significance")')
# Only risk model calls use the landmark age. Cross-sectional PheWAS keeps baseline age.
lines=s.splitlines()
for i,line in enumerate(lines):
    if any(x in line for x in ['extract_joint_risk(dat1', 'plot_risk(dat1', 'plot_cuminc_joint(dat1']):
        lines[i]=line.replace('vars.basic','vars.risk')
s='\n'.join(lines)+'\n'
# Geography retains baseline age for descriptive regressions, landmark age for Cox.
sub('setdiff(covs, c("center", "age", "sex.f", "tdi"))', 'setdiff(covs, c("center", "age.landmark", "sex.f", "tdi"))')
sub('score_z, age, sex.f, tdi,', 'score_z, age.landmark, sex.f, tdi,')
sub('~ score_z + age + sex.f + tdi', '~ score_z + age.landmark + sex.f + tdi')
# Append model metadata to absolute-risk data, without changing the fitted model.
sub('\tres\n}\n\nap_effect_from_assoc', '''\tres$N_model <- nrow(d)
  res$Events_model <- sum(d[[Y_event]])
  res$adjustment <- paste(covs2, collapse = "; ")
  res$model <- "Cox; DASH + MAHA (no interaction)"
  res$prediction <- "Fixed profile: numeric means; categorical modes in model sample"
  res$time_origin <- "End of all recorded WebQ assessments (or enrolment if later)"
  res$N_complete9 <- sum(d$maha_n_components == 9)
  res
}

ap_effect_from_assoc''')
sub('risk10 = round(100 * (1 - ss$surv), 2)', '''risk10 = 100 * (1 - ss$surv),
      N_model = nrow(d2), Events_model = sum(d2$event),
      N_group = sum(d2$dash == nd$dash[i] & d2$maha == nd$maha[i]),
      Events_group = sum(d2$event[d2$dash == nd$dash[i] & d2$maha == nd$maha[i]]),
      adjustment = paste(covs, collapse = "; "), model = "Cox; DASH * MAHA interaction",
      prediction = "Fixed profile: numeric means; categorical modes in model sample",
      time_origin = "End of all recorded WebQ assessments (or enrolment if later)"''')
sub('dplyr::select(Outcome, dash, maha, risk10)', 'dplyr::select(Outcome, dash, maha, risk10, everything())')
sub('ukb_construct <- ukb_construct_profile(dat)', '''ukb_construct <- ukb_construct_profile(dat)
group_audit <- bind_rows(lapply(diet.inc, function(nm) {
  dat %>% filter(!is.na(.data[[paste0("diet.", nm, ".3c")]])) %>%
    group_by(group = .data[[paste0("diet.", nm, ".3c")]]) %>%
    summarise(Diet = nm, N = n(), minimum_raw_score = min(.data[[paste0("diet.", nm, ".sum")]], na.rm = TRUE),
      maximum_raw_score = max(.data[[paste0("diet.", nm, ".sum")]], na.rm = TRUE), .groups = "drop")
}))
group_audit$rule <- "Unweighted rank quartiles among White participants with each score: Q1 low; Q2-Q3 middle; Q4 high; ties may split"
write_xlsx(group_audit, "Reviewer.group_cutpoints.xlsx")''')
p.write_text(s,encoding='utf-8')
p=root/'f/MAHA_publication.R'; s=p.read_text(encoding='utf-8')
s=s.replace('if (scope == "all") make_shared_score_figure()', '''shared_inputs <- c(aux_path("ukb", "FigS3.score_distribution_concordance.out.xlsx"),
  aux_path("nhanes", "FigS3.score_distribution_concordance.out.xlsx"),
  aux_path("chns", "FigS7.score_distribution_concordance.out.xlsx"))
if (all(file.exists(shared_inputs))) {
  make_shared_score_figure()
} else if (scope == "all") {
  stop("Shared score figure requires all three cohort outputs: ", paste(shared_inputs[!file.exists(shared_inputs)], collapse = "; "))
} else message("[publication] Shared score figure pending until all three cohorts have been run.")''')
p.write_text(s,encoding='utf-8')
