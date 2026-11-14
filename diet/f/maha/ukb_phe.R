# MAHA manuscript engine, reorganized from jielab/scripts 971a61b44b84500ab46e2f5592478505a0eae50e.
step_header("Data QC")
# Read and release the large recall file before loading the main phenotype table.
webq_raw <- readRDS(file.path(Sys.getenv("DIET_UKB_RDATA"), "diet0.rds"))
webq_window <- maha_webq_window(webq_raw)
rm(webq_raw); invisible(gc())
dat0 <- readRDS(file.path(indir, "Rdata", "all.rds"))
dat0 <- maha_join_scores(dat0, "eid")

dat <- dat0 %>%
	filter(ethnic.c == "White") %>%
	mutate(
		across(all_of(diet.sum.cols), std_10_90, .names = "{sub('\\\\.sum$', '', .col)}.pts"),
		across(all_of(diet.sum.cols), score_0_100, .names = "{sub('\\\\.sum$', '', .col)}.s100"),
		across(all_of(diet.sum.cols), score_q5, .names = "{sub('\\\\.sum$', '', .col)}.q5"),
		across(all_of(diet.sum.cols), score_0_100_q5, .names = "{sub('\\\\.sum$', '', .col)}.s100_q5"),
		across(all_of(diet.sum.cols), f3c, .names = "{sub('\\\\.sum$', '', .col)}.3c")
	)

# End of the entire recorded WebQ window also covers recalls whose food
# values may have survived occasion-level QC in the existing score builder.
rm(dat0); invisible(gc())
dat$eid <- as.character(dat$eid)
dat <- maha_set_landmark(dat, webq_window, diet.inc.s100)
ukb_time_audit <- list()
for (Y in Y.inc) {
	dat[grep(paste0("^", Y, "\\.([Y]?(t2e|r2e)|b2e|bi2e)$"), names(dat))] <- NULL
	dat <- maha_outcome_followup(dat, Y, date_follow_end)
	ukb_time_audit[[Y]] <- attr(dat, "maha_outcome_audit")
}
write_xlsx(bind_rows(ukb_time_audit), "Reviewer.followup_audit.xlsx")
writeLines(c(maha_followup_version,
	"Survival time origin: max(enrolment, last of all recorded WebQ dates).",
	"Outcome-specific events on/before landmark excluded; positive follow-up required.",
	"Missing/undated WebQ records excluded from survival analysis, retained in descriptive analyses.",
	"Cox age: enrolment age plus elapsed time to landmark. Other covariates remain enrolment measurements.",
	"Existing diet scores unchanged; scoring is not optimized using outcomes.",
	paste("Administrative censoring:", date_follow_end)), cohort_file("Reviewer.methods.txt"))

dat <- dat %>% mutate(across(any_of(diet.inc.pts), ~ mk_hml(., group_pct_th), .names = "{sub('\\\\.pts$', '', .col)}.hml"))

if (!all(ukb_component_cols %in% names(dat))) {
	comp <- readRDS(file.path(Sys.getenv("DIET_UKB_RDATA"), "diet.maha.rds")) %>% as.data.frame()
	if (!"eid" %in% names(comp)) stop("UKB diet.maha.rds does not contain eid.", call. = FALSE)
	miss_comp <- setdiff(ukb_component_cols, names(comp))
	if (length(miss_comp)) stop("UKB diet.maha.rds missing MAHA component columns: ", paste(miss_comp, collapse = ", "), call. = FALSE)
	comp <- comp %>% mutate(eid = as.character(eid)) %>% dplyr::select(eid, all_of(ukb_component_cols))
	dat <- dat %>% dplyr::select( - any_of(ukb_component_cols)) %>% left_join(comp, by = "eid")
}

dat$maha_n_components <- rowSums(!is.na(dat[ukb_component_cols]))
component_audit <- bind_rows(lapply(Y.inc, function(y) {
	cc <- complete.cases(dat[, c(paste0(y, c(".t2e", ".Yt2e")), covs, diet.inc.pts)])
	data.frame(Outcome = y, N = sum(cc), N_complete9 = sum(cc & dat$maha_n_components == 9),
		proportion_complete9 = if (sum(cc)) mean(dat$maha_n_components[cc] == 9) else NA_real_)
}))
write_xlsx(list(complete_case_survival = component_audit,
	score_sample = as.data.frame(table(dat$maha_n_components[!is.na(dat[[maha_sum]])]))),
	"Reviewer.component_completeness.xlsx")
group_audit <- bind_rows(lapply(diet.inc, function(nm) {
	dat %>% filter(!is.na(.data[[paste0("diet.", nm, ".3c")]])) %>%
		group_by(group = .data[[paste0("diet.", nm, ".3c")]]) %>%
		summarise(Diet = nm, N = n(), minimum_raw_score = min(.data[[paste0("diet.", nm, ".sum")]], na.rm = TRUE),
			maximum_raw_score = max(.data[[paste0("diet.", nm, ".sum")]], na.rm = TRUE), .groups = "drop")
}))
group_audit$rule <- "Unweighted rank quartiles among White participants with each score: Q1 low; Q2-Q3 middle; Q4 high; ties may split"
write_xlsx(group_audit, "Reviewer.group_cutpoints.xlsx")

maha_save_state("phe", c("dat"))
maha_write_phenotypes(dat, "eid")
