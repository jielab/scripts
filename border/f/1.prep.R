source(file.path(Sys.getenv("BORDER_FDIR"), "0.common.R"))
follow_end_source <- if (nzchar(Sys.getenv("BORDER_FOLLOW_END", ""))) "--follow-end" else default_follow_end_path
effective_follow_end <- follow_end()
# Freeze the resolved date for this preparation, including its metadata.
Sys.setenv(BORDER_FOLLOW_END = as.character(effective_follow_end))
cat("Follow-up end:", as.character(effective_follow_end), "from", follow_end_source, "\n")
dat <- read_all() %>% prep_covars() %>% make_birth_order()
message("[border] PROGRESS cohort: ",nrow(dat)," participants; eligible birth order: ",
        sum(dat$eligible_border,na.rm=TRUE),"; adopted==1 excluded: ",sum(dat$adopted_num==1,na.rm=TRUE))
input <- file.path(indir, "Rdata", "all.rds")
info <- file.info(input)
reconstructed <- OUTCOMES[!vapply(OUTCOMES, function(y) all(paste0(y, c(".t2e", ".Yt2e")) %in% names(dat)), logical(1))]
outcome_qc <- bind_rows(lapply(OUTCOMES, function(y) {
  source_cols <- if (y %in% reconstructed) paste(outcome_date_columns(dat, y), collapse=";") else "existing t2e/Yt2e"
  dat <<- ensure_outcome(dat, y)
  tibble(outcome=y, source=source_cols, follow_end=if(y %in% reconstructed) as.character(follow_end()) else "upstream all.rds",
         N=sum(!is.na(dat[[paste0(y,".t2e")]])), events=sum(dat[[paste0(y,".Yt2e")]] == 1, na.rm=TRUE))
}))
keep <- unique(c("eid","cnt_older","cnt_brother","cnt_sister","adopted","sibling_n","sibship_size","birth_rank","impossible_rank","eligible_border", ALL_EXPOSURES, covars,
                 unlist(lapply(OUTCOMES, function(y) paste0(y, c(".t2e", ".Yt2e"))))))
border <- dat[,keep,drop=FALSE]
attr(border, "border_metadata") <- list(input=normalizePath(input), input_size=info$size,
  input_mtime=as.numeric(info$mtime), outcomes=OUTCOMES, reconstructed=reconstructed,
  follow_end=as.character(effective_follow_end), follow_end_source=follow_end_source)
saveRDS(border, file.path(rawdir,"border.cohort.rds"))
qc <- bind_rows(
  tibble(metric="all.rds N", value=nrow(dat)),
  tibble(metric="complete three sibling counts", value=sum(complete.cases(dat[,c("cnt_older","cnt_brother","cnt_sister")]))),
  tibble(metric="adopted==1 excluded", value=sum(dat$adopted_num==1, na.rm=TRUE)),
  tibble(metric="adoption status unknown (retained if otherwise eligible)", value=sum(is.na(dat$adopted_num))),
  tibble(metric="impossible cnt_older > siblings", value=sum(dat$impossible_rank, na.rm=TRUE)),
  tibble(metric="eligible birth-order N", value=sum(dat$eligible_border, na.rm=TRUE)),
  bind_rows(lapply(ALL_EXPOSURES, function(x) tibble(metric=paste0(x," nonmissing"), value=sum(!is.na(dat[[x]])))))
)
write_csv2(qc,"cohort_qc.csv")
write_csv2(outcome_qc, "outcome_qc.csv")
writeLines(capture.output(sessionInfo()), file.path(rawdir,"sessionInfo.txt"))
cat("Saved ", file.path(rawdir,"border.cohort.rds"), "\n", sep="")
print(qc); print(outcome_qc)
