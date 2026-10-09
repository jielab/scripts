# Run from the repository root with Rscript --vanilla tests/test_shiny_evidence.R.
suppressPackageStartupMessages({library(shiny); library(DT); library(data.table)})
source("shiny/phyml.R")
source("shiny/evidence.R")
`%||%` <- function(x, y) if(is.null(x) || !length(x) || is.na(x) || !nzchar(x))y else x

# A tree-supported row stays visible even with poor/missing sequence QC.
report <- data.frame(
	locus_key = c("unsupported", "supported-poor", "supported-adequate", "supported-unknown"),
	locus_id = c("l0", "l1", "l2", "l3"), lineage = c("Neanderthal", "Neanderthal", "Denisovan", "Neanderthal"),
	call = c("tree_not_supported", rep("tree_supported", 3)), sequence_qc = c("adequate", "low", "adequate", NA),
	dataset_id = "1kg", genome_build = "GRCh37", chr = c("1", "3", "11", "X"),
	core_start = c(0, 100, 200, 300), core_end = c(90, 190, 290, 390),
	n_candidate_haplotypes = 1, n_candidate_copies = 2, ibdmix_status = "not_run")
stopifnot(identical(gu_phyml_supported(report)$locus_id, c("l1", "l2", "l3")),
	nrow(gu_phyml_supported(data.frame())) == 0,
	nrow(gu_phyml_supported(transform(report, call = NA_character_))) == 0)
dual <- report; dual$nonrisk_tree_pass <- c(1, 0, 0, NA)
stopifnot(identical(gu_phyml_visible(dual)$locus_id, c("l0", "l1", "l2", "l3")),
	identical(gu_phyml_visible(dual, "nonrisk")$locus_id, "l0"),
	identical(gu_phyml_visible(dual, "risk")$locus_id, c("l1", "l2", "l3")),
	nrow(gu_phyml_visible(dual, "all")) == 4)
validation <- data.frame(
	locus_key = c(rep("supported-poor", 5), "supported-adequate", "supported-unknown"),
	lineage = c(rep("Neanderthal", 5), "Denisovan", "Neanderthal"), method = "ibdmix",
	sample_id = c("A", "A", "B", "C", "C", "D", "E"), method_complete = c(1,1,1,1,1,1,0),
	comparison_available = c(1,1,1,0,0,0,0), evidence_eligible = 1, overlap_pass = c(1,1,0,1,1,1,1))
overview <- gu_phyml_overview(report, data.frame(), validation)
stopifnot(identical(overview$ibdmix_risk_support, c("未运行", "1 / 2（部分可评估）", "未评估", "未运行")))

# Filtering must not confuse overview row indices with full-report indices.
root <- tempfile("gu-evidence-test-"); dir.create(root)
write.table(report, file.path(root, "phyml_locus_report.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
for(name in c("phyml_haplotype_report", "phyml_copy_validation"))
	write.table(report[0, ], file.path(root, paste0(name, ".tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
testServer(function(input, output, session) {
	external <- reactiveVal(NULL)
	selection <- gu_phyml_report_server(input, output, session, reactive("1kg"), reactive("GRCh37"),
		root, external_record = external)
}, {
	session$setInputs(report_lineage = "all")
	stopifnot(selection()$locus_id == "l1")
	session$setInputs(report_locus_open = 2L)
	stopifnot(selection()$locus_id == "l2")
	session$setInputs(report_lineage = "Neanderthal", report_locus_select = 2L)
	stopifnot(selection()$locus_id == "l3")
	external("unsupported|Neanderthal"); session$flushReact()
	stopifnot(selection()$locus_id == "l0") # All-loci detail remains accessible.
	session$setInputs(report_lineage = "no-results")
	stopifnot(nrow(selection()) == 0)
})

con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
segments <- data.frame(
	segment_code = c("a","a","a","a","old","x","failed","missing","other-method","ukb","b38","den"),
	dataset_id = c(rep("1kg",9),"ukb","1kg","1kg"), genome_build = c(rep("GRCh37",10),"GRCh38","GRCh37"),
	sample_id = c("001","001","001","002",rep("003",8)), method = c(rep("ibdmix",8),"trace",rep("ibdmix",3)),
	source = c("Altai","Altai","Vindija","Altai","Altai.2013",rep("Altai",6),"Denisova"),
	source_class = c(rep("Neanderthal",11),"Denisovan"), reference_role = "scientific_reference",
	chr = c(rep("1",5),"X","2","3",rep("1",4)), score = seq_len(12) + 4)
DBI::dbWriteTable(con, "segments", segments)
DBI::dbExecute(con, "CREATE VIEW segments_evidence AS SELECT * FROM segments WHERE source NOT IN ('Altai.2013','Denisova.2013')")
catalog <- unique(segments[,c("segment_code","dataset_id","genome_build","source_class","chr")])
catalog$start <- seq_len(nrow(catalog)) * 100000; catalog$end <- catalog$start + 60000; catalog$length_bp <- 60000
DBI::dbWriteTable(con, "segment_catalog", catalog)
runs <- data.frame(dataset_id = c("1kg","1kg","1kg","ukb","1kg"),
	genome_build = c(rep("GRCh37",4),"GRCh38"), chr = c("1","X","2","1","1"),
	method = "ibdmix", status = c("complete","complete","failed","complete","complete"), evidence_eligible = c(1,0,1,1,1))
DBI::dbWriteTable(con, "method_runs", runs)
Q <- function(sql, params) DBI::dbGetQuery(con, sql, params = params)
d <- gu_ibdmix_evidence_loci(Q, "1kg", "GRCh37")
stopifnot(setequal(d$segment_code, c("a","den")), !anyDuplicated(d$segment_code),
	d$n_carriers[d$segment_code == "a"] == 2L, d$n_references[d$segment_code == "a"] == 2L,
	d$n_calls[d$segment_code == "a"] == 4L, d$max_lod[d$segment_code == "a"] == 8,
	gu_ibdmix_evidence_loci(Q, "ukb", "GRCh37")$segment_code == "ukb",
	gu_ibdmix_evidence_loci(Q, "1kg", "GRCh38")$segment_code == "b38",
	nrow(gu_ibdmix_evidence_loci(Q, "unknown", "GRCh37")) == 0)

# Reuse only the same database version, build and dataset.
db_path <- file.path(root, "db-stamp"); writeLines("v1", db_path)
n_queries <- 0L
countQ <- function(sql, params) {n_queries <<- n_queries + 1L; Q(sql, params)}
cache_dir <- file.path(root, "cache")
for(i in 1:2) gu_ibdmix_evidence_cached(countQ, "1kg", "GRCh37", db_path, cache_dir)
stopifnot(n_queries == 1L)
Sys.setFileTime(db_path, Sys.time() + 10)
invisible(gu_ibdmix_evidence_cached(countQ, "1kg", "GRCh37", db_path, cache_dir))
stopifnot(n_queries == 2L)

testServer(function(input, output, session) {
	selected <- reactiveVal(NULL); opened <- reactiveVal(NULL)
	loci <- gu_ibdmix_evidence_server(input, output, session, Q, reactive("1kg"), reactive("GRCh37"),
		reactive("v1"), db_path,
		function(r, method)selected(r), function(r, method)opened(r))
}, {
	session$setInputs(density_lineage = "all", ibdmix_evidence_chr = "ALL", ibdmix_evidence_carriers = 1)
	stopifnot(nrow(loci()) == 2)
	session$setInputs(ibdmix_evidence_open = "a")
	stopifnot(opened()$start == catalog$start[catalog$segment_code == "a"])
	session$setInputs(ibdmix_evidence_carriers = 2)
	stopifnot(identical(loci()$segment_code, "a"))
	session$setInputs(density_lineage = "Denisovan")
	stopifnot(nrow(loci()) == 0)
	session$setInputs(ibdmix_evidence_carriers = 1, ibdmix_evidence_select = "den")
	stopifnot(identical(loci()$segment_code, "den"), selected()$source_class == "Denisovan")
	session$setInputs(ibdmix_evidence_chr = "X")
	stopifnot(nrow(loci()) == 0)
})
DBI::dbDisconnect(con)
unlink(root, recursive = TRUE)
cat("Shiny evidence checks passed.\n")
