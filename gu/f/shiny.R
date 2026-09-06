# GU Shiny entry point; UI, server and review panels are maintained together.
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
app_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1L]]), mustWork = TRUE))
if ("--review" %in% commandArgs(trailingOnly = TRUE)) {
#!/usr/bin/env Rscript
# Standalone Shiny host for the prepared, self-contained evidence review page.
# Uses the existing shiny package only. Does not modify any environment or database.
args <- commandArgs(trailingOnly = TRUE)
value <- function(flag, default = NULL) {
  ix <- match(flag, args)
  if (is.na(ix)) return(default)
  if (ix == length(args)) stop("Missing value for ", flag)
  args[[ix + 1L]]
}
data_dir <- value("--data")
if (is.null(data_dir)) stop("Usage: Rscript run_shiny.R --data /mnt/d/analysis/gu/review [--port 3839] [--host 127.0.0.1]")
port <- suppressWarnings(as.integer(value("--port", "3839")))
if (is.na(port) || port < 1L || port > 65535L) stop("Invalid port")
host <- value("--host", "127.0.0.1")
if (!requireNamespace("shiny", quietly = TRUE)) stop("The existing R session does not have shiny. Activate your working GU R environment; this script will not install packages.")
data_dir <- normalizePath(data_dir, mustWork = TRUE)
if (!file.exists(file.path(data_dir, "review.html"))) stop("review.html missing. Run prepare_review.py first.")
shiny::addResourcePath("gu_review", data_dir)
ui <- shiny::fluidPage(
  shiny::tags$head(shiny::tags$style("html,body{margin:0;padding:0}.container-fluid{padding:0}")),
  shiny::tags$iframe(src = "gu_review/review.html", title = "GU evidence review",
                    style = "display:block;width:100%;height:100vh;border:0")
)
app <- shiny::shinyApp(ui, server = function(input, output, session) {})
shiny::runApp(app, host = host, port = port, launch.browser = FALSE)

  quit(status = 0)
}

# Shared data access
suppressPackageStartupMessages({
  library(shiny); library(bslib); library(plotly); library(DT); library(DBI); library(RSQLite); library(data.table); library(scales); library(ape)
})

`%||%` <- function(x,y) if (is.null(x) || length(x)==0 || is.na(x) || !nzchar(x)) y else x
if (!exists("app_dir", inherits = FALSE)) {
  ofile <- tryCatch(sys.frame(1)$ofile, error=function(e) NULL)
  app_dir <- normalizePath(dirname(ofile %||% getwd()), mustWork = FALSE)
}

default_db <- normalizePath(file.path(Sys.getenv("GU_ANALYSIS_ROOT", "/mnt/d/analysis/gu"), "gu.sqlite"), mustWork = FALSE)
db_path <- Sys.getenv("GU_SQLITE", unset = default_db)
if (!file.exists(db_path)) {
  stop("GU SQLite database not found: ", db_path, "\nRun: ./gu.sh final")
}
db_path <- normalizePath(db_path, mustWork = TRUE)
normalize_root <- dirname(db_path)
.gu_resolve_artifact <- function(path) {
  value <- as.character(path[[1]])
  if (is.na(value) || !nzchar(value)) return(value)
  is_absolute <- grepl("^(/|[A-Za-z]:[/\\\\])", value)
  normalizePath(if (is_absolute) value else file.path(normalize_root, value), mustWork = FALSE)
}

chr_order <- c(as.character(1:22), "X", "Y", "MT")
fmt_bp <- function(x) ifelse(x >= 1e6, paste0(round(x/1e6,2)," Mb"), ifelse(x>=1e3,paste0(round(x/1e3,1)," kb"),paste0(x," bp")))
q <- function(con, sql, params = NULL) {
  if (is.null(params) || length(params) == 0L) DBI::dbGetQuery(con, sql)
  else DBI::dbGetQuery(con, sql, params = params)
}
scalar_q <- function(con, sql, params=NULL) { z <- q(con,sql,params); if(nrow(z)) z[[1]][1] else NA }

reduce_intervals <- function(d) {
  if (!nrow(d)) return(d[, .(start,end)])
  x <- as.data.table(d)[order(start,end), .(start=as.numeric(start),end=as.numeric(end))]
  out <- vector("list", nrow(x)); k <- 0L; s <- x$start[1]; e <- x$end[1]
  if (nrow(x)>1) for(i in 2:nrow(x)) {
    if (x$start[i] <= e) e <- max(e,x$end[i]) else { k<-k+1L; out[[k]]<-c(s,e); s<-x$start[i]; e<-x$end[i] }
  }
  k<-k+1L; out[[k]]<-c(s,e)
  m <- do.call(rbind,out[seq_len(k)])
  data.table(start=as.numeric(m[,1]), end=as.numeric(m[,2]))
}

interval_bp <- function(d) { r <- reduce_intervals(as.data.table(d)); if(!nrow(r)) 0 else sum(r$end-r$start) }
intersection_bp <- function(a,b) {
  a <- reduce_intervals(as.data.table(a)); b <- reduce_intervals(as.data.table(b)); i<-1L;j<-1L;z<-0
  while(i<=nrow(a) && j<=nrow(b)) { z<-z+max(0,min(a$end[i],b$end[j])-max(a$start[i],b$start[j])); if(a$end[i]<b$end[j]) i<-i+1L else j<-j+1L }
  z
}

# GU evidence review add-on (read-only generated HTML; no database mutation).
.gu_review_root <- Sys.getenv("GU_PHYML_REPORT_DIR", file.path(dirname(db_path), "review"))
source(file.path(app_dir,"shiny_phyml_report.R"),local=TRUE)
if (file.exists(file.path(.gu_review_root, "review.html"))) {
  shiny::addResourcePath("gu_evidence_review", .gu_review_root)
}

# GU_DUAL_LEAD_V2_BEGIN

# GU_DUAL_LEAD_V2_END

# Dual-lead panels
# Read-only native Shiny panels. Uses packages already required by the GU viewer.
# No GWAS query, no VCF/caller rerun, no environment or SQLite changes.
gu_dual_lead_server <- function(input, output, session, selected_locus, review_root) {
  read_table <- function(path) {
    if (!file.exists(path)) return(data.frame())
    tryCatch(utils::read.delim(path, sep = "\t", header = TRUE,
              check.names = FALSE, colClasses = "character", na.strings = c("", "NA"),
              quote = "", comment.char = "", fileEncoding = "UTF-8"),
             error = function(e) data.frame())
  }
  # Poll small result tables, never the VCF or the large segment table.
  table_reactive <- function(name) shiny::reactiveFileReader(
    intervalMillis = 3000, session = session,
    filePath = file.path(review_root, paste0(name, ".tsv")), readFunc = read_table)
  lead <- table_reactive("lead_comparison")
  summary <- table_reactive("locus_tag_summary")
  pop <- table_reactive("tag_population_metrics")
  lookup <- table_reactive("gwas_lookup")
  value <- function(row, name, default = "—") {
    if (!name %in% names(row) || !nrow(row)) return(default)
    x <- row[[name]][[1L]]
    if (is.na(x) || !nzchar(as.character(x))) default else as.character(x)
  }
  coordinate_key <- shiny::reactive({
    r <- selected_locus()
    if (is.null(r) || !nrow(r)) return(NA_character_)
    # input_start/end, not the LD-selected range or the fourth-column label.
    st <- value(r, "input_start", value(r, "core_start", NA_character_))
    en <- value(r, "input_end", value(r, "core_end", NA_character_))
    if (is.na(st) || is.na(en)) return(NA_character_)
    ch <- sub("^chr", "", value(r, "chr"), ignore.case = TRUE)
    if (ch == "23") ch <- "X"
    paste0(value(r, "dataset_id"), "|", value(r, "genome_build"), "|", ch, ":", st, "-", en)
  })
  subset_locus <- function(d) {
    k <- coordinate_key()
    if (!nrow(d) || !"locus_key" %in% names(d) || is.na(k)) return(data.frame())
    d[!is.na(d$locus_key) & d$locus_key == k, , drop = FALSE]
  }
  num <- function(x) {
    y <- suppressWarnings(as.numeric(x))
    if (length(y) != 1L || is.na(y)) "—" else sprintf("%.3f", y)
  }
  output$dual_lead_cards <- shiny::renderUI({
    s <- subset_locus(summary()); d <- subset_locus(lead())
    if (!nrow(s)) return(shiny::tags$p(class = "text-muted",
      "Dual-lead results are not available. Run gu/f/prepare_review.py after normalize; existing regional results are unchanged."))
    if (!nrow(d)) return(shiny::tagList(
      shiny::tags$p(shiny::tags$b("Input lead: "), value(s, "input_lead_snp")),
      shiny::tags$p(shiny::tags$b("Tag status: "), value(s, "tag_status"), " · ", value(s, "tag_detail")),
      shiny::tags$p("No best SNP has been invented from summary-only data.")))
    d <- d[d$target_definition == "best_haplotype", , drop = FALSE]
    box <- function(role, title) {
      z <- d[d$role == role, , drop = FALSE]
      shiny::tags$div(class = "col-md-6",
        shiny::tags$div(class = "border rounded p-3 h-100",
          shiny::tags$h5(title),
          shiny::tags$p(shiny::tags$b(value(z, "snp_id", value(z, "variant_key", value(z, "requested_snp"))))),
          shiny::tags$p("Position (1-based): ", value(z, "pos_1based"),
                        " · tag allele: ", value(z, "tag_allele")),
          shiny::tags$p("SNP–haplotype r²: ", num(value(z, "r2")),
                        " · PPV: ", num(value(z, "ppv")),
                        " · sensitivity: ", num(value(z, "sensitivity"))),
          shiny::tags$small(value(z, "status"))))
    }
    shiny::tagList(
      shiny::tags$p(shiny::tags$b(value(s, "locus_key")), " · best haplotype: ",
                    value(s, "best_haplotype_id"), " · ", value(s, "best_haplotype_tier")),
      shiny::tags$div(class = "row g-3", box("input_lead", "Input-file lead SNP"),
                     box("best_tag", "1KG best-matched haplotype · best tag SNP")),
      shiny::tags$p(class = "text-muted mt-2",
        "The input SNP is a comparator, not an inclusion rule. A tag allele is not automatically a GWAS risk allele. A strong in-sample proxy is not proof of introgression or cross-population validity."))
  })
  output$dual_lead_table <- DT::renderDT({
    d <- subset_locus(lead())
    if (!nrow(d)) return(DT::datatable(data.frame(Status = "No dual-lead table yet"), rownames = FALSE, options = list(dom = "t")))
    keep <- intersect(c("target_definition", "role", "snp_id", "variant_key", "tag_allele",
                        "r2", "ppv", "sensitivity", "n_callable_copies", "n_target_callable",
                        "call_rate", "status", "n_equivalent_best_tags"), names(d))
    DT::datatable(d[, keep, drop = FALSE], rownames = FALSE, options = list(dom = "t", scrollX = TRUE))
  })
  output$dual_lead_population <- DT::renderDT({
    d <- subset_locus(pop())
    if (!nrow(d)) return(DT::datatable(data.frame(Status = "Population metrics unavailable"), rownames = FALSE, options = list(dom = "t")))
    d <- d[d$target_definition == "best_haplotype", , drop = FALSE]
    keep <- intersect(c("population", "role", "snp_id", "tag_allele", "r2", "ppv", "sensitivity",
                        "n_callable_copies", "n_target_callable", "metric_status", "selection_population"), names(d))
    DT::datatable(d[, keep, drop = FALSE], rownames = FALSE,
                  options = list(pageLength = 12, scrollX = TRUE))
  })
  output$download_dual_lead_gwas <- shiny::downloadHandler(
    filename = function() paste0("GU_", gsub("[^A-Za-z0-9._-]", "_", coordinate_key()), "_GWAS_lookup.tsv"),
    content = function(file) utils::write.table(subset_locus(lookup()), file, sep = "\t",
                  quote = FALSE, row.names = FALSE, na = ""))
}

# Application UI
ui <- page_navbar(
  title = "GU — Archaic Introgression Browser", id = "nav", fillable = FALSE,
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  header = tagList(
    tags$head(tags$style(HTML("
      .gu-browser-toolbar{display:flex;align-items:center;gap:.5rem;flex-wrap:wrap;margin-bottom:.6rem}
      .gu-browser-toolbar .gu-locus-label{font-weight:700;margin-right:auto}
      .gu-igv-frame{width:100%;height:330px;border:1px solid #ccd3da;border-radius:6px;background:#fff}
      .gu-header-filters{display:flex;gap:1rem;align-items:end;flex-wrap:wrap}
      .gu-method-note{border-left:4px solid #18bc9c;padding:.65rem .9rem;background:#f4fbf9;margin:.4rem 0 .8rem}
      .gu-hap-scroll{overflow-x:auto;max-width:100%;padding-bottom:.5rem}
      .gu-hap-table{border-collapse:collapse;font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;font-size:13px;line-height:1}
      .gu-hap-table th.gu-rowlab{position:sticky;left:0;z-index:3;background:#fff;text-align:right;white-space:nowrap;padding:5px 10px 5px 4px;border-right:1px solid #8d99a6;font-family:var(--bs-body-font-family);font-weight:500}
      .gu-hap-table th.gu-site{height:94px;min-width:20px;max-width:20px;padding:0;vertical-align:bottom;color:#607080;font-size:10px;font-weight:400}
      .gu-hap-table th.gu-site span{display:inline-block;transform:rotate(-65deg);transform-origin:bottom left;white-space:nowrap;margin-left:12px}
      .gu-hap-table td.gu-base{min-width:20px;width:20px;height:23px;text-align:center;vertical-align:middle;padding:0;border-right:1px solid #edf0f2;font-weight:800}
      .gu-hap-table td.gu-base-match{background:#dcf4e4;box-shadow:inset 0 -3px #198754}
      .gu-hap-table tr.gu-arch-last th,.gu-hap-table tr.gu-arch-last td{border-bottom:2px solid #34495e}
      .gu-hap-table tr.gu-control-first th,.gu-hap-table tr.gu-control-first td{border-top:2px solid #c0392b}
      .gu-hap-table tr.gu-control th.gu-rowlab{color:#a93226}
      .gu-hap-legend{display:flex;gap:.8rem;align-items:center;flex-wrap:wrap;margin:.2rem 0 .7rem}
      .gu-hap-legend span{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;font-weight:800}
      table.dataTable tbody tr{cursor:pointer}
    "))),
    div(class = "container-fluid py-2 gu-header-filters",
        selectInput("genome_build", "Genome build", choices = character(0), width = "220px"),
        selectInput("target_dataset", "Target dataset", choices = character(0), width = "220px"))
  ),
  nav_panel("Overview", value = "overview",
    card(
      card_header(tags$b("Compact genomic context — IGV-Web / UCSC")),
      layout_columns(
        selectInput("browser_locus", "PhyML locus", choices = character(0)),
        numericInput("browser_flank", "Context flank (bp)", 250000, min = 0, max = 5000000, step = 50000),
        selectInput("support_population", "Population", "ALL"),
        selectInput("support_denominator", "Curve denominator",
                    choices = c("PhyML-positive"="phyml_positive", "All PhyML-tested"="tested"),
                    selected = "phyml_positive"),
        tags$p(class = "text-muted mt-4",
               "IGV supplies gene context. GU evidence is shown in the aligned prevalence curves below; download the bedGraph to add those tracks to IGV."),
        col_widths = c(3, 2, 2, 2, 3)
      ),
      uiOutput("genome_browser")
    ),
    layout_columns(
      value_box("Samples tested", textOutput("n_samples", inline = TRUE)),
      value_box("Tree-clade carriers", textOutput("n_phyml_positive", inline = TRUE)),
      value_box("Independent overlap", textOutput("n_other_supported", inline = TRUE)),
      value_box("Selected-locus state", textOutput("locus_evidence_state", inline = TRUE)),
      value_box("Methods completed here", textOutput("n_methods", inline = TRUE))
    ),
    card(card_header("Locus evidence overview — click to update browser; double-click for sequence view"),
         DTOutput("overview_loci_table")),
    card(
      card_header("Selected-locus interpretation"),
      tags$div(class = "gu-method-note",
               "PhyML carriers are defined by a bootstrap-supported archaic/modern tree edge. Pairwise sequence identity is descriptive and is not used as the carrier rule. IBDmix, TRACE, and AS3 remain separate evidence streams; not run, unsupported, and exploratory are never converted to zero evidence."),
      uiOutput("locus_interpretation"),
      downloadButton("download_locus_evidence", "Download locus evidence")
    ),
    card(card_header("Independent method support"), DTOutput("locus_support_table")),
    accordion(
      id = "overview_details", open = FALSE,
      accordion_panel("Regional prevalence curves", plotlyOutput("support_trajectory", height = "430px")),
      accordion_panel("Individual support rows",
        layout_columns(
          checkboxInput("support_phyml_only", "Show tree-clade carriers only", TRUE),
          downloadButton("download_locus_tracks", "Download IGV bedGraph"),
          downloadButton("download_locus_support", "Download individual support"),
          col_widths = c(5, 3, 4)
        ),
        tags$p(class = "text-muted", "The interactive table is capped at 2,000 rows; the download contains every matching individual."),
        DTOutput("locus_sample_support"))
    )
  ),
  nav_panel("PhyML evidence", value = "evidence_review",
    gu_phyml_report_ui()
  ),
  nav_panel("PhyML loci", value = "phyml",
    layout_sidebar(
      sidebar(
        selectInput("locus_chr", "Chromosome", c("ALL", chr_order), "ALL"),
        selectInput("locus_status", "Run status", c("ALL", "pass"), "ALL"),
        numericInput("hap_n_match", "Matched haplotypes", 8, min = 1, max = 50),
        numericInput("hap_n_control", "Non-matched controls", 10, min = 0, max = 50),
        numericInput("hap_max_sites", "Maximum SNPs", 150, min = 20, max = 1000, step = 20),
        open = "open"
      ),
      card(card_header("Tree-defined modern–archaic haplotype evidence"), DTOutput("loci_table")),
      card(card_header("Selected locus"), uiOutput("haplotype_title")),
      # GU_DUAL_LEAD_V2_BEGIN
      tags$p("逐谱系与逐单倍型结果、双 SNP 指标及 IBDmix / TRACE 验证，请查看 PhyML evidence。"),
      # GU_DUAL_LEAD_V2_END
      card(
        card_header("Archaic references, candidate-clade modern haplotypes, and controls (A/C/G/T)"),
        tags$p(class = "text-muted",
               "Archaic references are shown first. A green cell background marks a modern base matching that row's best archaic reference; the letter colour identifies A/C/G/T."),
        uiOutput("haplotype_matrix")
      ),
      accordion(
        id = "phyml_details", open = FALSE,
        accordion_panel("Match statistics", DTOutput("haplotype_similarity")),
        accordion_panel("Raw haplotype groups", DTOutput("haplotype_table")),
        accordion_panel("Display provenance / technical details", verbatimTextOutput("haplotype_note"))
      ),
      card(card_header("PhyML candidate-clade tree for selected locus"),
           tags$p(class = "text-muted", "The shaded edge-side is selected directly from the tree: the selected candidate lineage and the run-specific bootstrap, purity, sensitivity and candidate-type criteria recorded in the evidence metadata. Without an explicit ancestral outgroup this is an unrooted ML topology; display rooting does not imply ancestry direction."),
           verbatimTextOutput("tree_summary"), plotOutput("phyml_tree", height = "520px")),
      card(card_header("Samples carrying selected haplotypes"), DTOutput("haplotype_samples"))
    )
  ),
  nav_panel("Segments", value = "segments",
    layout_sidebar(
      sidebar(selectInput("seg_method", "Method", "ALL"), selectInput("seg_chr", "Chromosome", c("ALL", chr_order), "ALL"),
              numericInput("seg_start", "Start (0-based)", 0, min = 0), numericInput("seg_end", "End", 200000000, min = 1),
              checkboxInput("show_reference", "Overlay published callset", TRUE),
              selectInput("reference_population", "Published population", "ALL"),
              sliderInput("seg_limit", "Maximum rows", 100, 10000, 2000, 100), actionButton("seg_go", "Query"), open = "open"),
      card(card_header("Regional segment landscape"), plotlyOutput("segment_plot", height = "470px")),
      card(card_header("Segment calls — click a row to update the Overview browser"), DTOutput("segment_table")),
      card(card_header("Published callset overlay — external_reference, never model input"), DTOutput("published_callset_table")),
      card(card_header("GU-AS3 overlap / reciprocal-overlap by 1KG population"), DTOutput("reference_overlap_table"))
    )
  ),
  nav_panel("Individuals", value = "individuals",
    card(card_header("Per-sample burden"),
         layout_columns(textInput("sample_search", "Sample contains", ""), selectInput("sample_method", "Method", "ALL"),
                        selectInput("burden_type", "Burden definition", c("raw_call", "nonredundant_union", "consensus_catalog"), "nonredundant_union"), col_widths = c(4, 4, 4)),
         DTOutput("sample_table"))
  ),
  nav_panel("Data", value = "data",
    card(card_header("Database"), verbatimTextOutput("db_info")),
    card(card_header("Exports"), tags$p("Published callsets are marked external_reference and are used only for overlays/validation; they are never AS3 model inputs or individual predictions. The SQLite database and normalized TSV files are computation artifacts."), downloadButton("download_segments", "Download filtered segments"))
  )
)

# Server helpers and reactive outputs
server <- local({
.gu_truth <- function(x) {
  toupper(trimws(as.character(x))) %in% c("1", "TRUE", "T", "YES", "Y", "PASS")
}

.gu_browser_urls <- function(target, flank = 250000L) {
  build <- as.character(target$genome_build[[1]])
  ids <- switch(build,
    GRCh37 = list(igv = "hg19", ucsc = "hg19"),
    GRCh38 = list(igv = "hg38", ucsc = "hg38"),
    CHM13 = list(igv = "chm13v2.0", ucsc = "hs1"),
    list(igv = "hg19", ucsc = "hg19")
  )
  chrom <- sub("^chr", "", as.character(target$chr[[1]]), ignore.case = TRUE)
  start <- max(1L, suppressWarnings(as.integer(target$start[[1]])) - as.integer(flank))
  end <- suppressWarnings(as.integer(target$end[[1]])) + as.integer(flank)
  locus <- paste0("chr", chrom, ":", start, "-", end)
  list(
    label = as.character(target$label[[1]] %||% "Selected region"),
    build = build, locus = locus,
    igv = paste0("https://igv.org/app/?genome=", ids$igv, "&locus=", utils::URLencode(locus, reserved = TRUE)),
    ucsc = paste0("https://genome.ucsc.edu/cgi-bin/hgTracks?db=", ids$ucsc, "&position=", utils::URLencode(locus, reserved = TRUE))
  )
}

.gu_cap_sites <- function(n, max_sites) {
  if (n <= max_sites) return(seq_len(n))
  unique(pmax(1L, pmin(n, as.integer(round(seq(1, n, length.out = max_sites))))))
}

.gu_phyml_paths <- function(locus) {
  raw <- .gu_resolve_artifact(locus$raw_file)
  if (is.na(raw) || !nzchar(raw)) stop("The normalized locus has no raw_file provenance.")
  final_dir <- dirname(raw)
  run_dir <- if (basename(final_dir) == "final") dirname(final_dir) else final_dir
  locus_id <- as.character(locus$locus_id[[1]])
  flat_locus_dir <- file.path(run_dir, "loci")
  nested_locus_dir <- file.path(flat_locus_dir, locus_id)
  flat_artifacts <- file.path(flat_locus_dir, c("sites.tsv", "archaic.tsv"))
  locus_dir <- if (any(file.exists(flat_artifacts))) flat_locus_dir else nested_locus_dir
  list(
    final = final_dir,
    run = run_dir,
    locus = locus_dir,
    sites = file.path(locus_dir, "sites.tsv"),
    archaic = file.path(locus_dir, "archaic.tsv"),
    haplotypes = file.path(final_dir, "haplotypes.tsv"),
    samples = file.path(final_dir, "haplotype_samples.tsv")
  )
}

.gu_read_phyml_view <- function(locus, n_match = 8L, n_control = 10L, max_sites = 150L) {
  paths <- .gu_phyml_paths(locus)
  missing <- c(paths$sites, paths$archaic, paths$haplotypes)
  missing <- missing[!file.exists(missing) | file.info(missing)$size <= 0]
  if (length(missing)) stop("Missing PhyML sequence artifact(s): ", paste(missing, collapse = "; "))

  sites <- data.table::fread(paths$sites, showProgress = FALSE)
  archaic <- data.table::fread(paths$archaic, showProgress = FALSE)
  hap <- data.table::fread(paths$haplotypes, showProgress = FALSE)
  target_locus_id <- as.character(locus$locus_id[[1]])
  if ("locus_id" %in% names(hap)) hap <- hap[as.character(hap$locus_id) == target_locus_id]
  required_sites <- c("chr", "pos", "ref", "alt")
  if (!all(required_sites %in% names(sites)) || !all(c("archaic", "seq") %in% names(archaic)) ||
      !all(c("hap_id", "seq") %in% names(hap)))
    stop("Malformed PhyML sites/archaic/haplotype artifact for ", target_locus_id)
  if (!nrow(sites) || !nrow(archaic) || !nrow(hap)) stop("No displayable sequence rows for ", target_locus_id)

  sequence_lengths <- c(nchar(as.character(archaic$seq)), nchar(as.character(hap$seq)))
  sequence_lengths <- sequence_lengths[is.finite(sequence_lengths)]
  n_sites <- min(c(nrow(sites), sequence_lengths))
  if (!is.finite(n_sites) || n_sites < 1L) stop("Sequence/site lengths do not overlap for ", target_locus_id)
  sites <- sites[seq_len(n_sites)]
  idx <- .gu_cap_sites(n_sites, max(20L, as.integer(max_sites)))

  if (!"direct_match_pass" %in% names(hap)) hap[, direct_match_pass := FALSE]
  if (!"prop_match" %in% names(hap)) hap[, prop_match := NA_real_]
  if (!"n" %in% names(hap)) hap[, n := NA_integer_]
  if (!"best_archaic" %in% names(hap)) hap[, best_archaic := NA_character_]
  if (!"best_lineage" %in% names(hap)) hap[, best_lineage := NA_character_]
  if (!"n_compared" %in% names(hap)) hap[, n_compared := NA_integer_]
  if (!"n_match" %in% names(hap)) hap[, n_match := NA_integer_]
  candidate_tips <- unique(strsplit(as.character(locus$candidate_clade_tips[[1]] %||% ""), ",", fixed=TRUE)[[1]])
  candidate_ok <- .gu_truth(locus$candidate_clade_pass[[1]])
  hap[, `:=`(.pass = candidate_ok & as.character(hap_id) %in% candidate_tips, .prop = suppressWarnings(as.numeric(prop_match)),
             .n = suppressWarnings(as.numeric(n)))]
  data.table::setorder(hap, -.pass, -.prop, -.n, hap_id, na.last = TRUE)
  matched <- head(hap[.pass == TRUE], max(1L, as.integer(n_match)))
  remaining <- hap[!hap_id %in% matched$hap_id]
  controls <- head(remaining[.pass == FALSE], max(0L, as.integer(n_control)))

  split_seq <- function(x) substring(as.character(x), idx, idx)
  arch_rows <- lapply(seq_len(nrow(archaic)), function(i) list(
    kind = "archaic", id = as.character(archaic$archaic[[i]]),
    label = paste0("[ARCH] ", archaic$archaic[[i]], if ("lineage" %in% names(archaic)) paste0(" · ", archaic$lineage[[i]]) else ""),
    bases = split_seq(archaic$seq[[i]]), best_archaic = as.character(archaic$archaic[[i]]),
    tooltip = paste0("Archaic reference: ", archaic$archaic[[i]])
  ))
  modern_rows <- function(d, kind) lapply(seq_len(nrow(d)), function(i) {
    pct <- suppressWarnings(as.numeric(d$prop_match[[i]]))
    label <- paste0(d$hap_id[[i]], " · n=", d$n[[i]], " · ", d$best_archaic[[i]],
                    if (is.finite(pct)) sprintf(" · %.1f%%", 100 * pct) else "")
    list(kind = kind, id = as.character(d$hap_id[[i]]), label = label,
         bases = split_seq(d$seq[[i]]), best_archaic = as.character(d$best_archaic[[i]]),
         tooltip = paste0("Haplotype ", d$hap_id[[i]], "\nCopies: ", d$n[[i]],
                          "\nBest archaic: ", d$best_archaic[[i]],
                          "\nCompared/matched: ", d$n_compared[[i]], "/", d$n_match[[i]],
                          if (is.finite(pct)) sprintf("\nProportion matched: %.2f%%", 100 * pct) else ""))
  })
  rows <- c(arch_rows, modern_rows(matched, "matched"), modern_rows(controls, "control"))
  list(
    locus = locus, paths = paths, sites = sites, idx = idx, rows = rows,
    archaic = archaic, matched = matched, controls = controls,
    n_sites = n_sites, n_display = length(idx)
  )
}

.gu_descendant_tips <- function(tree, node) {
  tips <- integer(); frontier <- as.integer(node)
  while (length(frontier)) {
    children <- tree$edge[tree$edge[, 1] %in% frontier, 2]
    tips <- c(tips, children[children <= ape::Ntip(tree)])
    frontier <- unique(children[children > ape::Ntip(tree)])
  }
  unique(tips)
}

.gu_candidate_edge_for_tree <- function(tree, targets) {
  targets <- intersect(unique(as.character(targets)), tree$tip.label)
  n_tip <- ape::Ntip(tree)
  if (length(targets) < 2L || n_tip < 3L) return(NULL)
  internal_children <- unique(tree$edge[tree$edge[, 2] > n_tip, 2])
  best <- NULL
  key_is_better <- function(left, right) {
    if (is.null(right)) return(TRUE)
    for (i in seq_along(left)) {
      if (left[[i]] < right[[i]]) return(TRUE)
      if (left[[i]] > right[[i]]) return(FALSE)
    }
    FALSE
  }
  for (node in internal_children) {
    down <- .gu_descendant_tips(tree, node)
    support <- if (length(tree$node.label)) suppressWarnings(as.numeric(tree$node.label[[node - n_tip]])) else NA_real_
    support_rank <- if (is.finite(support)) support else -1
    sides <- list(descendant=down, complement=setdiff(seq_len(n_tip), down))
    for (side_name in names(sides)) {
      side <- sides[[side_name]]; side_labels <- tree$tip.label[side]
      if (!length(side) || length(side) >= n_tip || !all(targets %in% side_labels)) next
      key <- c(length(side), -support_rank, if (identical(side_name,"complement")) 0 else 1, node)
      if (key_is_better(key, if (is.null(best)) NULL else best$key)) {
        best <- list(key=key, tips=side_labels, support=support, side=side_name)
      }
    }
  }
  best
}

.gu_redraw_phylogram_edges <- function(tree, plot_state, color, width) {
  edge <- tree$edge
  for (i in seq_len(nrow(edge))) {
    parent <- edge[i, 1]; child <- edge[i, 2]
    segments(plot_state$xx[parent], plot_state$yy[child],
             plot_state$xx[child], plot_state$yy[child], col=color, lwd=width)
  }
  for (parent in unique(edge[, 1])) {
    children <- edge[edge[, 1] == parent, 2]
    segments(plot_state$xx[parent], min(plot_state$yy[children]),
             plot_state$xx[parent], max(plot_state$yy[children]), col=color, lwd=width)
  }
}

function(input, output, session) {
  disconnect <- function(con) {
    if (!is.null(con) && DBI::dbIsValid(con)) try(DBI::dbDisconnect(con), silent = TRUE)
  }
  initial_con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  con_rv <- reactiveVal(initial_con)
  required_schema <- list(
    loci=c("best_archaic","best_lineage","n_compared","n_match","prop_match","tree_has_ancestral_outgroup","candidate_clade_pass","candidate_clade_tips"),
    locus_evidence=c("phyml_qc","evidence_state","region_coverage_fraction"),
    locus_method_support=c("method_available","evidence_eligible","availability_note")
  )
  schema_errors <- character()
  for (table in names(required_schema)) {
    fields <- if (DBI::dbExistsTable(initial_con,table)) DBI::dbListFields(initial_con,table) else character()
    missing <- setdiff(required_schema[[table]],fields)
    if (length(missing)) schema_errors <- c(schema_errors,paste0(table,".",missing))
  }
  if (length(schema_errors)) {
    disconnect(initial_con)
    stop("GU SQLite schema is older than this Shiny app; rerun ./gu.sh final (missing ",paste(schema_errors,collapse=", "),").")
  }
  session$onSessionEnded(function() disconnect(isolate(con_rv())))

  # normalize_results.py atomically replaces gu.sqlite. Reconnect each Shiny
  # session when that inode changes so a running app sees the new database.
  db_stamp <- reactivePoll(
    1000, session,
    checkFunc = function() {
      info <- file.info(db_path)
      paste(as.numeric(info$mtime), as.numeric(info$ctime), info$size, sep = ":")
    },
    valueFunc = function() {
      info <- file.info(db_path)
      paste(as.numeric(info$mtime), as.numeric(info$ctime), info$size, sep = ":")
    }
  )
  observeEvent(db_stamp(), ignoreInit = TRUE, {
    new_con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
    old_con <- con_rv()
    con_rv(new_con)
    disconnect(old_con)
  })
  Q <- function(sql, params = NULL) q(con_rv(), sql, params)
  SQ <- function(sql, params = NULL) scalar_q(con_rv(), sql, params)
  selected_or <- function(value, choices, fallback = "ALL") {
    if (!is.null(value) && length(value) && value[[1]] %in% choices) value[[1]] else fallback
  }
  has_summary <- function() DBI::dbExistsTable(con_rv(), "database_summary")
  summary_value <- function(metric, fallback_sql) {
    if (has_summary()) {
      x <- SQ("SELECT value FROM database_summary WHERE dataset_id=? AND genome_build=? AND metric=?", list(current_dataset(), current_build(), metric))
      return(if (is.na(x)) 0 else x)
    }
    SQ(fallback_sql, list(current_dataset(), current_build()))
  }

  builds <- reactive({
    db_stamp()
    if (has_summary()) {
      d <- Q("SELECT genome_build,MAX(CASE WHEN metric IN ('segment_calls','phyml_loci') AND value>0 THEN 1 ELSE 0 END) AS has_results FROM database_summary WHERE dataset_id!='reference' GROUP BY genome_build ORDER BY has_results DESC,CASE genome_build WHEN 'GRCh37' THEN 1 WHEN 'GRCh38' THEN 2 WHEN 'CHM13' THEN 3 ELSE 9 END,genome_build")
    } else {
      candidates <- c("GRCh37", "GRCh38", "CHM13")
      d <- do.call(rbind, lapply(candidates, function(b) {
        z <- Q("SELECT EXISTS(SELECT 1 FROM segments WHERE genome_build=? LIMIT 1) AS has_segments,EXISTS(SELECT 1 FROM loci WHERE genome_build=? LIMIT 1) AS has_loci,EXISTS(SELECT 1 FROM reference_callsets WHERE genome_build=? LIMIT 1) AS has_reference", list(b, b, b))
        if (sum(z[1, ]) > 0) data.frame(genome_build = b, has_results = as.integer(z$has_segments || z$has_loci)) else NULL
      }))
      if (is.null(d)) d <- data.frame(genome_build = character(), has_results = integer())
      d <- d[order(-d$has_results, match(d$genome_build, candidates)), , drop = FALSE]
    }
    if (nrow(d)) unique(d[[1]]) else character()
  })
  observe({
    choices <- builds()
    selected <- isolate(input$genome_build)
    if (!length(choices)) choices <- "GRCh37"
    if (is.null(selected) || !selected %in% choices) selected <- choices[[1]]
    updateSelectInput(session, "genome_build", choices = choices, selected = selected)
  })
  current_build <- reactive({ req(input$genome_build); input$genome_build })
  datasets <- reactive({
    db_stamp(); req(current_build())
    d <- Q("SELECT dataset_id FROM segments WHERE genome_build=? UNION SELECT dataset_id FROM loci WHERE genome_build=? ORDER BY dataset_id", list(current_build(), current_build()))
    if (nrow(d)) d[[1]] else character()
  })
  observe({
    choices <- datasets(); selected <- isolate(input$target_dataset)
    if (!length(choices)) choices <- "1kg"
    if (is.null(selected) || !selected %in% choices) selected <- choices[[1]]
    updateSelectInput(session, "target_dataset", choices = choices, selected = selected)
  })
  current_dataset <- reactive({ req(input$target_dataset); input$target_dataset })

  active_locus_rowid <- reactiveVal(NA_integer_)
  browser_override <- reactiveVal(NULL)
  all_loci_data <- reactive({
    db_stamp(); req(current_build(), current_dataset())
    Q("SELECT l.rowid,l.locus_id,l.chr,l.start AS selected_start,l.end AS selected_end,l.input_start,l.input_end,l.analysis_start,l.analysis_end,l.flank_bp,l.anchor_pos,l.selection_method,l.ld_r2_threshold,l.n_search_sites,l.n_ld_sites,l.n_ancestral_sites,l.n_sites,l.best_archaic,l.best_lineage,l.n_compared,l.n_match,l.prop_match,l.source,l.source_class,l.status,l.direct_match_pass,l.n_carriers,l.tree_status,l.n_bootstrap_nodes,l.bootstrap_min,l.bootstrap_median,l.bootstrap_max,l.candidate_lineage,l.tree_has_ancestral_outgroup,l.candidate_clade_pass,l.candidate_clade_rule,l.candidate_clade_bootstrap,l.candidate_clade_node,l.candidate_clade_side,l.candidate_clade_n_tips,l.candidate_clade_modern_tips,l.candidate_clade_archaic_tips,l.candidate_clade_specificity,l.candidate_clade_tips,l.tree_newick,l.raw_file,l.genome_build,l.dataset_id,e.sequence_information,e.region_coverage_fraction,e.phyml_qc,e.evidence_state,e.evidence_summary FROM loci l LEFT JOIN locus_evidence e ON e.dataset_id=l.dataset_id AND e.genome_build=l.genome_build AND e.locus_id=l.locus_id WHERE l.dataset_id=? AND l.genome_build=? AND l.method='phyml' ORDER BY CASE WHEN l.chr GLOB '[0-9]*' THEN CAST(l.chr AS INT) ELSE 99 END,l.start,l.locus_id", list(current_dataset(), current_build()))
  })
  observe({
    d <- all_loci_data()
    choices <- if (nrow(d)) setNames(as.character(d$rowid), paste0(d$locus_id, " | chr", d$chr, ":", format(d$selected_start, big.mark = ","), "-", format(d$selected_end, big.mark = ","))) else character()
    active <- isolate(active_locus_rowid())
    if (!nrow(d)) {
      active_locus_rowid(NA_integer_)
      updateSelectInput(session, "browser_locus", choices = character())
    } else {
      if (!is.finite(active) || !active %in% d$rowid) active <- d$rowid[[1]]
      active_locus_rowid(as.integer(active))
      updateSelectInput(session, "browser_locus", choices = choices, selected = as.character(active))
    }
  })
  observeEvent(input$browser_locus, {
    z <- suppressWarnings(as.integer(input$browser_locus))
    if (length(z) && is.finite(z)) {
      active_locus_rowid(z)
      browser_override(NULL)
    }
  }, ignoreInit = TRUE)
  selected_locus <- reactive({
    d <- all_loci_data(); if (!nrow(d)) return(NULL)
    i <- match(active_locus_rowid(), d$rowid)
    if (is.na(i)) i <- 1L
    d[i, , drop = FALSE]
  })

  methods <- reactive({
    b <- current_build(); target <- current_dataset()
    d <- Q("SELECT method FROM segments WHERE dataset_id=? AND genome_build=? UNION SELECT method FROM loci WHERE dataset_id=? AND genome_build=? ORDER BY method", list(target, b, target, b))
    if (nrow(d)) d[[1]] else character()
  })
  segment_methods <- reactive({
    b <- current_build(); target <- current_dataset()
    d <- Q("SELECT DISTINCT method FROM segments WHERE dataset_id=? AND genome_build=? ORDER BY method", list(target, b))
    if (nrow(d)) d[[1]] else character()
  })
  observe({
    choices <- c("ALL", segment_methods())
    updateSelectInput(session, "seg_method", choices = choices, selected = selected_or(input$seg_method, choices))
    burden_methods <- if (identical(input$burden_type, "consensus_catalog")) c("ALL","consensus") else choices
    updateSelectInput(session, "sample_method", choices = burden_methods, selected = selected_or(input$sample_method, burden_methods))
  })
  observe({
    b <- current_build()
    d <- Q("SELECT DISTINCT population FROM reference_callsets WHERE genome_build=? ORDER BY population", list(b))
    choices <- c("ALL", if (nrow(d)) d[[1]] else character())
    updateSelectInput(session, "reference_population", choices=choices, selected=selected_or(input$reference_population,choices))
  })
  observe({
    b <- current_build()
    d <- Q("SELECT DISTINCT status FROM loci WHERE dataset_id=? AND genome_build=? AND method='phyml' AND status IS NOT NULL ORDER BY status", list(current_dataset(), b))
    choices <- c("ALL", if (nrow(d)) d[[1]] else character())
    updateSelectInput(session, "locus_status", choices = choices, selected = selected_or(input$locus_status, choices))
  })

  output$n_samples <- renderText({
    r <- selected_locus(); if (is.null(r)) return("0")
    format(SQ("SELECT COUNT(*) FROM locus_sample_support WHERE dataset_id=? AND genome_build=? AND locus_id=?", list(current_dataset(), current_build(), r$locus_id[[1]])), big.mark = ",")
  })
  output$n_segments <- renderText({
    format(summary_value("segment_calls", "SELECT COUNT(*) FROM segments WHERE dataset_id=? AND genome_build=?"), big.mark = ",")
  })
  output$n_phyml_positive <- renderText({
    r <- selected_locus(); if (is.null(r)) return("0")
    format(SQ("SELECT COUNT(*) FROM locus_sample_support WHERE dataset_id=? AND genome_build=? AND locus_id=? AND phyml=1", list(current_dataset(), current_build(), r$locus_id[[1]])), big.mark = ",")
  })
  output$n_other_supported <- renderText({
    r <- selected_locus(); if (is.null(r)) return("0")
    format(SQ("SELECT COUNT(*) FROM locus_sample_support WHERE dataset_id=? AND genome_build=? AND locus_id=? AND phyml=1 AND (ibdmix=1 OR trace=1 OR as3=1)", list(current_dataset(), current_build(), r$locus_id[[1]])), big.mark = ",")
  })
  output$n_loci <- renderText({
    format(summary_value("phyml_loci", "SELECT COUNT(*) FROM loci WHERE dataset_id=? AND genome_build=? AND method='phyml'"), big.mark = ",")
  })
  output$n_methods <- renderText({
    r <- selected_locus(); if (is.null(r)) return("0")
    SQ("SELECT COUNT(*) FROM locus_method_support WHERE dataset_id=? AND genome_build=? AND locus_id=? AND population='ALL' AND method_available=1", list(current_dataset(), current_build(), r$locus_id[[1]]))
  })
  selected_evidence <- reactive({
    r <- selected_locus(); if (is.null(r)) return(data.frame())
    Q("SELECT * FROM locus_evidence WHERE dataset_id=? AND genome_build=? AND locus_id=?",
      list(current_dataset(), current_build(), r$locus_id[[1]]))
  })
  output$locus_evidence_score <- renderText({
    d <- selected_evidence(); if (!nrow(d) || is.na(d$evidence_score[[1]])) return("NA")
    sprintf("%.1f / 100", d$evidence_score[[1]])
  })
  output$locus_evidence_tier <- renderText({
    d <- selected_evidence(); if (!nrow(d)) return("not scored")
    gsub("_", " ", d$evidence_tier[[1]])
  })
  output$locus_evidence_state <- renderText({
    d <- selected_evidence(); if (!nrow(d)) return("not evaluated")
    gsub("_", " ", d$evidence_state[[1]])
  })
  output$locus_interpretation <- renderUI({
    d <- selected_evidence(); if (!nrow(d)) return(tags$div(class="alert alert-warning", "No evidence row; rerun ./gu.sh final"))
    coverage <- suppressWarnings(as.numeric(d$region_coverage_fraction[[1]]))
    tags$div(
      tags$h5(gsub("_", " ", d$evidence_state[[1]])),
      tags$p(d$evidence_summary[[1]]),
      tags$dl(class="row mb-0",
        tags$dt(class="col-sm-3", "PhyML QC"), tags$dd(class="col-sm-9", gsub("_", " ", d$phyml_qc[[1]])),
        tags$dt(class="col-sm-3", "Selected / input region"), tags$dd(class="col-sm-9", if (is.finite(coverage)) percent(coverage, accuracy=.1) else "NA"),
        tags$dt(class="col-sm-3", "Callable sites"), tags$dd(class="col-sm-9", paste0(d$n_match[[1]], " / ", d$n_compared[[1]])),
        tags$dt(class="col-sm-3", "Composite score"), tags$dd(class="col-sm-9", "Not calculated; native method statistics are not commensurate probabilities"))
    )
  })
  output$locus_evidence_table <- renderDT({
    d <- selected_evidence()
    if (!nrow(d)) return(datatable(data.frame(Message="No evidence row; rerun ./gu.sh final"),rownames=FALSE,options=list(dom="t")))
    show <- d[, intersect(c("locus_id","chr","source_class","best_archaic","best_lineage","n_tested","n_phyml_carriers","carrier_fraction","n_compared","n_match","prop_match","sequence_information","candidate_clade_bootstrap","candidate_clade_specificity","independent_methods_available","independent_methods_supporting","independent_method_names","tract_support_mean","ibdmix_max_lod","trace_max_posterior","as3_max_score","evidence_score","evidence_completeness","evidence_tier","probability_calibrated"),names(d)),drop=FALSE]
    for (column in intersect(c("carrier_fraction","prop_match","candidate_clade_specificity","tract_support_mean","evidence_completeness"),names(show))) {
      show[[column]] <- ifelse(is.na(show[[column]]),NA,percent(show[[column]],accuracy=.1))
    }
    datatable(show,rownames=FALSE,options=list(dom="t",scrollX=TRUE))
  })
  output$download_locus_evidence <- downloadHandler(
    filename=function(){r<-selected_locus();paste0(current_dataset(),"-",current_build(),"-",r$locus_id[[1]],"-evidence.tsv")},
    content=function(file)data.table::fwrite(selected_evidence(),file,sep="\t")
  )
  output$n_reference_calls <- renderText({
    format(SQ("SELECT COUNT(*) FROM reference_callsets WHERE genome_build=? AND reference_role='external_reference'", list(current_build())), big.mark=",")
  })
  output$overview_loci_table <- renderDT({
    d <- all_loci_data()
    if (!nrow(d)) return(datatable(data.frame(Message="No PhyML loci"), rownames=FALSE, options=list(dom="t")))
    roles <- c(rs35044562="positive control", rs10735079="positive control", rs8176719="negative control", rs7412="negative control", rs143054933="focus")
    support <- Q("SELECT locus_id,method,method_available,evidence_eligible,availability_note,n_tested,n_carriers,n_phyml_carriers,n_phyml_supported,phyml_support_fraction FROM locus_method_support WHERE dataset_id=? AND genome_build=? AND population='ALL' AND method IN ('ibdmix','trace','as3')", list(current_dataset(),current_build()))
    method_label <- function(lid, method) {
      z <- support[support$locus_id==lid & support$method==method,,drop=FALSE]
      if (!nrow(z) || z$method_available[[1]]!=1) {
        if (method=="as3" && d$chr[d$locus_id==lid][[1]]=="X") return("unsupported")
        return("not run")
      }
      prefix <- if (z$evidence_eligible[[1]]==1) "" else "exploratory · "
      if (z$n_phyml_carriers[[1]]>0) {
        rate <- suppressWarnings(as.numeric(z$phyml_support_fraction[[1]]))
        return(paste0(prefix,z$n_phyml_supported[[1]],"/",z$n_phyml_carriers[[1]],if(is.finite(rate)) paste0(" (",percent(rate,accuracy=.1),")") else ""))
      }
      paste0(prefix,z$n_carriers[[1]]," calls / ",z$n_tested[[1]]," tested")
    }
    show <- data.frame(
      role=unname(ifelse(d$locus_id %in% names(roles),roles[d$locus_id],"study locus")),
      locus=d$locus_id,
      region=paste0("chr",d$chr,":",format(d$selected_start,big.mark=","),"–",format(d$selected_end,big.mark=",")),
      phyml_qc=gsub("_"," ",d$phyml_qc),
      tree=paste0(ifelse(d$candidate_clade_pass==1,"PASS","no pass")," · ",d$candidate_lineage," · BS ",ifelse(is.na(d$candidate_clade_bootstrap),"NA",d$candidate_clade_bootstrap)),
      tree_carriers=paste0(d$n_carriers," / ",vapply(d$locus_id,function(x) SQ("SELECT COUNT(*) FROM locus_sample_support WHERE dataset_id=? AND genome_build=? AND locus_id=?",list(current_dataset(),current_build(),x)),numeric(1))),
      IBDmix=vapply(d$locus_id,method_label,character(1),method="ibdmix"),
      TRACE=vapply(d$locus_id,method_label,character(1),method="trace"),
      AS3=vapply(d$locus_id,method_label,character(1),method="as3"),
      conclusion=gsub("_"," ",d$evidence_state),check.names=FALSE
    )
    selected <- match(active_locus_rowid(), d$rowid)
    datatable(show, selection = list(mode = "single", selected = selected, target = "row"), rownames = FALSE,
              options = list(paging = FALSE, scrollX = TRUE, dom = "t"),
              callback = DT::JS("table.on('dblclick','tbody tr',function(){var i=table.row(this).index();Shiny.setInputValue('overview_locus_dblclick',i+1,{priority:'event'});});"))
  })
  observeEvent(input$overview_loci_table_rows_selected, {
    d <- all_loci_data(); i <- input$overview_loci_table_rows_selected
    if (length(i) && i >= 1L && i <= nrow(d)) {
      active_locus_rowid(as.integer(d$rowid[[i]])); browser_override(NULL)
      updateSelectInput(session, "browser_locus", selected = as.character(d$rowid[[i]]))
    }
  }, ignoreInit = TRUE)
  observeEvent(input$overview_locus_dblclick, {
    d <- all_loci_data(); i <- suppressWarnings(as.integer(input$overview_locus_dblclick))
    if (is.finite(i) && i >= 1L && i <= nrow(d)) {
      active_locus_rowid(as.integer(d$rowid[[i]])); browser_override(NULL)
      bslib::nav_select("nav", "phyml", session = session)
    }
  }, ignoreInit = TRUE)

  observe({
    r <- selected_locus()
    choices <- "ALL"
    if (!is.null(r)) {
      d <- Q("SELECT DISTINCT population FROM locus_method_support WHERE dataset_id=? AND genome_build=? AND locus_id=? ORDER BY CASE population WHEN 'ALL' THEN 0 ELSE 1 END,population",
             list(current_dataset(), current_build(), r$locus_id[[1]]))
      if (nrow(d)) choices <- d[[1]]
    }
    updateSelectInput(session, "support_population", choices = choices,
                      selected = selected_or(input$support_population, choices))
  })
  current_population <- reactive(input$support_population %||% "ALL")
  locus_support_data <- reactive({
    r <- selected_locus(); if (is.null(r)) return(data.frame())
    Q("SELECT method,source_class,method_available,evidence_eligible,availability_note,n_tested,n_carriers,carrier_fraction,ci_low,ci_high,n_phyml_carriers,n_phyml_supported,phyml_support_fraction FROM locus_method_support WHERE dataset_id=? AND genome_build=? AND locus_id=? AND population=? ORDER BY CASE method WHEN 'phyml' THEN 1 WHEN 'ibdmix' THEN 2 WHEN 'trace' THEN 3 WHEN 'as3' THEN 4 ELSE 9 END",
      list(current_dataset(), current_build(), r$locus_id[[1]], current_population()))
  })
  locus_trajectory_data <- reactive({
    r <- selected_locus(); if (is.null(r)) return(data.frame())
    mode <- input$support_denominator %||% "phyml_positive"
    Q("SELECT method,source_class,bin_start,bin_end,bin_mid,n_denominator,n_carriers,prevalence,ci_low,ci_high FROM locus_trajectory WHERE dataset_id=? AND genome_build=? AND locus_id=? AND population=? AND denominator_mode=? ORDER BY method,bin_start",
      list(current_dataset(), current_build(), r$locus_id[[1]], current_population(), mode))
  })
  output$locus_support_table <- renderDT({
    d <- locus_support_data()
    if (!nrow(d)) return(datatable(data.frame(Message = "No normalized support data; rerun ./gu.sh final"), rownames = FALSE, options = list(dom = "t")))
    d$availability <- ifelse(d$method_available == 1 & d$evidence_eligible == 1, "completed",
                             ifelse(d$method_available == 1, "completed · exploratory", "not run / unsupported"))
    # Not run is missing evidence, not a measured zero prevalence.
    inactive <- is.na(d$method_available) | d$method_available != 1
    missing_cols <- intersect(c("n_carriers", "carrier_fraction", "ci_low", "ci_high",
                                "n_phyml_supported", "phyml_support_fraction"), names(d))
    d[inactive, missing_cols] <- NA
    d$carrier_prevalence <- ifelse(is.na(d$carrier_fraction), NA, percent(d$carrier_fraction, accuracy = 0.1))
    d$prevalence_ci95 <- ifelse(is.na(d$ci_low), NA, paste0(percent(d$ci_low, accuracy = 0.1), " – ", percent(d$ci_high, accuracy = 0.1)))
    d$phyml_carriers_supported <- ifelse(inactive | d$n_phyml_carriers == 0, "not evaluable", paste0(d$n_phyml_supported, " / ", d$n_phyml_carriers))
    d$phyml_support_rate <- ifelse(is.na(d$phyml_support_fraction), NA, percent(d$phyml_support_fraction, accuracy = 0.1))
    show <- d[, c("method","source_class","availability","availability_note","n_tested","n_carriers","carrier_prevalence","prevalence_ci95","phyml_carriers_supported","phyml_support_rate"), drop = FALSE]
    datatable(show, rownames = FALSE, options = list(dom = "t", scrollX = TRUE))
  })
  output$support_trajectory <- renderPlotly({
    d <- locus_trajectory_data(); r <- selected_locus()
    if (!nrow(d) || is.null(r)) return(plotly_empty() %>% layout(annotations = list(list(text = "No locus trajectory; rerun normalize after PhyML and tract callers", showarrow = FALSE))))
    palette <- c(phyml="#2c3e50", ibdmix="#e67e22", trace="#8e44ad", as3="#16883f")
    p <- plot_ly()
    for (method in intersect(names(palette), unique(d$method))) {
      z <- d[d$method == method, , drop = FALSE]
      rgba <- switch(method, phyml="rgba(44,62,80,0.12)", ibdmix="rgba(230,126,34,0.12)", trace="rgba(142,68,173,0.12)", as3="rgba(22,136,63,0.12)")
      z$hover <- paste0(toupper(method), " · ", z$source_class, "<br>chr", r$chr[[1]], ":", format(z$bin_start,big.mark=","), "-", format(z$bin_end,big.mark=","),
                        "<br>carriers: ", z$n_carriers, "/", z$n_denominator, "<br>prevalence: ", percent(z$prevalence,accuracy=.01),
                        "<br>95% CI: ", percent(z$ci_low,accuracy=.01), "–", percent(z$ci_high,accuracy=.01))
      p <- add_ribbons(p, data=z, x=~bin_mid, ymin=~ci_low, ymax=~ci_high, name=paste(method,"95% CI"),
                       fillcolor=rgba, line=list(color="transparent"), hoverinfo="skip", showlegend=FALSE, inherit=FALSE)
      p <- add_lines(p, data=z, x=~bin_mid, y=~prevalence, name=method, text=~hover, hoverinfo="text",
                     line=list(color=unname(palette[[method]]), width=if (method=="phyml") 4 else 2, dash=if (method=="phyml") "dash" else "solid"), inherit=FALSE)
    }
    denominator_label <- if (identical(input$support_denominator,"tested")) "all PhyML-tested" else "PhyML-positive"
    layout(p, hovermode="x unified", xaxis=list(title=paste0("chr",r$chr[[1]]," position (bp)", " · ", current_population(), " · ", denominator_label)),
           yaxis=list(title="Carrier prevalence",tickformat=".1%",range=c(0,max(.02,d$ci_high,na.rm=TRUE))),
           shapes=list(list(type="rect",x0=r$selected_start[[1]],x1=r$selected_end[[1]],y0=0,y1=1,yref="paper",fillcolor="rgba(24,188,156,.08)",line=list(width=0),layer="below")),
           legend=list(orientation="h",x=0,y=1.08),margin=list(t=55))
  })
  locus_sample_query <- function(limit = TRUE) {
    r <- selected_locus(); if (is.null(r)) return(data.frame())
    sql <- "SELECT sample_id,population,super_population,phyml,ibdmix,trace,as3,n_methods,methods_support,matched_haplotypes,best_lineages,matched_source_classes,max_prop_match,max_segment_score,max_posterior FROM locus_sample_support WHERE dataset_id=? AND genome_build=? AND locus_id=?"
    params <- list(current_dataset(),current_build(),r$locus_id[[1]])
    if (current_population() != "ALL") { sql <- paste0(sql," AND population=?"); params <- c(params,current_population()) }
    if (isTRUE(input$support_phyml_only)) sql <- paste0(sql," AND phyml=1")
    sql <- paste0(sql," ORDER BY n_methods DESC,phyml DESC,max_prop_match DESC,sample_id")
    if (isTRUE(limit)) sql <- paste0(sql," LIMIT 2000")
    list(sql=sql,params=params)
  }
  locus_sample_data <- reactive({
    query <- locus_sample_query(TRUE)
    if (is.data.frame(query)) return(query)
    Q(query$sql,query$params)
  })
  output$support_concordance <- renderPlotly({
    r <- selected_locus(); if (is.null(r)) return(plotly_empty())
    sql <- "SELECT methods_support AS support_pattern,n_methods,COUNT(*) AS n_samples FROM locus_sample_support WHERE dataset_id=? AND genome_build=? AND locus_id=? AND phyml=1"
    params <- list(current_dataset(),current_build(),r$locus_id[[1]])
    if (current_population() != "ALL") { sql<-paste0(sql," AND population=?");params<-c(params,current_population()) }
    d <- Q(paste0(sql," GROUP BY support_pattern,n_methods ORDER BY n_methods DESC,n_samples DESC"),params)
    if (!nrow(d)) return(plotly_empty())
    d$support_pattern <- factor(d$support_pattern,levels=rev(d$support_pattern))
    plot_ly(d,x=~n_samples,y=~support_pattern,color=~factor(n_methods),type="bar",orientation="h",
            text=~paste0(support_pattern,"<br>samples: ",format(n_samples,big.mark=",")),hoverinfo="text") %>%
      layout(xaxis=list(title="Individuals"),yaxis=list(title=""),showlegend=FALSE,margin=list(l=135))
  })
  output$locus_sample_support <- renderDT({
    d <- locus_sample_data()
    if (!nrow(d)) return(datatable(data.frame(Message="No matching individuals"),rownames=FALSE,options=list(dom="t")))
    if ("max_prop_match" %in% names(d)) d$max_prop_match <- ifelse(is.na(d$max_prop_match),NA,percent(d$max_prop_match,accuracy=.1))
    datatable(d,rownames=FALSE,options=list(pageLength=20,scrollX=TRUE))
  })
  output$download_locus_support <- downloadHandler(
    filename=function(){r<-selected_locus();paste0(current_dataset(),"-",current_build(),"-",r$locus_id[[1]],"-individual-support.tsv")},
    content=function(file){
      query<-locus_sample_query(FALSE)
      if(is.data.frame(query)){data.table::fwrite(query,file,sep="\t");return(invisible(NULL))}
      result<-DBI::dbSendQuery(con_rv(),query$sql);on.exit(DBI::dbClearResult(result),add=TRUE)
      DBI::dbBind(result,query$params)
      first<-TRUE
      repeat{
        chunk<-DBI::dbFetch(result,n=50000)
        if(!nrow(chunk)){
          if(first)data.table::fwrite(chunk,file,sep="\t")
          break
        }
        data.table::fwrite(chunk,file,sep="\t",append=!first,col.names=first)
        first<-FALSE
      }
    }
  )
  output$download_locus_tracks <- downloadHandler(
    filename=function(){r<-selected_locus();paste0(current_dataset(),"-",current_build(),"-",r$locus_id[[1]],"-prevalence.bedGraph")},
    content=function(file){
      r<-selected_locus();d<-locus_trajectory_data();con<-base::file(file,"wt");on.exit(close(con),add=TRUE)
      for(method in unique(d$method)){
        z<-d[d$method==method,,drop=FALSE]
        writeLines(paste0("track type=bedGraph name=\"GU_",method,"_prevalence\" description=\"",current_dataset()," ",r$locus_id[[1]]," ",current_population()," ",input$support_denominator %||% "phyml_positive","\""),con)
        write.table(data.frame(paste0("chr",r$chr[[1]]),z$bin_start,z$bin_end,round(z$prevalence,6)),con,sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
      }
    }
  )

  loci_data <- reactive({
    req(input$locus_chr, input$locus_status)
    sql <- "SELECT l.rowid,l.locus_id,l.chr,l.start AS selected_start,l.end AS selected_end,l.input_start,l.input_end,l.analysis_start,l.analysis_end,l.flank_bp,l.anchor_pos,l.selection_method,l.ld_r2_threshold,l.n_search_sites,l.n_ld_sites,l.n_ancestral_sites,l.n_sites,l.best_archaic,l.best_lineage,l.n_compared,l.n_match,l.prop_match,l.source,l.source_class,l.status,l.direct_match_pass,l.n_carriers,l.tree_status,l.n_bootstrap_nodes,l.bootstrap_min,l.bootstrap_median,l.bootstrap_max,l.candidate_lineage,l.tree_has_ancestral_outgroup,l.candidate_clade_pass,l.candidate_clade_rule,l.candidate_clade_bootstrap,l.candidate_clade_node,l.candidate_clade_side,l.candidate_clade_n_tips,l.candidate_clade_modern_tips,l.candidate_clade_archaic_tips,l.candidate_clade_specificity,l.candidate_clade_tips,l.tree_newick,l.raw_file,l.genome_build,l.dataset_id,e.sequence_information,e.region_coverage_fraction,e.phyml_qc,e.evidence_state,e.evidence_summary FROM loci l LEFT JOIN locus_evidence e ON e.dataset_id=l.dataset_id AND e.genome_build=l.genome_build AND e.locus_id=l.locus_id WHERE l.dataset_id=? AND l.genome_build=? AND l.method='phyml'"
    params <- list(current_dataset(),current_build())
    if (input$locus_chr != "ALL") { sql <- paste0(sql, " AND l.chr=?"); params <- c(params, input$locus_chr) }
    if (input$locus_status != "ALL") { sql <- paste0(sql, " AND l.status=?"); params <- c(params, input$locus_status) }
    Q(paste0(sql, " ORDER BY CASE WHEN l.chr GLOB '[0-9]*' THEN CAST(l.chr AS INT) ELSE 99 END,l.start"), params)
  })
  # GU_DUAL_LEAD_V2_BEGIN
  gu_phyml_report_server(input, output, session, current_dataset, current_build, .gu_review_root)
  # GU_DUAL_LEAD_V2_END
  output$loci_table <- renderDT({
    d <- loci_data()
    selected <- match(active_locus_rowid(), d$rowid)
    datatable(d[, setdiff(names(d), c("tree_newick", "raw_file")), drop = FALSE],
              selection = list(mode = "single", selected = selected, target = "row"),
              rownames = FALSE, options = list(pageLength = 15, scrollX = TRUE))
  })
  observeEvent(input$loci_table_rows_selected, {
    d <- loci_data(); i <- input$loci_table_rows_selected
    if (length(i) && i >= 1L && i <= nrow(d)) {
      active_locus_rowid(as.integer(d$rowid[[i]])); browser_override(NULL)
      updateSelectInput(session, "browser_locus", selected = as.character(d$rowid[[i]]))
    }
  }, ignoreInit = TRUE)
  phyml_file <- function(name) {
    r <- selected_locus(); if (is.null(r)) return(NULL)
    file.path(dirname(.gu_resolve_artifact(r$raw_file)), name)
  }
  haplotype_view <- reactive({
    r <- selected_locus(); req(!is.null(r))
    tryCatch(
      .gu_read_phyml_view(r, input$hap_n_match %||% 8L, input$hap_n_control %||% 10L, input$hap_max_sites %||% 150L),
      error = function(e) structure(list(message = conditionMessage(e)), class = "gu_haplotype_error")
    )
  })
  output$haplotype_title <- renderUI({
    r <- selected_locus()
    if (is.null(r)) return(tags$div(class = "alert alert-info", "No PhyML locus is available for this genome build."))
    tags$div(
      tags$h5(paste0(r$locus_id[[1]], " · selected tree region chr", r$chr[[1]], ":",
                     format(r$selected_start[[1]], big.mark = ","), "-", format(r$selected_end[[1]], big.mark = ","))),
      tags$p(class = "mb-0", paste0("Input BED: ", format(r$input_start[[1]], big.mark = ","), "-",
                                      format(r$input_end[[1]], big.mark = ","), "; search interval: ",
                                      format(r$analysis_start[[1]], big.mark = ","), "-", format(r$analysis_end[[1]], big.mark = ","),
                                      "; anchor: ", format(r$anchor_pos[[1]], big.mark = ","),
                                      "; LD r²: ", r$ld_r2_threshold[[1]], "; source: ", r$source[[1]], "; status: ", r$status[[1]]))
    )
  })
  output$haplotype_matrix <- renderUI({
    v <- haplotype_view()
    if (inherits(v, "gu_haplotype_error")) return(tags$div(class = "alert alert-warning", v$message))
    base_colours <- c(A = "#16883f", C = "#1769d2", G = "#e67e22", T = "#d62728", N = "#9aa3aa")
    positions <- suppressWarnings(as.integer(v$sites$pos[v$idx]))
    refs <- toupper(as.character(v$sites$ref[v$idx])); alts <- toupper(as.character(v$sites$alt[v$idx]))
    arch_rows <- v$rows[vapply(v$rows, function(x) identical(x$kind, "archaic"), logical(1))]
    arch_lookup <- setNames(lapply(arch_rows, `[[`, "bases"), vapply(arch_rows, `[[`, character(1), "id"))
    header <- tags$tr(
      tags$th(class = "gu-rowlab", "Sequence"),
      lapply(seq_along(positions), function(j) tags$th(class = "gu-site", title = paste0("REF/ALT: ", refs[[j]], "/", alts[[j]]),
                                                        tags$span(format(positions[[j]], big.mark = ","))))
    )
    n_arch <- length(arch_rows)
    row_tags <- lapply(seq_along(v$rows), function(i) {
      row <- v$rows[[i]]; is_arch <- identical(row$kind, "archaic"); is_control <- identical(row$kind, "control")
      best <- arch_lookup[[row$best_archaic]]
      cells <- lapply(seq_along(row$bases), function(j) {
        base <- toupper(row$bases[[j]]); if (!base %in% names(base_colours)) base <- "N"
        is_match <- !is_arch && !is.null(best) && base %in% c("A", "C", "G", "T") && identical(base, best[[j]])
        title <- paste0(row$tooltip, "\nSite: chr", v$sites$chr[v$idx[j]], ":", format(positions[[j]], big.mark = ","),
                        "\nAllele: ", base, "\nREF/ALT: ", refs[[j]], "/", alts[[j]],
                        if (is_match) paste0("\nMatches best archaic: ", row$best_archaic) else "")
        tags$td(class = paste("gu-base", if (is_match) "gu-base-match" else ""),
                style = paste0("color:", base_colours[[base]]), title = title, if (base == "N") "" else base)
      })
      cls <- c(if (i == n_arch) "gu-arch-last", if (is_control) "gu-control",
               if (is_control && !any(vapply(v$rows[seq_len(max(0, i - 1L))], function(x) identical(x$kind, "control"), logical(1)))) "gu-control-first")
      tags$tr(class = paste(cls[nzchar(cls)], collapse = " "), tags$th(class = "gu-rowlab", title = row$tooltip, row$label), cells)
    })
    tags$div(
      tags$div(class = "gu-hap-legend",
               tags$span(style = paste0("color:", base_colours[["A"]]), "A"),
               tags$span(style = paste0("color:", base_colours[["C"]]), "C"),
               tags$span(style = paste0("color:", base_colours[["G"]]), "G"),
               tags$span(style = paste0("color:", base_colours[["T"]]), "T"),
               tags$span(style = "background:#dcf4e4;box-shadow:inset 0 -3px #198754;padding:3px 7px;color:#34495e", "modern = best archaic match"),
               tags$span(style = "color:#a93226", "red label = non-matched control")),
      tags$p(class = "text-muted", paste0("chr", v$sites$chr[[1]], ": ", format(min(positions), big.mark = ","), "-",
                                           format(max(positions), big.mark = ","), "; showing ", v$n_display, " of ", v$n_sites,
                                           " callable SNPs. Hover cells and row labels for provenance.")),
      tags$div(class = "gu-hap-scroll", tags$table(class = "gu-hap-table", tags$thead(header), tags$tbody(row_tags)))
    )
  })
  output$haplotype_similarity <- renderDT({
    v <- haplotype_view()
    if (inherits(v, "gu_haplotype_error")) return(datatable(data.frame(Message = v$message), rownames = FALSE, options = list(dom = "t")))
    matched <- data.table::copy(v$matched); controls <- data.table::copy(v$controls)
    if (nrow(matched)) matched[, display_group := "tree candidate"]
    if (nrow(controls)) controls[, display_group := "outside candidate clade"]
    d <- data.table::rbindlist(list(matched, controls), fill = TRUE)
    keep <- intersect(c("display_group", "hap_id", "n", "best_archaic", "best_lineage", "n_compared", "n_match", "prop_match", "direct_match_pass", "copies"), names(d))
    datatable(d[, ..keep], rownames = FALSE, options = list(pageLength = 20, scrollX = TRUE))
  })
  output$haplotype_note <- renderText({
    v <- haplotype_view()
    if (inherits(v, "gu_haplotype_error")) return(v$message)
    paste0("Run directory: ", v$paths$run, "\nSites: ", v$paths$sites, "\nArchaic sequences: ", v$paths$archaic,
           "\nNormalized haplotypes: ", v$paths$haplotypes, "\nCallable sites: ", v$n_sites, "; displayed: ", v$n_display,
           "\nDisplayed candidate groups: ", nrow(v$matched), "; outside-clade controls: ", nrow(v$controls),
           "\nOnly A/C/G/T are printed; missing or ambiguous bases are blank. Site subsampling is deterministic and evenly spaced.")
  })
  output$haplotype_table <- renderDT({
    r <- selected_locus()
    if (is.null(r)) return(datatable(data.frame(Message = "Select a locus"), rownames = FALSE))
    f <- phyml_file("haplotypes.tsv")
    if (!file.exists(f)) return(datatable(data.frame(Message = paste("Missing", f)), rownames = FALSE))
    d <- data.table::fread(f); d <- d[locus_id == r$locus_id[[1]]]
    keep <- intersect(c("hap_id", "n", "best_archaic", "best_lineage", "n_compared", "n_match", "prop_match", "direct_match_pass", "archaic_match", "seq"), names(d))
    datatable(d[, ..keep], rownames = FALSE, options = list(scrollX = TRUE, pageLength = 15))
  })
  output$haplotype_samples <- renderDT({
    r <- selected_locus()
    if (is.null(r)) return(datatable(data.frame(Message = "Select a locus"), rownames = FALSE))
    f <- phyml_file("haplotype_samples.tsv")
    if (!file.exists(f)) return(datatable(data.frame(Message = paste("Missing", f)), rownames = FALSE))
    d <- data.table::fread(f); d <- d[locus_id == r$locus_id[[1]]]
    tips <- unique(strsplit(as.character(r$candidate_clade_tips[[1]] %||% ""), ",", fixed=TRUE)[[1]])
    if (.gu_truth(r$candidate_clade_pass[[1]]) && "hap_id" %in% names(d)) d <- d[as.character(hap_id) %in% tips] else d <- d[0]
    datatable(d, rownames = FALSE, options = list(pageLength = 15))
  })
  output$tree_summary <- renderText({
    r <- selected_locus()
    if (is.null(r)) return("Select a locus")
    paste0("status=", r$tree_status, "; INFO/AA ancestral outgroup=", ifelse(r$tree_has_ancestral_outgroup==1,"present","absent"),
           "; clade pass=", r$candidate_clade_pass,
           "; rule=", r$candidate_clade_rule, "; candidate lineage=", r$candidate_lineage,
           "; candidate-clade bootstrap=", r$candidate_clade_bootstrap,
           "; clade tips=", r$candidate_clade_n_tips,
           " (modern=", r$candidate_clade_modern_tips, ", archaic=", r$candidate_clade_archaic_tips, ")",
           "\nAll-tree bootstrap min/median/max=", r$bootstrap_min, "/", r$bootstrap_median, "/", r$bootstrap_max,
           ". Candidate bootstrap is branch repeatability, not an introgression probability.")
  })
  output$phyml_tree <- renderPlot({
    r <- selected_locus(); req(!is.null(r), nzchar(r$tree_newick[[1]]))
    tree <- ape::read.tree(text = r$tree_newick[[1]])
    targets <- unique(strsplit(as.character(r$candidate_clade_tips[[1]] %||% ""), ",", fixed=TRUE)[[1]])
    targets <- intersect(targets[!is.na(targets) & nzchar(targets)], tree$tip.label)
    candidate_side <- targets
    outside <- setdiff(tree$tip.label,candidate_side)
    if (length(candidate_side) >= 2L && length(outside)) {
      outgroup <- if ("Ancestral" %in% outside) "Ancestral" else outside[[1]]
      tree <- tryCatch(ape::root(ape::unroot(tree), outgroup=outgroup, resolve.root=FALSE), error=function(e) tree)
    }
    tree <- ape::ladderize(tree)
    archaic <- grepl("Altai|Chagyr|Vindija|Denisova|Neander", tree$tip.label, ignore.case = TRUE)
    target_idx <- match(targets,tree$tip.label); target_idx <- target_idx[!is.na(target_idx)]
    candidate_side_idx <- match(candidate_side,tree$tip.label); candidate_side_idx <- candidate_side_idx[!is.na(candidate_side_idx)]
    candidate_node <- if (length(candidate_side_idx)>=2L) ape::getMRCA(tree,candidate_side_idx) else NA_integer_
    candidate_desc <- if (is.finite(candidate_node)) .gu_descendant_tips(tree,candidate_node) else integer()
    tip_col <- rep("#68757f",ape::Ntip(tree)); tip_col[target_idx] <- "#b8322a"; tip_col[archaic] <- "#111111"
    show_labels <- ape::Ntip(tree)<=500L
    par(mar=c(1.5,1,1,8),xpd=NA)
    plot(tree,type="phylogram",direction="rightwards",show.node.label=FALSE,show.tip.label=show_labels,
         tip.color=tip_col,cex=max(.25,min(.8,25/sqrt(ape::Ntip(tree)))),edge.color="#42484d",no.margin=FALSE)
    if(length(candidate_desc) && is.finite(candidate_node)) {
      pp <- get("last_plot.phylo",envir=.PlotPhyloEnv)
      rect(pp$xx[candidate_node],min(pp$yy[candidate_desc])-.4,max(pp$xx,na.rm=TRUE)*1.015,max(pp$yy[candidate_desc])+.4,
           col=adjustcolor("#ef9a9a",alpha.f=.28),border=NA)
      .gu_redraw_phylogram_edges(tree,pp,"#42484d",1)
    }
    ape::tiplabels(pch=16,col=tip_col,cex=.35,frame="none")
    if(!show_labels && length(target_idx)) ape::tiplabels(text=tree$tip.label[target_idx],tip=target_idx,frame="none",adj=c(-.05,.5),cex=.65,col="#b8322a")
    candidate_bootstrap <- suppressWarnings(as.numeric(r$candidate_clade_bootstrap[[1]]))
    if(is.finite(candidate_node) && is.finite(candidate_bootstrap)) ape::nodelabels(candidate_bootstrap,node=candidate_node,frame="rect",bg="white",cex=.7,col="#8e1b17")
    try(ape::add.scale.bar(cex=.65,lwd=.8),silent=TRUE)
  })

  browser_target <- reactive({
    z <- browser_override()
    if (!is.null(z)) return(z)
    r <- selected_locus(); if (is.null(r)) return(NULL)
    data.frame(chr = r$chr[[1]], start = r$selected_start[[1]], end = r$selected_end[[1]],
               genome_build = r$genome_build[[1]], label = paste0("PhyML · ", r$locus_id[[1]]), stringsAsFactors = FALSE)
  })
  output$genome_browser <- renderUI({
    r <- browser_target()
    if (is.null(r)) return(tags$div(class = "alert alert-info", "No locus or segment is available for this build."))
    u <- .gu_browser_urls(r, as.integer(input$browser_flank %||% 250000L))
    tags$div(
      tags$div(class = "gu-browser-toolbar",
               tags$span(class = "gu-locus-label", paste0(u$label, " · ", u$build, " · ", u$locus)),
               tags$a("Open IGV ↗", href = u$igv, target = "_blank", rel = "noopener", class = "btn btn-sm btn-primary"),
               tags$a("Open UCSC ↗", href = u$ucsc, target = "_blank", rel = "noopener", class = "btn btn-sm btn-outline-primary")),
      # IGV-Web loads workers and browser-managed resources dynamically.  A
      # sandboxed iframe breaks those paths in several Chrome releases; retain
      # the same unsandboxed embed contract as the previously working GU UI.
      tags$iframe(src = u$igv, class = "gu-igv-frame", title = "IGV-Web genome browser")
    )
  })
  observeEvent(list(current_build(), current_dataset()), browser_override(NULL), ignoreInit = TRUE)

  segment_data <- eventReactive(list(input$seg_go, input$genome_build, input$target_dataset), ignoreNULL = FALSE, {
    req(input$seg_method, input$seg_chr, input$seg_limit, input$seg_end > input$seg_start)
    sql <- "SELECT sample_id,method,source_class AS source,chr,start,end,length_bp,haplotype,score,posterior,locus_id,batch_id FROM segments WHERE dataset_id=? AND genome_build=? AND end>? AND start<?"
    params <- list(current_dataset(), current_build(), input$seg_start, input$seg_end)
    if (input$seg_method != "ALL") { sql <- paste0(sql, " AND method=?"); params <- c(params, input$seg_method) }
    if (input$seg_chr != "ALL") { sql <- paste0(sql, " AND chr=?"); params <- c(params, input$seg_chr) }
    Q(paste0(sql, " ORDER BY chr,start LIMIT ", as.integer(input$seg_limit)), params)
  })
  reference_data <- eventReactive(list(input$seg_go,input$genome_build,input$show_reference,input$reference_population), ignoreNULL=FALSE, {
    if (!isTRUE(input$show_reference)) return(data.frame())
    req(input$seg_chr, input$reference_population, input$seg_limit, input$seg_end > input$seg_start)
    sql <- "SELECT dataset_id,population,source_class AS source,reference_role,chr,start,end,end-start AS length_bp,raw_file FROM reference_callsets WHERE genome_build=? AND end>? AND start<?"
    params <- list(current_build(),input$seg_start,input$seg_end)
    if (input$seg_chr != "ALL") { sql<-paste0(sql," AND chr=?");params<-c(params,input$seg_chr) }
    if (input$reference_population != "ALL") { sql<-paste0(sql," AND population=?");params<-c(params,input$reference_population) }
    Q(paste0(sql," ORDER BY population,chr,start LIMIT ",as.integer(input$seg_limit)),params)
  })
  output$segment_table <- renderDT(datatable(segment_data(), selection = "single", rownames = FALSE, options = list(scrollX = TRUE, pageLength = 20)))
  observeEvent(input$segment_table_rows_selected, {
    d <- segment_data(); i <- input$segment_table_rows_selected
    if (length(i) && i >= 1L && i <= nrow(d)) {
      browser_override(data.frame(chr = d$chr[[i]], start = d$start[[i]], end = d$end[[i]],
                                  genome_build = current_build(),
                                  label = paste0(d$method[[i]], " · ", d$sample_id[[i]] %||% "sample", " · ", d$source[[i]]),
                                  stringsAsFactors = FALSE))
    }
  }, ignoreInit = TRUE)
  output$published_callset_table <- renderDT({
    d<-reference_data(); if(!nrow(d))return(datatable(data.frame(Message="No matching external reference intervals"),rownames=FALSE))
    datatable(d,rownames=FALSE,options=list(scrollX=TRUE,pageLength=15))
  })
  output$reference_overlap_table <- renderDT({
    req(input$seg_chr, input$reference_population, input$seg_end > input$seg_start)
    sql <- "SELECT reference_population AS population,COUNT(*) AS n_overlaps,COUNT(DISTINCT sample_id) AS n_samples,ROUND(AVG(result_overlap),4) AS mean_result_overlap,ROUND(AVG(reference_overlap),4) AS mean_reference_overlap,ROUND(AVG(reciprocal_overlap),4) AS mean_reciprocal_overlap,SUM(overlap_bp) AS overlap_bp FROM reference_callset_overlaps WHERE result_dataset_id=? AND genome_build=? AND segment_end>? AND segment_start<?"
    params<-list(current_dataset(),current_build(),input$seg_start,input$seg_end)
    if(input$seg_chr!="ALL"){sql<-paste0(sql," AND chr=?");params<-c(params,input$seg_chr)}
    if(input$reference_population!="ALL"){sql<-paste0(sql," AND reference_population=?");params<-c(params,input$reference_population)}
    d<-Q(paste0(sql," GROUP BY reference_population ORDER BY reference_population"),params)
    datatable(d,rownames=FALSE,options=list(pageLength=15))
  })
  output$segment_plot <- renderPlotly({
    d <- segment_data(); r<-reference_data(); if (!nrow(d) && !nrow(r)) return(plotly_empty())
    p<-plot_ly()
    if(nrow(d)){
      d$row<-seq_len(nrow(d));d$label<-paste(d$sample_id,d$method,d$source,sep=" | ")
      p<-add_segments(p,data=d,x=~start,xend=~end,y=~row,yend=~row,color=~method,text=~label,hoverinfo="text",line=list(width=5))
    }
    if(nrow(r)){
      offset<-if(nrow(d))nrow(d) else 0L;r$row<-offset+seq_len(nrow(r));r$label<-paste("external_reference",r$dataset_id,r$population,r$source,sep=" | ")
      p<-add_segments(p,data=r,x=~start,xend=~end,y=~row,yend=~row,name="published callset",text=~label,hoverinfo="text",line=list(width=5,color="#d35400",dash="dash"),inherit=FALSE)
    }
    layout(p,xaxis=list(title="Position (bp)"),yaxis=list(title="Calls + external reference",showticklabels=FALSE))
  })
  sample_data <- reactive({
    req(input$burden_type, input$sample_method)
    sql <- "SELECT b.sample_id,p.population,p.super_population,b.method,b.source_class AS source,b.burden_type,b.n_input_segments,b.n_merged_intervals,b.n_chromosomes,b.total_bp,b.dosage_bp FROM sample_burden b LEFT JOIN sample_populations p ON p.dataset_id=b.dataset_id AND p.sample_id=b.sample_id WHERE b.dataset_id=? AND b.genome_build=? AND b.burden_type=?"
    params <- list(current_dataset(),current_build(),input$burden_type)
    if (isTRUE(nzchar(input$sample_search))) { sql <- paste0(sql, " AND b.sample_id LIKE ?"); params <- c(params, paste0("%", input$sample_search, "%")) }
    if (input$sample_method != "ALL") { sql <- paste0(sql, " AND b.method=?"); params <- c(params, input$sample_method) }
    Q(paste0(sql, " ORDER BY b.total_bp DESC LIMIT 10000"), params)
  })
  output$sample_table <- renderDT(datatable(sample_data(), rownames = FALSE, options = list(pageLength = 25)))
  output$db_info <- renderText({
    paste("Database:", db_path, "\nSize:", format(file.info(db_path)$size, big.mark = ","), "bytes\nTarget dataset:", current_dataset(), "\nGenome build:", current_build(), "\nMethods:", paste(methods(), collapse = ", "), "\nExternal reference intervals:", SQ("SELECT COUNT(*) FROM reference_callsets WHERE genome_build=? AND reference_role='external_reference'",list(current_build())))
  })
  output$download_segments <- downloadHandler(
    filename = function() paste0("gu-segments-", current_dataset(), "-", current_build(), "-", Sys.Date(), ".tsv"),
    content = function(file) data.table::fwrite(segment_data(), file, sep = "\t")
  )
}

})

shiny::runApp(shiny::shinyApp(ui, server),
  host = Sys.getenv("GU_SHINY_HOST", "127.0.0.1"),
  port = as.integer(Sys.getenv("GU_SHINY_PORT", "3838")), launch.browser = interactive())
