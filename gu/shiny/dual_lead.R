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

