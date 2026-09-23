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
    CHM13 = list(igv = "hs1", ucsc = "hs1"),
    list(igv = "hg19", ucsc = "hg19")
  )
  chrom <- sub("^chr", "", as.character(target$chr[[1]]), ignore.case = TRUE)
  start <- max(1L, suppressWarnings(as.integer(target$start[[1]])) + 1L - as.integer(flank))
  end <- suppressWarnings(as.integer(target$end[[1]])) + as.integer(flank)
  if(end-start+1 > IGV_region_max) {
    mid <- (start+end)/2
    start <- max(1,floor(mid-IGV_region_max/2)+1)
    end <- start+IGV_region_max-1
  }
  locus <- paste0("chr", chrom, ":", start, "-", end)
  list(
    label = as.character(target$label[[1]] %||% "Selected region"),
    build = build, locus = locus, genome = ids$igv,
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
  if("role" %in% names(hap))hap[, .pass := role=="risk"]
  data.table::setorder(hap, -.pass, -.n, hap_id, na.last = TRUE)
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
    d <- if (has_summary()) {
      Q("SELECT DISTINCT dataset_id FROM database_summary WHERE genome_build=? AND dataset_id!='reference' AND metric IN ('segment_calls','phyml_loci') AND value>0 ORDER BY dataset_id", list(current_build()))
    } else Q("SELECT dataset_id FROM segments WHERE genome_build=? UNION SELECT dataset_id FROM loci WHERE genome_build=? ORDER BY dataset_id", list(current_build(), current_build()))
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
  active_lead_record <- reactiveVal(NULL)
  linked_region <- reactiveVal(NULL)
  selected_matching_record <- reactiveVal(NULL)
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
      if(!identical(isolate(active_locus_rowid()),z)) {
        active_locus_rowid(z); browser_override(NULL)
        r<-isolate(selected_locus());linked_region(data.frame(chr=r$chr,start=r$selected_start,end=r$selected_end))
      }
    }
  }, ignoreInit = TRUE)
  observe({
    d <- all_loci_data()
    choices <- if(nrow(d)) setNames(as.character(d$rowid),paste0(d$locus_id," · chr",d$chr)) else character()
    updateSelectInput(session,"phyml_locus",choices=choices,selected=as.character(active_locus_rowid()))
  })
  observeEvent(input$phyml_locus, {
    id <- suppressWarnings(as.integer(input$phyml_locus)); d <- all_loci_data()
    i <- match(id,d$rowid); req(length(i)==1L,!is.na(i))
    if(!identical(active_locus_rowid(),id)) {
      active_locus_rowid(id); browser_override(NULL)
      linked_region(data.frame(chr=d$chr[[i]],start=d$selected_start[[i]],end=d$selected_end[[i]]))
      updateSelectInput(session,"browser_locus",selected=as.character(id))
    }
  },ignoreInit=TRUE)
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
  density_root <- Sys.getenv("GU_DENSITY_DIR", file.path(dirname(db_path),"normalize","density"))
  density_stamp <- reactivePoll(3000,session,
    checkFunc=function() {
      paths <- Sys.glob(file.path(density_root,"*","current.tsv"))
      paste(paths,as.numeric(file.info(paths)$mtime),collapse=";")
    },valueFunc=function() Sys.time())
  density_data <- reactive({
    slug <- gsub("[^A-Za-z0-9_.-]", "_", paste(current_dataset(),current_build(),sep="."))
    root <- file.path(Sys.getenv("GU_DENSITY_DIR", file.path(dirname(db_path),"normalize","density")),slug)
    pointer <- file.path(root,"current.tsv")
    # Database replacement invalidates this cache; no interval scan occurs in R.
    db_stamp(); density_stamp()
    validate(need(file.exists(pointer),"尚无 IBDmix 全基因组密度缓存；运行 gu.sh shiny 可生成。"))
    version <- data.table::fread(pointer)$directory[[1]]
    path <- file.path(root,version)
    manifest <- jsonlite::fromJSON(file.path(path,"manifest.json"))
    validate(need(!is.null(manifest$schema) && manifest$schema>=11,"请重新运行 gu.sh shiny，生成五参考及 Altai + Denisova、All Five 合并统计。"))
    info<-file.info(db_path)
    validate(need(!is.null(manifest$mtime_ns) && !is.null(manifest$size) &&
      abs(as.numeric(info$mtime)*1e9-manifest$mtime_ns)<1e6 && info$size==manifest$size,
      "密度缓存与当前数据库不匹配；请在 IBDmix 完成后运行 gu.sh final，再更新 Shiny。"))
    samples <- as.data.frame(data.table::fread(file.path(path,"samples.tsv")))
    bins <- as.data.frame(data.table::fread(file.path(path,"bins.tsv")))
    read_matrix <- function(filename) {
      con <- gzfile(file.path(path,filename),open="rt")
      on.exit(close(con))
      mat <- as.matrix(data.table::fread(text=paste(readLines(con,warn=FALSE),collapse="\n"),header=FALSE,na.strings=c("nan","NA")))
      validate(need(nrow(mat)==nrow(samples) && ncol(mat)==nrow(bins),"密度缓存维度不一致，请重新生成。"))
      mat
    }
    matrices<-list(Neanderthal=read_matrix("matrix.tsv.gz"),Denisovan=read_matrix("denisovan_matrix.tsv.gz"))
    summary_file <- file.path(path,"neanderthal_summary.tsv")
    archaic_summary_file <- file.path(path,"archaic_summary.tsv")
    list(samples=samples,bins=bins,z=matrices$Neanderthal,matrices=matrices,
         summary=if(file.exists(summary_file))as.data.frame(data.table::fread(summary_file,na.strings=c("","NA"))) else NULL,
         archaic_summary=if(file.exists(archaic_summary_file))as.data.frame(data.table::fread(archaic_summary_file,na.strings=c("","NA"))) else NULL,
         manifest=manifest)
  })
  ibdmix_summary <- reactive({
    d<-density_data()
    validate(need(!is.null(d$summary),"尚无 Neanderthal 汇总缓存；重新运行 gu.sh shiny 即可生成。"))
    gu_ibdmix_summary(d$samples,d$summary)
  })
  output$ibdmix_summary_map <- renderPlot({
    gu_ibdmix_map(ibdmix_summary(),.gu_summary_land,.gu_summary_locations)
  },res=120,alt="Static world map of average Neanderthal sequence coverage by 1000 Genomes population. Red pie sectors show percentages; N/A denotes untested populations.")
  output$ibdmix_summary_text <- renderUI({
    gu_ibdmix_summary_text(ibdmix_summary(),density_data()$manifest)
  })
  output$ibdmix_summary_table <- renderUI({
    d<-density_data()
    validate(need(!is.null(d$archaic_summary),"请重新运行 gu.sh shiny，生成五参考及两种合并统计的汇总缓存。"))
    gu_ibdmix_summary_table(ibdmix_summary(),d$archaic_summary)
  })
  density_selection <- reactiveVal(NULL)
  density_view <- reactiveVal(NULL) # NULL = All; explicit region = shared viewport
  density_geometry <- reactive(gu_density_geometry(density_data()$bins))
  output$density_chromosomes <- renderUI({
    r<-density_view(); active<-if(is.null(r)) "All" else as.character(r$chr[[1]])
    tagList(lapply(c("All",as.character(1:22),"X"),function(ch)
      tags$button(type="button",class=paste("btn btn-sm",if(ch==active)"btn-primary" else "btn-link"),
        onclick=sprintf("Shiny.setInputValue('density_chr','%s',{priority:'event'})",ch),ch)),
      actionButton("density_zoom_out","− Zoom out",class="btn-sm"))
  })
  set_overview_region <- function(r) {
    linked_region(r[,c("chr","start","end"),drop=FALSE])
    browser_override(NULL)
  }
  observeEvent(input$density_chr, {
    ch<-input$density_chr; g<-density_geometry()
    if(ch=="All") {density_view(NULL);density_selection(NULL);linked_region(NULL);return()}
    req(ch %in% names(g$lengths))
    set_overview_region(data.frame(chr=ch,start=0,end=g$lengths[[ch]]))
  })
  observeEvent(linked_region(), {
    r<-linked_region();req(!is.null(r))
    density_view(r);density_selection(r)
  },ignoreNULL=TRUE)
  observeEvent(input$density_zoom_out, {
    r<-density_view();req(!is.null(r));g<-density_geometry();ch<-as.character(r$chr[[1]])
    req(ch %in% names(g$lengths));len<-g$lengths[[ch]]
    if(r$end-r$start>=len) {density_view(NULL);density_selection(NULL);linked_region(NULL);return()}
    width<-min(len,2*(r$end-r$start));mid<-(r$start+r$end)/2
    start<-max(0,min(len-width,floor(mid-width/2)))
    set_overview_region(data.frame(chr=ch,start=start,end=start+width))
  })
  density_plot_data <- reactive({
    d<-density_data();r<-density_view()
    if(!is.null(r) && r$end[[1]]-r$start[[1]]<=1e7) {
      regional<-lapply(names(d$matrices),function(lineage) {
        ld<-d;ld$z<-d$matrices[[lineage]]
        gu_density_region(ld,r,Q,current_dataset(),current_build(),lineage)
      })
      d$bins<-regional[[1]]$bins
      d$matrices<-setNames(lapply(regional,function(x)x$z),names(d$matrices))
      return(d)
    }
    if(!is.null(r)) {
      ix<-which(as.character(d$bins$chr)==as.character(r$chr[[1]]) & d$bins$end>r$start[[1]] & d$bins$start<r$end[[1]])
      d$bins<-d$bins[ix,,drop=FALSE];d$matrices<-lapply(d$matrices,function(z)z[,ix,drop=FALSE])
    }
    d
  })
  output$introgression_density <- renderPlotly({
    d<-density_plot_data();bins<-d$bins;r<-density_view();g<-density_geometry()
    populations<-c("Ref","AFR","EAS","EUR","SAS","AMR")
    colours<-c(Neanderthal="#2878B5",Denisovan="#8B5A2B")
    choice<-input$density_lineage %||% "all"
    lineages<-if(choice=="all")names(colours)else intersect(choice,names(colours))
    chromosomes<-intersect(g$chromosomes,unique(as.character(bins$chr)))
    offsets<-if(is.null(r))g$offsets else setNames(rep(0,length(g$offsets)),names(g$offsets))
    x<-(offsets[as.character(bins$chr)]+(bins$start+bins$end)/2)/1e6
    panels<-lapply(seq_along(populations),function(i) {
      pop<-populations[[i]]
      keep<-if(pop=="Ref")d$samples$population %in% "YRI" else d$samples$super_population %in% pop
      p<-plot_ly(source="density")
      unavailable<-character()
      for(lineage in lineages) {
        z<-d$matrices[[lineage]][keep,,drop=FALSE]
        means<-if(nrow(z))colMeans(z,na.rm=TRUE) else rep(NA_real_,nrow(bins))
        means[!is.finite(means)]<-NA_real_
        colour<-if(pop=="Ref")"#444444"else colours[[lineage]]
        label<-if(lineage=="Denisovan")"Denisova"else lineage
        if(!any(is.finite(means)))unavailable<-c(unavailable,label)
        for(ch in chromosomes) {
        ix<-which(as.character(bins$chr)==ch)
        keys<-paste(ch,bins$start[ix],bins$end[ix],sep=":")
        p<-add_trace(p,x=x[ix],y=means[ix],customdata=keys,type="scatter",mode="lines+markers",
          line=list(color=colour,width=1.5,dash=if(pop=="Ref" && lineage=="Denisovan")"dot"else "solid"),marker=list(color=colour,size=4),
          name=label,legendgroup=lineage,showlegend=(i==2L && ch==chromosomes[[1]]),connectgaps=FALSE,
          text=paste0(if(pop=="Ref")"Ref (YRI)" else pop," · ",label," · chr",ch,":",bins$start[ix]+1,"–",bins$end[ix]),
          hovertemplate="%{text}<br>Mean coverage: %{y:.3f}%<extra></extra>")
      }
      }
      axis<-if(is.null(r)) list(title="Chromosome",tickmode="array",tickvals=(g$offsets+g$lengths/2)/1e6,ticktext=g$chromosomes,range=c(0,sum(g$lengths)/1e6),showgrid=FALSE) else
        list(title=paste0("chr",r$chr[[1]]," position (Mb)"),range=c(r$start[[1]],r$end[[1]])/1e6,showgrid=FALSE)
      layout(p,yaxis=list(title=paste0(if(pop=="Ref")"Ref · YRI" else pop," (%)"),rangemode="tozero"),xaxis=axis,
        annotations=if(length(unavailable))list(list(x=1,y=1,xref="paper",yref="paper",text=paste0(paste(unavailable,collapse=" / ")," · N/A"),showarrow=FALSE,xanchor="right",yanchor="top",font=list(size=10,color="#888888")))else list())
    })
    p<-subplot(panels,nrows=length(populations),shareX=TRUE,titleY=TRUE,margin=.018) %>%
      layout(margin=list(l=85,r=25,b=55,t=35),legend=list(orientation="h",x=0,y=1.04),dragmode="zoom",paper_bgcolor="#ffffff",plot_bgcolor="#ffffff") %>%
      config(displaylogo=FALSE,doubleClick=FALSE,scrollZoom=TRUE,modeBarButtonsToRemove=c("select2d","lasso2d"))
    key<-if(is.null(r))NULL else paste(r$chr[[1]],r$start[[1]],r$end[[1]],sep=":")
    htmlwidgets::onRender(p,sprintf("function(el){window.guDensityBind(el,%s);}",jsonlite::toJSON(key,auto_unbox=TRUE,null="null")))
  })
  observeEvent(input$density_bin_click, {
    fields<-strsplit(input$density_bin_click,":",fixed=TRUE)[[1]]
    req(length(fields)==3L)
    r<-data.frame(chr=fields[1],start=as.numeric(fields[2]),end=as.numeric(fields[3]))
    req(all(is.finite(c(r$start,r$end))),r$end>r$start)
    set_overview_region(r)
  })
  observeEvent(input$density_bin_open, {
    fields<-strsplit(input$density_bin_open,":",fixed=TRUE)[[1]]
    req(length(fields)==3L)
    r<-data.frame(chr=fields[1],start=as.numeric(fields[2]),end=as.numeric(fields[3]),super_population="ALL")
    g<-density_geometry()
    req(r$chr %in% names(g$lengths),all(is.finite(c(r$start,r$end))),r$start>=0,r$end>r$start,r$end<=g$lengths[[r$chr]])
    open_matching_record(r,"ibdmix")
  })
  observeEvent(input$density_range, {
    e<-input$density_range;r<-isolate(density_view());g<-density_geometry()
    if(isTRUE(e$reset)) {
      if(is.null(r))return()
      set_overview_region(data.frame(chr=r$chr,start=0,end=g$lengths[[as.character(r$chr[[1]])]]));return()
    }
    lo<-as.numeric(e$start)*1e6;hi<-as.numeric(e$end)*1e6
    req(length(lo)==1,length(hi)==1,is.finite(lo),is.finite(hi),hi>lo)
    if(is.null(r)) {
      # A shared viewport has one chromosome. A cross-chromosome drag selects
      # the chromosome at its midpoint, clipped to that chromosome's bounds.
      mid<-max(0,min(sum(g$lengths)-1,(lo+hi)/2))
      ch<-names(g$offsets)[which(g$offsets<=mid & g$offsets+g$lengths>mid)[1]]
      lo<-lo-g$offsets[[ch]];hi<-hi-g$offsets[[ch]]
    } else ch<-as.character(r$chr[[1]])
    lo<-max(0,min(g$lengths[[ch]]-1,floor(lo)));hi<-max(lo+1,min(g$lengths[[ch]],ceiling(hi)))
    set_overview_region(data.frame(chr=ch,start=lo,end=hi))
  })
  observeEvent(input$igv_region, {
    e<-input$igv_region;ch<-sub("^chr","",e$chr);g<-density_geometry()
    req(ch %in% names(g$lengths),is.finite(e$start),is.finite(e$end),e$end>e$start)
    start<-max(0,min(g$lengths[[ch]]-1,floor(e$start)));end<-min(g$lengths[[ch]],ceiling(e$end))
    req(end>start)
    set_overview_region(data.frame(chr=ch,start=start,end=end))
  })
  for(method in c("ibdmix","trace","as3")) gu_method_server(method,Q,current_dataset,current_build)
  observeEvent(list(current_build(),current_dataset()),{density_view(NULL);density_selection(NULL)},ignoreInit=TRUE)
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
  report_selection <- gu_phyml_report_server(input, output, session, current_dataset, current_build, .gu_review_root,
    active_lead_record, linked_region,
    external_lineage=reactive(sub("^region_","",input$phyml_tree_choice)),
    external_locus=selected_locus)
  observeEvent(list(input$report_locus_select,input$report_locus_open), {
    req(!is.null(input$report_locus_select) || !is.null(input$report_locus_open))
    r <- report_selection(); req(nrow(r)>0)
    d <- all_loci_data(); i <- match(r$locus_id[[1]],d$locus_id); req(!is.na(i))
    active_lead_record(r$record_id[[1]])
    active_locus_rowid(as.integer(d$rowid[[i]])); browser_override(NULL)
    updateSelectInput(session,"phyml_tree_choice",selected=r$lineage[[1]])
    linked_region(data.frame(chr=d$chr[[i]],start=d$selected_start[[i]],end=d$selected_end[[i]]))
    selected_matching_record(NULL);density_selection(NULL)
    updateSelectInput(session,"browser_locus",selected=as.character(d$rowid[[i]]))
  })
  observeEvent(input$report_locus_open, {
    bslib::nav_select("nav","phyml",session=session)
  })
  # GU_DUAL_LEAD_V2_END
  output$loci_table <- renderDT({
    d <- loci_data()
    selected <- match(active_locus_rowid(), d$rowid)
    show <- data.frame(locus=d$locus_id, chr=d$chr,
      region=paste0(format(d$selected_start,big.mark=","),"–",format(d$selected_end,big.mark=",")),
      SNPs=d$n_sites, lineage=d$candidate_lineage, bootstrap=d$candidate_clade_bootstrap,
      carriers=d$n_carriers, conclusion=gsub("_"," ",d$evidence_state),check.names=FALSE)
    datatable(show,
              selection = list(mode = "single", selected = selected, target = "row"),
              rownames = FALSE, options = list(paging = FALSE, scrollX = TRUE, dom="t"))
  })
  observeEvent(input$loci_table_rows_selected, {
    d <- loci_data(); i <- input$loci_table_rows_selected
    if (length(i) && i >= 1L && i <= nrow(d)) {
      if(identical(active_locus_rowid(),as.integer(d$rowid[[i]])))return()
      active_locus_rowid(as.integer(d$rowid[[i]])); browser_override(NULL)
      updateSelectInput(session, "browser_locus", selected = as.character(d$rowid[[i]]))
      linked_region(data.frame(chr=d$chr[[i]],start=d$selected_start[[i]],end=d$selected_end[[i]]))
    }
  }, ignoreInit = TRUE)
  phyml_file <- function(name) {
    r <- selected_locus(); if (is.null(r)) return(NULL)
    file.path(dirname(.gu_resolve_artifact(r$raw_file)), name)
  }
  haplotype_view <- reactive({
    r <- selected_tree(); req(!is.null(r))
    tryCatch(
      .gu_read_phyml_view(r, input$hap_n_match %||% 8L, input$hap_n_control %||% 10L, input$hap_max_sites %||% 150L),
      error = function(e) structure(list(message = conditionMessage(e)), class = "gu_haplotype_error")
    )
  })
  output$haplotype_title <- renderUI({
    r <- selected_locus()
    if (is.null(r)) return(tags$span("暂无 PhyML 位点"))
    tags$div(tags$b(r$locus_id[[1]]), " · chr", r$chr[[1]], ":",
      format(r$selected_start[[1]]+1,big.mark=",",scientific=FALSE), "–",
      format(r$selected_end[[1]],big.mark=",",scientific=FALSE))
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
               tags$span(style = "color:#a93226", "红色标签 = 非风险对照（序列面板）")),
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
    if (nrow(matched)) matched[, display_group := "GWAS risk"]
    if (nrow(controls)) controls[, display_group := "nonrisk / unresolved"]
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
    r <- selected_tree()
    if (is.null(r)) return(datatable(data.frame(Message = "Select a locus"), rownames = FALSE))
    f <- phyml_file("haplotype_samples.tsv")
    if (!file.exists(f)) return(datatable(data.frame(Message = paste("Missing", f)), rownames = FALSE))
    d <- data.table::fread(f); d <- d[locus_id == r$locus_id[[1]]]
    tips <- unique(strsplit(as.character(r$candidate_clade_tips[[1]] %||% ""), ",", fixed=TRUE)[[1]])
    if("role" %in% names(d)) d <- d[role == "risk"] else if (.gu_truth(r$candidate_clade_pass[[1]]) && "hap_id" %in% names(d)) d <- d[as.character(hap_id) %in% tips] else d <- d[0]
    datatable(d, rownames = FALSE, options = list(pageLength = 15))
  })
  available_locus_trees <- reactive({
    r <- selected_locus(); if (is.null(r)) return(data.frame())
    read_trees <- function(name) {
      path <- phyml_file(name)
      if (is.null(path) || !file.exists(path)) return(data.frame())
      d <- as.data.frame(data.table::fread(path, na.strings=c("", "NA")))
      if (!nrow(d)) return(d)
      d[d$locus_id == r$locus_id[[1]], , drop=FALSE]
    }
    d <- read_trees("evidence_trees.tsv")
    if (nrow(d)) d$choice_id <- d$expected_lineage
    region <- read_trees("region_common_trees.tsv")
    if (nrow(region)) {
      region$choice_id <- paste0("region_",region$lineage_filter)
      region$expected_lineage <- region$lineage_filter
      region$candidate_clade_pass <- 0L
      region$candidate_clade_bootstrap <- NA_real_
      region$tree_call_reason <- "regional_tree_for_inspection"
      region$tree_has_ancestral_outgroup <- as.integer(grepl("Ancestral",region$tree_newick,fixed=TRUE))
      d <- as.data.frame(data.table::rbindlist(list(d,region),fill=TRUE))
    }
    d
  })
  observeEvent(available_locus_trees(), {
    d <- available_locus_trees(); r <- selected_locus()
    if (!nrow(d)) {
      updateSelectInput(session, "phyml_tree_choice", choices=c("Locus representative"="representative"))
    } else {
      bs <- d$candidate_clade_bootstrap
      mixed <- d$tree_call_reason=="supported_mixed_archaic_edge_lineage_unresolved"
      if ("mixed_lineage_bootstrap" %in% names(d)) bs[which(mixed)] <- d$mixed_lineage_bootstrap[which(mixed)]
      labels <- paste0(d$expected_lineage, " · ", ifelse(d$candidate_clade_pass==1,"supported",ifelse(d$tree_call_reason=="supported_mixed_archaic_edge_lineage_unresolved","mixed archaic; lineage unresolved","unconfirmed")),
                       " · BS ", ifelse(is.na(bs),"NA",bs))
      regional <- grepl("^region_",d$choice_id)
      labels[regional] <- paste0("Full-region tree · ",d$expected_lineage[regional]," · inspect topology")
      report <- isolate(report_selection())
      lineage <- if(nrow(report) && identical(report$locus_id[[1]],r$locus_id[[1]])) report$lineage[[1]] else r$candidate_lineage[[1]]
      chosen <- match(lineage, d$choice_id)
      if (is.na(chosen)) chosen <- 1L
      updateSelectInput(session, "phyml_tree_choice", choices=setNames(d$choice_id,labels), selected=d$choice_id[[chosen]])
    }
  })
  selected_tree <- reactive({
    r <- selected_locus(); if (is.null(r)) return(NULL)
    d <- available_locus_trees()
    choice <- input$phyml_tree_choice %||% r$candidate_lineage[[1]]
    i <- match(choice, d$choice_id)
    if (!is.na(i)) {
      t <- d[i,,drop=FALSE]
      for (field in intersect(names(r), names(t))) r[[field]] <- t[[field]]
      r$candidate_lineage <- t$expected_lineage
      if (grepl("^region_",choice)) {
        r$candidate_clade_pass <- 0L
        r$candidate_clade_tips <- ""
        r$candidate_edge_tips <- ""
        r$candidate_clade_rule <- "regional_tree_for_inspection; not an automatic introgression call"
        return(r)
      }
      r$candidate_clade_modern_tips <- t$n_candidate_tips_in_clade
      r$candidate_clade_archaic_tips <- t$n_expected_archaic_tips_in_clade
      r$candidate_clade_tips <- paste(na.omit(c(t$candidate_tips_in_clade,t$expected_archaic_tips_in_clade)),collapse=",")
      r$candidate_clade_rule <- t$tree_call_reason
      r$candidate_edge_tips <- paste(na.omit(c(t$candidate_tips_in_clade,t$control_tips_in_clade,t$candidate_context_tips_in_clade,t$expected_archaic_tips_in_clade)),collapse=",")
      r$mixed_lineage_bootstrap <- if ("mixed_lineage_bootstrap" %in% names(t)) t$mixed_lineage_bootstrap else NA_real_
      if (!.gu_truth(t$candidate_clade_pass[[1]]) && "mixed_lineage_edge_tips" %in% names(t) && !is.na(t$mixed_lineage_edge_tips[[1]])) {
        r$candidate_edge_tips <- t$mixed_lineage_edge_tips
        r$candidate_clade_tips <- paste(na.omit(c(t$mixed_lineage_candidate_tips,t$mixed_lineage_archaic_tips)),collapse=",")
      }
    }
    r
  })
  output$tree_summary <- renderText({
    r <- selected_tree()
    if (is.null(r)) return("Select a locus")
    if (grepl("regional_tree_for_inspection",r$candidate_clade_rule[[1]] %||% "",fixed=TRUE)) {
      return(paste0("Full-region tree · ",r$candidate_lineage," · status=",r$tree_status,
                    "; ancestral outgroup=",ifelse(r$tree_has_ancestral_outgroup==1,"present","absent"),
                    ". Topology for inspection; not an automatic introgression call."))
    }
    paste0("status=", r$tree_status, "; INFO/AA ancestral outgroup=", ifelse(r$tree_has_ancestral_outgroup==1,"present","absent"),
           "; clade pass=", r$candidate_clade_pass,
           "; rule=", r$candidate_clade_rule, "; candidate lineage=", r$candidate_lineage,
           "; candidate-clade bootstrap=", r$candidate_clade_bootstrap,
           "; clade tips=", r$candidate_clade_n_tips,
           " (modern=", r$candidate_clade_modern_tips, ", archaic=", r$candidate_clade_archaic_tips, ")",
           "\nMixed-lineage edge bootstrap=", if ("mixed_lineage_bootstrap" %in% names(r)) r$mixed_lineage_bootstrap else NA,
           " (does not establish a lineage-specific introgression call)",
           "\nAll-tree bootstrap min/median/max=", r$bootstrap_min, "/", r$bootstrap_median, "/", r$bootstrap_max,
           ". Candidate bootstrap is branch repeatability, not an introgression probability.")
  })
  panel_b_bundle <- reactive({
    r <- selected_tree(); req(!is.null(r))
    validate(need(!is.na(r$tree_newick[[1]]) && nzchar(r$tree_newick[[1]]), "No phylogeny was constructed for this locus; see the evidence reason."))
    gu_b_bundle(r,dirname(phyml_file("evidence_trees.tsv")),r$locus_id[[1]])
  })
  output$phyml_tree <- renderPlot({
    gu_draw_panel_b(panel_b_bundle(),as.integer(input$tree_display_min_copies %||% "11"))
  },res=130)
  output$download_panel_b_pdf <- downloadHandler(
    filename=function()paste0(selected_locus()$locus_id[[1]],".phylogeny.panelB.pdf"),
    content=function(file) {
      grDevices::pdf(file,width=10,height=10.4,useDingbats=FALSE)
      on.exit(grDevices::dev.off())
      gu_draw_panel_b(panel_b_bundle(),as.integer(input$tree_display_min_copies %||% "11"))
    })
  browser_target <- reactive({
    r <- linked_region()
    if(!is.null(r)) return(data.frame(r,genome_build=current_build(),label="Selected region"))
    z <- browser_override()
    if (!is.null(z)) return(z)
    r <- selected_locus(); if (is.null(r)) return(NULL)
    data.frame(chr = r$chr[[1]], start = r$selected_start[[1]], end = r$selected_end[[1]],
               genome_build = r$genome_build[[1]], label = paste0("PhyML · ", r$locus_id[[1]]), stringsAsFactors = FALSE)
  })
  matching_region <- reactive({
    r<-linked_region();if(!is.null(r))return(r)
    r<-selected_locus();if(is.null(r))return(NULL)
    data.frame(chr=r$chr[[1]],start=r$selected_start[[1]],end=r$selected_end[[1]])
  })
  select_matching_record <- function(r,method) {
    selected_matching_record(r)
    browser_override(data.frame(chr=r$chr[[1]],start=r$start[[1]],end=r$end[[1]],genome_build=current_build(),label=paste0(method," · ",r$super_population[[1]] %||% "ALL")))
    linked_region(r[,c("chr","start","end"),drop=FALSE])
    session$sendCustomMessage("gu-region-selected",list(chr=as.character(r$chr[[1]]),start=r$start[[1]],end=r$end[[1]]))
  }
  open_matching_record <- function(r,method) {
    select_matching_record(r,method)
    updateSelectInput(session,paste0(method,"-chr"),selected=as.character(r$chr[[1]]))
    updateNumericInput(session,paste0(method,"-start"),value=r$start[[1]])
    updateNumericInput(session,paste0(method,"-end"),value=r$end[[1]])
    if("super_population" %in% names(r))updateSelectInput(session,paste0(method,"-population"),selected=r$super_population[[1]])
    bslib::nav_select("nav",method,session=session)
  }
  for(method in c("trace","as3"))gu_matching_server(paste0("overview_",method),method,Q,current_dataset,current_build,matching_region,select_matching_record,open_matching_record)
  output$genome_browser <- renderUI({
    r <- browser_target()
    if (is.null(r)) return(tags$div(class = "alert alert-info", "No locus or segment is available for this build."))
    u <- .gu_browser_urls(r, if(is.null(linked_region())) as.integer(input$browser_flank %||% 250000L) else 0L)
    session$sendCustomMessage("gu-igv-region",list(genome=u$genome,locus=u$locus,reference=.gu_local_references[[u$genome]]))
    tags$div(
      tags$div(class = "gu-browser-toolbar",
               tags$span(class = "gu-locus-label", paste0(u$label, " · ", u$build, " · ", u$locus)),
               tags$a("Open IGV ↗", href = u$igv, target = "_blank", rel = "noopener", class = "btn btn-sm btn-primary"),
               tags$a("Open UCSC ↗", href = u$ucsc, target = "_blank", rel = "noopener", class = "btn btn-sm btn-outline-primary"))
    )
  })
  observeEvent(list(current_build(), current_dataset()), {browser_override(NULL);linked_region(NULL);selected_matching_record(NULL)}, ignoreInit = TRUE)

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
