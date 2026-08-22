source("global.R", local = TRUE, encoding = "UTF-8")
if (!interactive()) sink(stderr(), type = "output")

# ---- Helpers ----
.build_choices <- function(paths, names_vec) {
  x <- paste(paths, names_vec, sep = "$")
  names(x) <- names_vec
  x
}

gwas_choices <- function(upload_df = NULL) {
  if (!is.null(upload_df) && nrow(upload_df)) {
    return(.build_choices(upload_df$datapath, upload_df$name))
  } else {
    if (!dir.exists(PATHS$GWAS_DIR)) return(character(0))
    gf <- list.files(PATHS$GWAS_DIR, full.names = TRUE)
    gf <- gf[grepl("\\.(vcf(\\.gz)?|txt(\\.gz)?|tsv(\\.gz)?|gz)$", gf, ignore.case = TRUE)]
    if (!length(gf)) return(character(0))
    return(.build_choices(gf, basename(gf)))
  }
}

.extract_path <- function(choice_vec) {
  # input like "/path/file$filename.txt" -> returns "/path/file"
  if (is.null(choice_vec) || !length(choice_vec)) return(character(0))
  sapply(strsplit(choice_vec, "\\$"), `[`, 1)
}

# Read GWAS sumstats robustly using your deal_vcf/deal_vcf_old if present
.read_sumstats <- function(fp) {
  stopifnot(length(fp) == 1)
  ext <- tolower(basename(fp))
  g <- NULL
  if (grepl("\\.vcf(\\.gz)?$", ext)) {
    if (exists("deal_vcf_old")) {
      g <- deal_vcf_old(fp)
    } else {
      g <- data.table::fread(fp, nThread = 1, data.table = FALSE)
    }
  } else {
    if (exists("deal_vcf")) {
      g <- deal_vcf(fp)
    } else {
      g <- data.table::fread(fp, sep = "\t", header = TRUE, nThread = 1, data.table = FALSE)
    }
  }
  # Normalize columns: need "snp" and "p"
  nm <- names(g)
  if (!"snp" %in% nm && "SNP" %in% nm) g$snp <- g$SNP
  if (!"p"   %in% nm && "P"   %in% nm) g$p   <- g$P
  if (!"snp" %in% names(g) || !"p" %in% names(g)) {
    stop("GWAS file missing required columns (snp/P or SNP/P). File: ", fp)
  }
  g <- g[, c("snp","p")]
  g$p <- suppressWarnings(as.numeric(g$p))
  g <- g[is.finite(g$p) & g$p > 0, , drop = FALSE]
  g
}

# Pairwise correlation on -log10(P) using intersecting SNPs (downsample if huge)
.pair_cor <- function(df1, df2, max_snps = 200000L) {
  snps <- intersect(df1$snp, df2$snp)
  n <- length(snps)
  if (n < 100L) return(NA_real_)
  if (n > max_snps) snps <- snps[sample.int(n, max_snps)]
  v1 <- -log10(df1$p[match(snps, df1$snp)])
  v2 <- -log10(df2$p[match(snps, df2$snp)])
  if (sum(is.finite(v1) & is.finite(v2)) < 100L) return(NA_real_)
  suppressWarnings(cor(v1, v2, use = "complete.obs", method = "pearson"))
}

# ---- Server ----
shinyServer(function(input, output, session) {

  # Correlation: keep matrix here and draw directly (no files)
  cor_mat <- reactiveVal(NULL)

  output$report_C2 <- renderUI({
    if (is.null(cor_mat())) {
      return(h4("Click 'Run' to compute and show the correlation heatmap."))
    }
    plotOutput("plotc2", height = 500)
  })

  output$plotc2 <- renderPlot({
    req(cor_mat())
    corrplot::corrplot(cor_mat(), method = "shade", tl.col = "black")
  })

  # File pickers using the shared helper
  output$files_C21 <- renderUI({
    files <- gwas_choices(input$gwasFilecs)
    if (!length(files)) return(helpText("No GWAS files found in data/gwas"))
    checkboxGroupInput("gwasfiles21", "GWAS to LDSC (X-axis)", choices = files, selected = files[1])
  })
  output$files_C22 <- renderUI({
    files <- gwas_choices(input$gwasFilecs)
    if (!length(files)) return(helpText("No GWAS files found in data/gwas"))
    sel <- if (length(files) > 1) files[2] else files[1]
    checkboxGroupInput("gwasfiles22", "GWAS to LDSC (Y-axis)", choices = files, selected = sel)
  })
  output$files_C31 <- renderUI({
    files <- gwas_choices(input$gwasFilecs)
    if (!length(files)) return(helpText("No GWAS files found in data/gwas"))
    radioButtons("gwasfiles31", "exposure GWAS file", choices = files, selected = files[1])
  })
  output$files_C32 <- renderUI({
    files <- gwas_choices(input$gwasFilecs)
    if (!length(files)) return(helpText("No GWAS files found in data/gwas"))
    sel <- if (length(files) > 1) files[2] else files[1]
    radioButtons("gwasfiles32", "outcome GWAS file", choices = files, selected = sel)
  })
  output$files_C41 <- renderUI({
    files <- gwas_choices(input$gwasFilecs)
    if (!length(files)) return(helpText("No GWAS files found in data/gwas"))
    checkboxGroupInput("gwasfiles4", "input GWAS files to Colocalization", choices = files, selected = files[1])
  })
  output$files_C6 <- renderUI({
    files <- gwas_choices(input$gwasFilecs)
    if (!length(files)) return(helpText("No GWAS files found in data/gwas"))
    checkboxGroupInput("gwasfiles6", "GWAS files to be integrated", choices = files, selected = files[1])
  })

  # ---- Lightweight, file-free correlation runner for the Correlation tab ----
  fun_runc2 <- function() {
    x_sel <- input$gwasfiles21
    y_sel <- input$gwasfiles22
    if (is.null(x_sel) || is.null(y_sel) || !length(x_sel) || !length(y_sel)) {
      showNotification("Please select at least one GWAS on each axis.", type = "warning")
      return(invisible(NULL))
    }

    x_paths <- .extract_path(x_sel)
    y_paths <- .extract_path(y_sel)
    labs_x  <- ifelse(grepl("\\$", x_sel), sapply(strsplit(x_sel, "\\$"), `[`, 2), basename(x_paths))
    labs_y  <- ifelse(grepl("\\$", y_sel), sapply(strsplit(y_sel, "\\$"), `[`, 2), basename(y_paths))

    withProgress(message = "Computing correlations", value = 0, {
      incProgress(0.05)

      # Read all selected GWAS once
      all_paths <- unique(c(x_paths, y_paths))
      g_list <- vector("list", length(all_paths))
      names(g_list) <- all_paths

      for (i in seq_along(all_paths)) {
        g_list[[i]] <- tryCatch(
          .read_sumstats(all_paths[i]),
          error = function(e) { showNotification(paste("Read failed:", basename(all_paths[i])), type="error"); NULL }
        )
        incProgress(0.4/length(all_paths))
      }
      # Filter out failures
      ok <- vapply(g_list, function(z) !is.null(z) && nrow(z) > 0, logical(1))
      if (!all(ok)) {
        g_list <- g_list[ok]
        x_keep <- x_paths %in% names(g_list)
        y_keep <- y_paths %in% names(g_list)
        x_paths <- x_paths[x_keep]; labs_x <- labs_x[x_keep]
        y_paths <- y_paths[y_keep]; labs_y <- labs_y[y_keep]
      }
      if (!length(x_paths) || !length(y_paths)) {
        showNotification("No valid GWAS after reading files.", type = "error")
        return(invisible(NULL))
      }

      # Build correlation matrix (X by Y)
      M <- matrix(NA_real_, nrow = length(labs_x), ncol = length(labs_y),
                  dimnames = list(labs_x, labs_y))
      for (i in seq_along(x_paths)) {
        for (j in seq_along(y_paths)) {
          M[i, j] <- .pair_cor(g_list[[x_paths[i]]], g_list[[y_paths[j]]])
        }
        incProgress(0.55/length(x_paths))
      }

      # If X and Y sets are identical, show a symmetric matrix
      if (setequal(labs_x, labs_y) && length(labs_x) == length(labs_y)) {
        # Mirror into a square matrix
        full <- matrix(NA_real_, nrow = length(labs_x), ncol = length(labs_x),
                       dimnames = list(labs_x, labs_x))
        full[rownames(M), colnames(M)] <- M
        full[colnames(M), rownames(M)] <- t(M)
        diag(full) <- 1
        cor_mat(full)
      } else {
        cor_mat(M)
      }

      output$report_C2 <- renderUI({ plotOutput("plotc2", height = 500) })
      incProgress(0.4)
    })
  }

  # ---- Other tabs: keep your existing functions if they exist ----
  observeEvent(input$runC2, { tryCatch(fun_runc2(), error = function(e) showNotification(e$message, type="error")) })
  observeEvent(input$runC3, { tryCatch(fun_runc3(), error = function(e) NULL) })
  observeEvent(input$runC4, { tryCatch(fun_runc4(), error = function(e) NULL) })
  observeEvent(input$runC5, { tryCatch(fun_runc5(), error = function(e) NULL) })
  observeEvent(input$runC6, { tryCatch(fun_runc6(), error = function(e) NULL) })
  observeEvent(input$analyze, {
    # Run correlation by default; others only if you've defined those helpers elsewhere
    tryCatch(fun_runc2(), error = function(e) NULL)
    for (f in list(fun_runc3, fun_runc4, fun_runc5, fun_runc6)) {
      if (exists(deparse(substitute(f)), mode = "function")) tryCatch(f(), error = function(e) NULL)
    }
  })
})
