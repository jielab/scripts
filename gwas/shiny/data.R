# GWAS and BED cache functions. GWAS positions are 1-based; block BED
# edges remain 0-based half-open. No source GWAS or published BED is rewritten.
bplot_chr <- function(x) {
  x <- toupper(sub("^CHR", "", toupper(as.character(x))))
  x[x == "23"] <- "X"
  x
}

bplot_lengths <- function(build) {
  values <- if (build == "37") c(249250621,243199373,198022430,191154276,180915260,
    171115067,159138663,146364022,141213431,135534747,135006516,133851895,
    115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,
    48129895,51304566,155270560) else c(248956422,242193529,198295559,190214555,
    181538259,170805979,159345973,145138636,138394717,133797422,135086622,
    133275309,114364328,107043718,101991189,90338345,83257441,80373285,
    58617616,64444167,46709983,50818468,156040895)
  stats::setNames(values, c(as.character(1:22), "X"))
}

bplot_identity <- function(path) {
  info <- file.info(path)
  if (is.na(info$size)) stop("Missing file: ", path)
  list(path = normalizePath(path), size = info$size, mtime = as.numeric(info$mtime))
}

bplot_chain <- function(opt, build) file.path(opt$chain_dir,
  if (opt$source_build == "37" && build == "38") "hg19ToHg38.over.chain.gz" else "hg38ToHg19.over.chain.gz")

bplot_signature <- function(opt, i, build) {
  list(version = 1L, source = bplot_identity(opt$tracks$file[i]),
       source_build = opt$source_build, target_build = build, max_points = opt$max_points,
       chain = if (build != opt$source_build) bplot_identity(bplot_chain(opt, build)) else NULL,
       helper = unname(tools::md5sum(file.path(opt$app_dir, "data.R"))),
       liftover = if (build != opt$source_build) unname(tools::md5sum(file.path(opt$app_dir, "..", "f", "gwas_liftover.py"))) else NULL)
}

bplot_cache_dir <- function(opt, i, build) file.path(opt$output_dir, "cache", paste0("grch", build), sprintf("track_%02d", i))

bplot_cache_meta <- function(opt, i, build) {
  directory <- bplot_cache_dir(opt, i, build)
  if (!file.exists(file.path(directory, "meta.rds"))) return(NULL)
  tryCatch({
    meta <- readRDS(file.path(directory, "meta.rds"))
    if (!identical(meta$signature, bplot_signature(opt, i, build))) return(NULL)
    sizes <- file.info(file.path(directory, meta$cache_files))$size
    if (anyNA(sizes) || !identical(unname(sizes), meta$cache_sizes)) return(NULL)
    meta
  }, error = function(e) NULL)
}

bplot_atomic_rds <- function(object, path) {
  temp <- paste0(path, ".", Sys.getpid(), ".tmp")
  on.exit(unlink(temp), add = TRUE)
  saveRDS(object, temp, compress = FALSE)
  if (!file.rename(temp, path)) stop("Cannot publish cache: ", path)
}

bplot_thin <- function(data, budget) {
  if (nrow(data) <= budget) return(data)
  # Retain the strongest association in every horizontal bin, supplemented by
  # evenly sampled background points. A zoom queries the full cached data again.
  bins <- max(1L, floor(budget / 2))
  span <- max(data$POS) - min(data$POS) + 1
  bin <- floor((data$POS - min(data$POS)) / span * bins)
  peaks <- data.table::data.table(bin = bin, logp = data$LOGP)[,
    .(row = .I[which.max(logp)]), by = bin]$row
  background <- as.integer(seq(1, nrow(data), length.out = budget - length(peaks)))
  data[sort(unique(c(peaks, background)))]
}

bplot_read_gwas <- function(path, build) {
  command <- paste(if (nzchar(Sys.which("pigz"))) "pigz -dc -p 1 --" else "gzip -cd --", shQuote(path))
  compressed <- grepl("\\.(gz|bgz)$", path)
  header <- if (compressed) names(data.table::fread(cmd = command, nrows = 0L)) else names(data.table::fread(path, nrows = 0L))
  if (!all(c("CHR", "POS", "P") %in% header)) stop("GWAS requires CHR, POS and P: ", path)
  keep <- intersect(c("SNP", "CHR", "POS", "P", "LOG10P"), header)
  data <- if (compressed) data.table::fread(cmd = command, select = keep, showProgress = FALSE) else
    data.table::fread(path, select = keep, showProgress = FALSE)
  total <- nrow(data)
  data[, CHR := bplot_chr(CHR)]
  data[, `:=`(POS = suppressWarnings(as.numeric(POS)), P = suppressWarnings(as.numeric(P)))]
  lengths <- bplot_lengths(build)
  data <- data[CHR %in% names(lengths) & is.finite(POS) & POS >= 1 & POS == floor(POS) &
                 POS <= unname(lengths[CHR]) & is.finite(P) & P >= 0 & P <= 1]
  if (!"SNP" %in% names(data)) data[, SNP := paste0("chr", CHR, ":", POS)]
  data[, SNP := as.character(SNP)]
  data[, LOGP := -log10(pmax(P, 1e-300))]
  # Standardized LOG10P retains precision when a tiny P underflows to zero.
  if ("LOG10P" %in% names(data)) {
    data[, LOG10P := suppressWarnings(as.numeric(LOG10P))]
    data[P <= 1e-300 & is.finite(LOG10P) & LOG10P >= 300, LOGP := LOG10P]
    data[, LOG10P := NULL]
  }
  data.table::setkey(data, CHR, POS)
  list(data = data, input = total, excluded = total - nrow(data))
}

bplot_prepare <- function(opt, build, progress_file) {
  progress <- function(stage, track = 0L) {
    temp <- paste0(progress_file, ".tmp")
    jsonlite::write_json(list(stage = stage, track = track, total = nrow(opt$tracks)), temp, auto_unbox = TRUE)
    file.rename(temp, progress_file)
    cat(stage, "\n"); flush.console()
  }
  for (i in seq_len(nrow(opt$tracks))) {
    name <- opt$tracks$trait[i]
    progress(paste0(i, "/", nrow(opt$tracks), " · ", name, " · GRCh", build), i)
    if (!is.null(bplot_cache_meta(opt, i, build))) next
    directory <- bplot_cache_dir(opt, i, build)
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    signature <- bplot_signature(opt, i, build)
    path <- opt$tracks$file[i]
    if (build != opt$source_build) {
      progress(paste0(i, "/", nrow(opt$tracks), " · liftOver ", name, " → GRCh", build), i)
      path <- file.path(directory, paste0(name, ".gz"))
      command <- c(file.path(opt$app_dir, "..", "f", "gwas_liftover.py"), "--input", opt$tracks$file[i],
        "--output", path, "--qc-prefix", file.path(directory, "qc", name),
        "--source-build", opt$source_build, "--target-build", build,
        "--chain", bplot_chain(opt, build), "--liftOver", opt$liftover_bin)
      rc <- system2(Sys.which("python3"), shQuote(command))
      if (rc != 0L) stop("liftOver failed for ", name, "; see preparation log")
    }
    parsed <- bplot_read_gwas(path, build)
    data <- parsed$data
    if (!nrow(data)) stop("No valid variants on 1-22/X: ", name)
    overview <- list()
    counts <- numeric()
    for (chr in names(bplot_lengths(build))) {
      part <- data[.(chr)]
      part <- part[!is.na(POS)]
      if (!nrow(part)) next
      counts[chr] <- nrow(part)
      bplot_atomic_rds(part, file.path(directory, paste0("chr", chr, ".rds")))
      overview[[chr]] <- bplot_thin(part, max(2L, floor(opt$max_points / 23)))
    }
    bplot_atomic_rds(data.table::rbindlist(overview), file.path(directory, "overview.rds"))
    files <- c("overview.rds", paste0("chr", names(counts), ".rds"))
    meta <- list(signature = signature, counts = counts, total = nrow(data), excluded = parsed$excluded,
                 cache_files = files, cache_sizes = unname(file.info(file.path(directory, files))$size))
    bplot_atomic_rds(meta, file.path(directory, "meta.rds"))
    rm(data, parsed, part, overview); gc(verbose = FALSE)
  }
  progress("ready", nrow(opt$tracks))
}

# IDs are stable, one-based genomic row numbers within each race/build BED.
# BED [START, END) contains a 1-based GWAS POS when START < POS <= END.
block_cache <- new.env(parent = emptyenv())
bplot_blocks <- function(directory, race, build, chr) {
  path <- file.path(directory, paste0(race, ".", build, ".bed"))
  empty <- data.table::data.table(CHR = character(), START = numeric(), END = numeric(), ID = integer())
  if (!file.exists(path)) return(list(edges = numeric(), data = empty, note = paste0(race, "：无对应 block 文件")))
  identity <- bplot_identity(path)
  cached <- block_cache[[path]]
  if (is.null(cached) || !identical(cached$identity, identity)) {
    raw <- data.table::fread(path, header = FALSE, fill = TRUE, colClasses = "character", showProgress = FALSE)
    if (ncol(raw) < 3L) stop("Invalid block BED: ", path)
    data <- data.table::data.table(CHR = bplot_chr(trimws(raw[[1]])),
      START = suppressWarnings(as.numeric(trimws(raw[[2]]))), END = suppressWarnings(as.numeric(trimws(raw[[3]]))))
    data <- data[CHR %in% names(bplot_lengths(build)) & is.finite(START) & is.finite(END) & START >= 0 & END > START]
    data[, rank := match(CHR, names(bplot_lengths(build)))]
    data.table::setorder(data, rank, START, END)
    data[, `:=`(rank = NULL, ID = seq_len(.N))]
    block_cache[[path]] <- list(identity = identity, data = data)
  } else data <- cached$data
  note <- if (build == "37" && race %in% c("EAS", "SAS")) "ASN 来源" else paste0(race, " blocks")
  if (chr != "All") data <- data[CHR == chr]
  if (!nrow(data)) note <- paste0(note, " · chr", chr, " 无边界记录")
  list(edges = if (chr == "All") numeric() else sort(unique(c(data$START, data$END))), data = data, note = note)
}

bplot_block_ids <- function(data, blocks) {
  ids <- rep(NA_integer_, nrow(data))
  for (chr in unique(data$CHR)) {
    b <- blocks[CHR == chr]
    if (!nrow(b)) next
    rows <- which(data$CHR == chr)
    j <- findInterval(data$POS[rows] - 1, b$START)
    valid <- j > 0L
    valid[valid] <- data$POS[rows[valid]] <= b$END[j[valid]]
    ids[rows[valid]] <- b$ID[j[valid]]
  }
  ids
}

bplot_triplet <- function(directory, race, build, chr, id) {
  blocks <- bplot_blocks(directory, race, build, chr)$data
  j <- match(as.integer(id), blocks$ID)
  if (is.na(j)) return(NULL)
  adjacent <- blocks[seq(max(1L, j - 1L), min(nrow(blocks), j + 1L))]
  list(race = race, id = as.integer(id), chr = chr, build = build,
       start = min(adjacent$START) + 1, end = max(adjacent$END),
       selected = c(blocks$START[j] + 1, blocks$END[j]), ids = adjacent$ID)
}

bplot_view <- function(opt, build, chr, range, show_blocks = TRUE) {
  lengths <- bplot_lengths(build)
  offsets <- stats::setNames(c(0, head(cumsum(lengths), -1)), names(lengths))
  lapply(seq_len(nrow(opt$tracks)), function(i) {
    directory <- bplot_cache_dir(opt, i, build)
    meta <- readRDS(file.path(directory, "meta.rds"))
    path <- file.path(directory, if (chr == "All") "overview.rds" else paste0("chr", chr, ".rds"))
    available <- chr == "All" || chr %in% names(meta$counts)
    data <- if (available) readRDS(path) else data.table::data.table(SNP = character(), CHR = character(), POS = numeric(), P = numeric(), LOGP = numeric())
    total <- if (chr == "All") meta$total else nrow(data)
    if (chr != "All" && nrow(data)) {
      data <- data[POS >= range[1] & POS <= range[2]]
      total <- nrow(data)
      data <- bplot_thin(data, opt$max_points)
    }
    data[, X := if (chr == "All") POS + unname(offsets[CHR]) else POS]
    blocks <- bplot_blocks(opt$block_dir, opt$tracks$race[i], build, chr)
    data[, ID := bplot_block_ids(data, blocks$data)]
    data[, BLOCK := ifelse(is.na(ID), "", paste("block:", opt$tracks$race[i], ID))]
    if (!show_blocks) blocks$edges <- numeric()
    blocks$edges <- blocks$edges[blocks$edges >= range[1] & blocks$edges <= range[2]]
    list(data = data, total = total, blocks = blocks, track = opt$tracks[i, ])
  })
}
