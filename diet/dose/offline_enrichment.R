# Local-only GMT over-representation analysis. No package or network dependency.
offline_gmt_enrichment <- function(genes, path) {
    if (!file.exists(path)) stop("Missing local GMT library: ", path)
    rows <- strsplit(readLines(path, warn = FALSE), "\t", fixed = TRUE)
    if (!length(rows) || any(lengths(rows) < 3L)) stop("Invalid/empty GMT library: ", path)
    terms <- vapply(rows, `[[`, character(1), 1L)
    if (anyDuplicated(terms) || any(!nzchar(terms))) stop("Duplicate/empty GMT terms: ", path)
    clean <- function(x) sort(unique(toupper(trimws(x[nzchar(trimws(x))]))))
    sets <- lapply(rows, function(x) clean(x[-c(1L, 2L)]))
    if (any(lengths(sets) == 0L)) stop("Empty GMT gene set: ", path)
    universe <- sort(unique(unlist(sets, use.names = FALSE)))
    genes <- clean(genes)
    query <- intersect(genes, universe)
    hits <- lapply(sets, intersect, y = query)
    k <- lengths(hits)
    m <- lengths(sets)
    n <- length(query)
    N <- length(universe)
    p <- if (n) phyper(k - 1L, m, N - m, n, lower.tail = FALSE) else rep(1, length(k))
    data.frame(
        Term = terms, Overlap = paste0(k, "/", m), P.value = p,
        Adjusted.P.value = p.adjust(p, method = "BH"),
        Genes = vapply(hits, paste, character(1), collapse = ";"),
        Background = "All annotated genes in this local GMT library",
        Background_n = N, Query_n = n, Query_unmapped_n = length(setdiff(genes, universe)),
        Library_file = normalizePath(path), Library_md5 = unname(tools::md5sum(path)),
        stringsAsFactors = FALSE
    )
}
