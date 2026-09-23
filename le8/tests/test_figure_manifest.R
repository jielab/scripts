# Run from the repository root: Rscript tests/test_figure_manifest.R
source("f/figure_policy.R")

local({
  root <- tempfile("figure-manifest-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  rawdir <- file.path(root, "output")
  incoming <- file.path(root, "incoming")
  dir.create(rawdir)
  dir.create(incoming)

  # PGS and the grouped figure renderer deliberately supply different metadata.
  old <- data.frame(file = "c3.Fig7.pgs_coloc_triangulation.png",
    group = "c3.Fig7.pgs_coloc_triangulation", panels = 3L,
    source = "PGS focus aggregate results")
  fresh <- data.frame(file = "c3.Fig1.posterior_evidence.png",
    group = "posterior_evidence", panels = 6L,
    sources = "c3.Fig1.evidence_triage.png", policy = "grouped panels")
  data.table::fwrite(old, file.path(rawdir, "figure_manifest.csv"))
  data.table::fwrite(fresh, file.path(incoming, "figure_manifest.csv"))
  audit <- data.frame(source = fresh$sources, group = fresh$group, status = "rendered")
  data.table::fwrite(audit, file.path(incoming, "figure_omission_audit.csv"))
  # Refresh only copies files; distinct bytes let us verify renumbering integrity.
  writeBin(charToRaw("existing PGS figure"), file.path(rawdir, old$file))
  writeBin(charToRaw("new coloc figure"), file.path(incoming, fresh$file))
  expected <- unname(tools::md5sum(c(file.path(rawdir, old$file),
    file.path(incoming, fresh$file))))

  merged <- le8_refresh_figure_files(rawdir, incoming)
  stopifnot(nrow(merged) == 2L,
    identical(merged$group, c(old$group, fresh$group)),
    identical(merged$source, c(old$source, NA_character_)),
    identical(merged$sources, c(NA_character_, fresh$sources)),
    identical(merged$policy, c(NA_character_, fresh$policy)),
    identical(unname(tools::md5sum(file.path(rawdir, merged$file))), expected),
    length(list.files(file.path(rawdir, "_history"))) >= 2L)

  # A repeated render replaces its theme, retaining the PGS figure and metadata.
  fresh <- fresh[, rev(names(fresh))]
  fresh$policy <- "updated policy"
  data.table::fwrite(fresh, file.path(incoming, "figure_manifest.csv"))
  repeated <- le8_refresh_figure_files(rawdir, incoming)
  stopifnot(nrow(repeated) == 2L,
    identical(repeated$file, merged$file),
    repeated$source[1] == old$source,
    repeated$policy[2] == "updated policy",
    identical(unname(tools::md5sum(file.path(rawdir, repeated$file))), expected))

  # Refresh without incoming figures and first-time rendering remain supported.
  unchanged <- le8_refresh_figure_files(rawdir)
  stopifnot(nrow(unchanged) == 2L, identical(unchanged$file, repeated$file))
  emptydir <- file.path(root, "empty")
  dir.create(emptydir)
  initial <- le8_refresh_figure_files(emptydir, incoming)
  stopifnot(nrow(initial) == 1L, initial$group == fresh$group,
    file.exists(file.path(emptydir, initial$file)))
})
cat("Figure manifest regression checks passed.\n")
