# GU Shiny server
# v6 haplotype viewer. It reads loci_avcf
# matrices plus the already-produced hap_match.tsv and haplotype_sample_map.tsv.
# v6 shows true 1KG IDs, matched rows plus negative controls, scale-safe match
# percentages, and IGV/UCSC genomic context. normalize also exports browser files.

.gu_expected_arch <- c("Altai", "Chagyr", "Vindija", "Denisova", "Denisova25")

.gu_canon_arch <- function(x) {
  y <- tolower(gsub("[^a-z0-9]", "", as.character(x)))
  ifelse(grepl("altai", y), "Altai",
    ifelse(grepl("chag", y), "Chagyr",
      ifelse(grepl("vind", y), "Vindija",
        ifelse(grepl("denisova25|den25", y), "Denisova25",
          ifelse(grepl("denis", y), "Denisova", as.character(x))))))
}

.gu_arch_label <- function(f) .gu_canon_arch(tools::file_path_sans_ext(basename(f)))

.gu_modern_base <- function(gt, hap, ref, alt) {
  g <- sub(":.*$", "", as.character(gt))
  g <- gsub("\\|", "/", g)
  a <- if (hap == 1L) sub("/.*$", "", g) else ifelse(grepl("/", g), sub("^.*/", "", g), g)
  out <- rep("N", length(a))
  out[a == "0"] <- ref[a == "0"]
  out[a == "1"] <- alt[a == "1"]
  direct <- toupper(a) %in% c("A","C","G","T")
  out[direct] <- toupper(a[direct])
  out
}

.gu_arch_base <- function(gt, ref, alt) {
  g <- sub(":.*$", "", as.character(gt))
  g <- gsub("\\|", "/", g)
  sp <- strsplit(g, "/", fixed = TRUE)
  vapply(seq_along(sp), function(i) {
    z <- sp[[i]]
    z <- z[!is.na(z) & nzchar(z) & z != "."]
    if (!length(z)) return("N")
    zu <- unique(z)
    if (length(zu) != 1L) return("N")
    a <- toupper(zu[1])
    if (a == "0") return(ref[i])
    if (a == "1") return(alt[i])
    if (a %in% c("A","C","G","T")) return(a)
    "N"
  }, character(1))
}

.gu_split_sources <- function(x) {
  z <- unlist(strsplit(as.character(x %||% ""), ";", fixed=TRUE), use.names=FALSE)
  z <- trimws(z[nzchar(trimws(z))])
  unique(.gu_canon_arch(z))
}

.gu_first_dir <- function(x) {
  x <- unique(x[!is.na(x) & nzchar(x)])
  z <- x[dir.exists(x)]
  if (length(z)) z[1] else NA_character_
}

.gu_locus_paths <- function(r) {
  raw <- as.character(r$raw_file[1] %||% "")
  root0 <- if (nzchar(raw)) dirname(dirname(raw)) else ""
  ar <- Sys.getenv("GU_ANALYSIS_ROOT", "/mnt/d/analysis/gu")
  m <- as.character(r$method[1] %||% "loci_avcf")
  tr <- as.character(r$trait[1] %||% "")
  roots <- c(root0, file.path(ar, m, tr), file.path(ar, paste0("loci_", sub("^loci_", "", m)), tr))
  root <- .gu_first_dir(roots)
  if (is.na(root)) stop("Cannot resolve the locus analysis directory. raw_file=", raw)
  id <- as.character(r$locus_id[1])
  mats <- c(file.path(root,"mat",id), file.path(root,"mat",tr,id))
  matdir <- .gu_first_dir(mats)
  if (is.na(matdir)) stop("Missing locus matrix directory for ", id, ": tried ", paste(mats, collapse="; "))
  list(
    root=root, id=id, matdir=matdir,
    kg=file.path(matdir,"kg.tsv"),
    samples=file.path(matdir,"kg.samples.tsv"),
    hap_report=file.path(root,"report","hap_match.tsv"),
    sample_report=file.path(root,"report","haplotype_sample_map.tsv"),
    core=file.path(root,"hap",paste0(id,".core.tsv"))
  )
}

.gu_read_locus_base <- function(r) {
  if (!grepl("avcf", as.character(r$method[1]), ignore.case=TRUE))
    stop("The sequence-level haplotype inspector currently uses loci_avcf output. Select a loci_avcf locus.")
  p <- .gu_locus_paths(r)
  if (!file.exists(p$kg)) stop("Missing kg.tsv: ", p$kg)
  kg <- data.table::fread(p$kg, header=FALSE, showProgress=FALSE)
  if (ncol(kg) < 6L) stop("Malformed kg.tsv (need chr,pos,ref,alt,aa + modern samples): ", p$kg)
  data.table::setnames(kg, c("chr","pos","ref","alt","aa",paste0("s",seq_len(ncol(kg)-5L))))
  kg[, `:=`(chr=as.character(chr), pos=as.integer(pos), ref=toupper(as.character(ref)), alt=toupper(as.character(alt)))]
  kg <- kg[nchar(ref)==1L & nchar(alt)==1L & ref %chin% c("A","C","G","T") & alt %chin% c("A","C","G","T")]
  if (!nrow(kg)) stop("No biallelic SNVs in locus matrix: ", p$kg)

  # Keep the normalized inherited interval when possible; fall back to the complete locus matrix.
  st <- suppressWarnings(as.numeric(r$start[1])); en <- suppressWarnings(as.numeric(r$end[1]))
  if (is.finite(st) && is.finite(en) && en > st) {
    z <- kg[pos >= st & pos <= en]
    if (nrow(z) >= 2L) kg <- z
  }

  samp <- setdiff(names(kg), c("chr","pos","ref","alt","aa"))
  sm <- data.table(sidx=samp, sample=samp)
  if (file.exists(p$samples)) {
    zz <- data.table::fread(p$samples, header=FALSE, showProgress=FALSE)
    if (nrow(zz)) {
      sm <- data.table(sidx=paste0("s",seq_len(nrow(zz))), sample=as.character(zz[[1]]))
      sm <- sm[sidx %in% samp]
    }
  }

  # Five high-coverage archaic panels. Missing files remain N so the panel is explicit.
  amat <- matrix("N", nrow=length(.gu_expected_arch), ncol=nrow(kg),
                 dimnames=list(.gu_expected_arch, as.character(kg$pos)))
  af <- list.files(p$matdir, pattern="\\.tsv$", full.names=TRUE)
  af <- af[!basename(af) %in% c("kg.tsv","kg.samples.tsv")]
  available_arch <- character()
  for (f in af) {
    lab <- .gu_arch_label(f)
    if (!lab %in% .gu_expected_arch) next
    a <- tryCatch(data.table::fread(f, header=FALSE, showProgress=FALSE), error=function(e) NULL)
    if (is.null(a) || ncol(a) < 5L || !nrow(a)) next
    a <- a[,1:5]
    data.table::setnames(a,c("chr","pos","ref","alt","gt"))
    a[, `:=`(pos=as.integer(pos), ref=toupper(as.character(ref)), alt=toupper(as.character(alt)))]
    a[, base := .gu_arch_base(gt, ref, alt)]
    ii <- match(kg$pos, a$pos)
    ok <- !is.na(ii)
    amat[lab,ok] <- a$base[ii[ok]]
    available_arch <- unique(c(available_arch,lab))
  }

  # Reproduce the key modern-polymorphism filter used in loci_avcf without constructing
  # all 2*N modern haplotype sequences.  Only 0/1 biallelic GTs are needed here.
  G <- as.matrix(kg[, ..samp])
  gv <- sub(":.*$", "", as.vector(G))
  gv <- gsub("\\|", "/", gv)
  G <- matrix(gv, nrow=nrow(kg), ncol=length(samp))
  m00 <- matrix(G %in% c("0/0"), nrow=nrow(G)); m11 <- matrix(G %in% c("1/1"), nrow=nrow(G))
  mhet <- matrix(G %in% c("0/1","1/0"), nrow=nrow(G))
  mh0 <- matrix(G %in% c("0"), nrow=nrow(G)); mh1 <- matrix(G %in% c("1"), nrow=nrow(G))
  n0 <- 2*rowSums(m00) + rowSums(mhet) + rowSums(mh0)
  n1 <- 2*rowSums(m11) + rowSums(mhet) + rowSums(mh1)
  # Display every 1KG-polymorphic biallelic SNP in the inherited interval.
  # Keep the stricter >=2 copies definition separately for the original AVCF
  # matching-site logic used by loci_avcf.
  modern_poly <- n0 > 0L & n1 > 0L
  avcf_poly <- n0 >= 2L & n1 >= 2L

  src <- .gu_split_sources(r$source[1])
  src <- intersect(src, .gu_expected_arch)
  if (!length(src)) src <- intersect(.gu_expected_arch, available_arch)
  arch_called <- if (length(src)) apply(amat[src,,drop=FALSE],2,function(z) all(z %in% c("A","C","G","T"))) else rep(TRUE,nrow(kg))
  matching <- avcf_poly & arch_called

  diagnostic <- rep(FALSE,nrow(kg))
  if (file.exists(p$core)) {
    co <- tryCatch(data.table::fread(p$core, showProgress=FALSE), error=function(e) NULL)
    if (!is.null(co) && nrow(co)) {
      pc <- intersect(c("pos","POS","bp","BP"), names(co))[1]
      dc <- intersect(c("is_diagnostic_archaic","diagnostic"), names(co))[1]
      if (!is.na(pc) && !is.na(dc)) {
        flag <- as.character(co[[dc]]) %in% c("TRUE","T","1","true")
        diagnostic <- kg$pos %in% suppressWarnings(as.integer(co[[pc]][flag]))
      }
    }
  }

  # loci_avcf has already calculated haplotype-to-archaic match statistics.
  # Read those values rather than recomputing them in Shiny.  The separate
  # haplotype_sample_map.tsv converts internal s1418_2-style copy IDs back to
  # the real 1000 Genomes sample and phased haplotype.
  hr <- data.table()
  if (file.exists(p$hap_report)) {
    hr <- tryCatch(data.table::fread(p$hap_report, showProgress=FALSE), error=function(e) data.table())
    if (nrow(hr)) {
      if ("trait" %in% names(hr) && !is.na(r$trait[1]) && nzchar(as.character(r$trait[1])))
        hr <- hr[as.character(trait)==as.character(r$trait[1])]
      if ("id" %in% names(hr)) hr <- hr[as.character(id)==as.character(r$locus_id[1])]
    }
  }

  cm <- data.table()
  if (nrow(hr) && "copies" %in% names(hr)) {
    match_cols <- names(hr)[grepl("_match$", names(hr), ignore.case=TRUE)]
    cm <- data.table::rbindlist(lapply(seq_len(nrow(hr)), function(i) {
      cp <- unlist(strsplit(as.character(hr$copies[i] %||% ""), ";", fixed=TRUE), use.names=FALSE)
      cp <- cp[nzchar(cp)]
      if (!length(cp)) return(NULL)
      z <- data.table(
        copy=cp,
        hap_id=if ("hap_id" %in% names(hr)) as.character(hr$hap_id[i]) else NA_character_,
        hap_n=if ("n" %in% names(hr)) suppressWarnings(as.integer(hr$n[i])) else NA_integer_,
        best_lineage=if ("best_lineage" %in% names(hr)) as.character(hr$best_lineage[i]) else NA_character_,
        matched_archaics=if ("matched_archaics" %in% names(hr)) as.character(hr$matched_archaics[i]) else NA_character_,
        best_arch=if ("best_arch" %in% names(hr)) as.character(hr$best_arch[i]) else NA_character_,
        best_match=if ("best_match" %in% names(hr)) suppressWarnings(as.numeric(hr$best_match[i])) else NA_real_,
        carry_risk=if ("carry_risk" %in% names(hr)) as.character(hr$carry_risk[i]) else NA_character_
      )
      for (an in .gu_expected_arch) {
        hit <- match_cols[vapply(match_cols, function(cc) {
          identical(.gu_canon_arch(sub("_match$", "", cc, ignore.case=TRUE)), an)
        }, logical(1))]
        z[[paste0("match_",an)]] <- if (length(hit)) suppressWarnings(as.numeric(hr[[hit[1]]][i])) else NA_real_
      }
      z
    }), fill=TRUE)
  }

  # Preferred mapping: report/haplotype_sample_map.tsv, generated by loci_avcf.
  # It is authoritative for the s-index -> real 1KG sample mapping.
  sr <- data.table()
  if (file.exists(p$sample_report)) {
    sr <- tryCatch(data.table::fread(p$sample_report, showProgress=FALSE), error=function(e) data.table())
    if (nrow(sr)) {
      if ("trait" %in% names(sr) && !is.na(r$trait[1]) && nzchar(as.character(r$trait[1])))
        sr <- sr[as.character(trait)==as.character(r$trait[1])]
      if ("id" %in% names(sr)) sr <- sr[as.character(id)==as.character(r$locus_id[1])]
      if ("sample" %in% names(sr)) data.table::setnames(sr,"sample","sample_id")
      if ("haplotype" %in% names(sr)) data.table::setnames(sr,"haplotype","sample_haplotype")
      keep <- intersect(c("copy","sample_id","sample_haplotype","pop","super_pop","sidx"), names(sr))
      if ("copy" %in% keep) sr <- unique(sr[, ..keep]) else sr <- data.table()
    }
  }

  if (nrow(cm)) {
    if (nrow(sr)) cm <- merge(cm, sr, by="copy", all.x=TRUE, sort=FALSE)
    if (!"sample_id" %in% names(cm)) cm[, sample_id := NA_character_]
    if (!"sample_haplotype" %in% names(cm)) cm[, sample_haplotype := NA_character_]
    if (!"pop" %in% names(cm)) cm[, pop := NA_character_]
    if (!"super_pop" %in% names(cm)) cm[, super_pop := NA_character_]
    cm[, sidx := sub("_([12])$", "", copy)]
    fallback_sample <- sm$sample[match(cm$sidx, sm$sidx)]
    bad <- is.na(cm$sample_id) | !nzchar(as.character(cm$sample_id))
    cm$sample_id[bad] <- fallback_sample[bad]
    bad2 <- is.na(cm$sample_id) | !nzchar(as.character(cm$sample_id))
    cm$sample_id[bad2] <- cm$sidx[bad2]
    fallback_hap <- suppressWarnings(as.integer(sub("^.*_([12])$", "\\1", cm$copy)))
    badh <- is.na(cm$sample_haplotype) | !nzchar(as.character(cm$sample_haplotype))
    cm$sample_haplotype[badh] <- fallback_hap[badh]
    cm[, best_arch_canon := .gu_canon_arch(best_arch)]
    # best_arch is haplotype-specific; matched_archaics is locus-level and therefore
    # must not be used to label every haplotype as a source match.
    cm[, source_hit := if (length(src)) best_arch_canon %in% src else TRUE]
    cm <- cm[!duplicated(copy)]
    cm <- cm[order(-source_hit, -best_match, -hap_n, na.last=TRUE)]
  }

  source_candidates <- if (nrow(cm) && length(src)) cm[source_hit == TRUE] else cm
  if (!nrow(source_candidates) && nrow(cm)) source_candidates <- cm

  list(row=r, paths=p, kg=kg, samples=sm, samp=samp, archaic=amat,
       available_arch=available_arch, source_arch=src,
       polymorphic=modern_poly, matching=matching, diagnostic=diagnostic,
       copy_meta=cm, source_candidates=source_candidates)
}

.gu_copy_matrix <- function(b, copies, idx=seq_len(nrow(b$kg))) {
  if (!length(copies)) return(matrix(character(),nrow=0,ncol=length(idx)))
  out <- matrix("N",nrow=length(copies),ncol=length(idx),dimnames=list(copies,as.character(b$kg$pos[idx])))
  for (i in seq_along(copies)) {
    m <- regexec("^(s[0-9]+)_([12])$", copies[i])
    z <- regmatches(copies[i],m)[[1]]
    if (length(z) != 3L || !z[2] %in% b$samp) next
    h <- as.integer(z[3])
    out[i,] <- .gu_modern_base(b$kg[[z[2]]][idx], h, b$kg$ref[idx], b$kg$alt[idx])
  }
  out
}

.gu_cap_idx <- function(idx,maxn,priority=integer()) {
  idx <- sort(unique(as.integer(idx)))
  maxn <- max(1L,as.integer(maxn))
  if (length(idx) <= maxn) return(idx)
  pr <- intersect(idx,unique(as.integer(priority)))
  if (length(pr) >= maxn) return(sort(pr[seq_len(maxn)]))
  rem <- setdiff(idx,pr); need <- maxn-length(pr)
  take <- unique(pmax(1L,pmin(length(rem),round(seq(1,length(rem),length.out=need)))))
  sort(unique(c(pr,rem[take])))[seq_len(min(maxn,length(unique(c(pr,rem[take])))))]
}

.gu_copy_label <- function(b, cp) {
  vapply(cp, function(one) {
    z <- if (nrow(b$copy_meta)) b$copy_meta[copy == one][1] else data.table()
    sidx <- sub("_([12])$", "", one)
    hap0 <- sub("^.*_([12])$", "\\1", one)
    sample0 <- if (nrow(z) && "sample_id" %in% names(z)) as.character(z$sample_id[1]) else NA_character_
    if (is.na(sample0) || !nzchar(sample0)) {
      sample0 <- b$samples$sample[match(sidx,b$samples$sidx)]
      if (is.na(sample0) || !nzchar(sample0)) sample0 <- paste0("UNMAPPED(",sidx,")")
    }
    if (nrow(z) && "sample_haplotype" %in% names(z) && !is.na(z$sample_haplotype[1]) && nzchar(as.character(z$sample_haplotype[1])))
      hap0 <- as.character(z$sample_haplotype[1])
    paste0(sample0," | hap",hap0)
  }, character(1))
}

.gu_pct <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (!length(x) || !is.finite(x[1])) return(NA_real_)
  z <- x[1]
  if (z < 0) return(NA_real_)
  if (z <= 1.000001) z <- 100*z
  # Reports exist in fraction, percent, and percent*100 form.
  while (z > 100.000001) z <- z/100
  min(100, z)
}

.gu_hover_one <- function(b, cp) {
  lab <- .gu_copy_label(b, cp)
  z <- if (nrow(b$copy_meta)) b$copy_meta[copy == cp][1] else data.table()
  if (!nrow(z)) return(lab)
  parts <- c(lab)
  if ("pop" %in% names(z) && !is.na(z$pop[1]) && nzchar(as.character(z$pop[1]))) {
    pp <- as.character(z$pop[1]); sp <- if ("super_pop" %in% names(z)) as.character(z$super_pop[1]) else ""
    parts <- c(parts, paste0("Population: ",pp, if (!is.na(sp) && nzchar(sp)) paste0(" / ",sp) else ""))
  }
  if ("hap_id" %in% names(z) && !is.na(z$hap_id[1])) {
    nn <- if ("hap_n" %in% names(z) && is.finite(suppressWarnings(as.numeric(z$hap_n[1])))) paste0("; n=",z$hap_n[1]) else ""
    parts <- c(parts,paste0("Haplotype class: ",z$hap_id[1],nn))
  }
  bm <- if ("best_match" %in% names(z)) .gu_pct(z$best_match[1]) else NA_real_
  ba <- if ("best_arch" %in% names(z)) as.character(z$best_arch[1]) else ""
  if (is.finite(bm)) parts <- c(parts,sprintf("Best archaic match: %s (%.1f%%)",ba,bm))
  mtxt <- character()
  for (an in .gu_expected_arch) {
    cc <- paste0("match_",an)
    if (cc %in% names(z)) {
      vv <- .gu_pct(z[[cc]][1])
      if (is.finite(vv)) mtxt <- c(mtxt,sprintf("%s %.1f%%",an,vv))
    }
  }
  if (length(mtxt)) parts <- c(parts,paste0("Matches: ",paste(mtxt,collapse="; ")))
  if ("carry_risk" %in% names(z) && !is.na(z$carry_risk[1])) parts <- c(parts,paste0("Carry risk allele: ",z$carry_risk[1]))
  parts <- c(parts,paste0("Reported locus source: ",as.character(b$row$source[1])))
  paste(parts,collapse="\n")
}

.gu_pick_matched <- function(b, n_each) {
  allcp <- c(paste0(b$samp,"_1"), paste0(b$samp,"_2"))
  x <- data.table::copy(b$source_candidates)
  if (!nrow(x)) return(character())
  x <- x[copy %in% allcp]
  if (!nrow(x)) return(character())
  x[, hap_key := ifelse(is.na(hap_id) | !nzchar(hap_id), copy, hap_id)]
  rf <- if ("carry_risk" %in% names(x)) tolower(as.character(x$carry_risk)) %in% c("true","t","1") else rep(FALSE,nrow(x))
  x[, risk_flag := rf]
  data.table::setorder(x, -risk_flag, -best_match, -hap_n, na.last=TRUE)
  # Prefer one representative copy from each distinct reported haplotype class,
  # with risk-carrying source-matching haplotypes first when they exist.
  rep1 <- x[!duplicated(hap_key), copy]
  out <- head(rep1, n_each)
  if (length(out) < n_each) {
    fill <- setdiff(x$copy, out)
    out <- c(out, head(fill, n_each-length(out)))
  }
  unique(out)
}

.gu_pick_controls <- function(b, matched, n=10L) {
  allcp <- c(paste0(b$samp,"_1"), paste0(b$samp,"_2"))
  reported <- if (nrow(b$copy_meta) && "copy" %in% names(b$copy_meta)) as.character(b$copy_meta$copy) else character()
  pool <- setdiff(allcp, unique(c(matched, reported)))
  if (!length(pool)) pool <- setdiff(allcp, matched)
  n <- min(max(0L,as.integer(n)),length(pool))
  if (!n) return(character())
  key <- paste0(as.character(b$row$locus_id[1]),"|",as.character(b$row$chr[1]),"|",as.character(b$row$start[1]))
  chars <- utf8ToInt(key); seed <- sum(chars * seq_along(chars)) %% .Machine$integer.max
  old <- if (exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)) get(".Random.seed",envir=.GlobalEnv) else NULL
  on.exit(if (is.null(old)) { if (exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)) rm(".Random.seed",envir=.GlobalEnv) } else assign(".Random.seed",old,envir=.GlobalEnv),add=TRUE)
  set.seed(seed); sample(pool,n,replace=FALSE)
}

.gu_browser_urls <- function(r, flank=250000L) {
  chr <- paste0("chr",sub("^chr","",as.character(r$chr[1]),ignore.case=TRUE))
  st <- max(1L,as.integer(r$start[1])-as.integer(flank)); en <- as.integer(r$end[1])+as.integer(flank)
  locus <- paste0(chr,":",st,"-",en)
  build <- if (grepl("38",Sys.getenv("GU_BUILD","37"))) "hg38" else "hg19"
  list(locus=locus, build=build,
       igv=paste0("https://igv.org/app/?genome=",build,"&locus=",utils::URLencode(locus,reserved=TRUE)),
       ucsc=paste0("https://genome.ucsc.edu/cgi-bin/hgTracks?db=",build,"&position=",utils::URLencode(locus,reserved=TRUE)))
}

function(input, output, session) {
  con <- dbConnect(SQLite(), db_path)
  session$onSessionEnded(function() dbDisconnect(con))
  methods <- q(con,"SELECT DISTINCT method FROM segments ORDER BY method")$method
  sources <- q(con,"SELECT DISTINCT source_class FROM segments ORDER BY source_class")$source_class
  updateSelectInput(session,"reg_method",choices=c("ALL",methods)); updateSelectInput(session,"data_method",choices=c("ALL",methods))

  metric <- function(name,fallback) {
    if (DBI::dbExistsTable(con,"overview_metrics")) scalar_q(con,"SELECT value FROM overview_metrics WHERE metric=?",list(name)) else scalar_q(con,fallback)
  }
  output$n_samples <- renderText(comma(metric("n_samples","SELECT COUNT(DISTINCT sample_id) FROM sample_burden")))
  output$n_segments <- renderText(comma(metric("n_segments","SELECT SUM(n_segments) FROM sample_burden")))
  output$n_codes <- renderText(comma(metric("n_codes","SELECT COUNT(*) FROM segment_catalog")))
  output$n_methods <- renderText(comma(metric("n_methods","SELECT COUNT(*) FROM (SELECT method FROM segments UNION SELECT method FROM loci)")))
  output$db_info <- renderText(paste("SQLite:",db_path,"\nSize:",format(file.info(db_path)$size,big.mark=",",scientific=FALSE),"bytes"))

  has_table <- function(x) DBI::dbExistsTable(con,x)
  active_locus_idx <- reactiveVal(1L)
  overview_loci <- reactive({
    q(con,"SELECT rowid AS locus_rowid,trait,locus_id,chr,start,end,source,status,n_carriers,raw_file,'loci_avcf' AS method FROM loci WHERE method='loci_avcf' ORDER BY CASE WHEN chr GLOB '[0-9]*' THEN CAST(chr AS INT) ELSE 99 END,start,locus_id")
  })
  output$overview_loci_detail <- renderDT({
    d <- overview_loci()
    show <- d[,c("trait","locus_id","chr","start","end","source","status","n_carriers"),drop=FALSE]
    sel <- active_locus_idx()
    if (!nrow(show) || !is.finite(sel) || sel < 1L || sel > nrow(show)) sel <- 1L
    datatable(
      show,
      selection=list(mode="single", selected=sel, target="row"),
      options=list(pageLength=8, scrollX=TRUE, dom='tip'),
      rownames=FALSE,
      callback=DT::JS(
        "table.on('click', 'tbody tr', function() {",
        "  var idx = table.row(this).index();",
        "  if (idx !== undefined && idx !== null) {",
        "    Shiny.setInputValue('hap_locus_click', idx + 1, {priority: 'event'});",
        "  }",
        "});",
        "table.on('dblclick', 'tbody tr', function() { var idx=table.row(this).index(); Shiny.setInputValue('locus_row_dblclick',idx+1,{priority:'event'}); });"
      )
    )
  })

  observeEvent(input$hap_locus_click, {
    d <- overview_loci()
    i <- suppressWarnings(as.integer(input$hap_locus_click))
    if (length(i) && is.finite(i) && i >= 1L && i <= nrow(d)) active_locus_idx(i)
  }, ignoreInit=TRUE)
  observeEvent(input$locus_row_dblclick, {
    i<-as.integer(input$locus_row_dblclick); d<-overview_loci()
    if (is.finite(i) && i>=1L && i<=nrow(d)) { active_locus_idx(i); bslib::nav_select("gu_nav","locus_viewer",session=session) }
  },ignoreInit=TRUE)
  observeEvent(input$locus_back,bslib::nav_select("gu_nav","overview",session=session),ignoreInit=TRUE)

  # Retain DT's native selection event as a second path (keyboard / accessibility).
  observeEvent(input$overview_loci_detail_rows_selected, {
    d <- overview_loci()
    i <- input$overview_loci_detail_rows_selected
    if (length(i) && i >= 1L && i <= nrow(d) && !identical(i, active_locus_idx())) active_locus_idx(i)
  }, ignoreInit=TRUE)

  selected_locus <- reactive({
    d <- overview_loci(); if (!nrow(d)) return(NULL)
    i <- active_locus_idx()
    if (!length(i) || !is.finite(i) || i < 1L || i > nrow(d)) i <- 1L
    d[i,,drop=FALSE]
  })

  active_bin <- reactiveVal(NULL)
  browser_target <- reactiveVal("locus")
  output$individual_section_title <- renderText({
    n <- metric("n_samples","SELECT COUNT(DISTINCT sample_id) FROM segments WHERE sample_id IS NOT NULL")
    paste0("Individual based results (N = ",comma(n),")")
  })
  observe({
    m<-input$overview_ind_method %||% "ibdmix"
    ss<-q(con,"SELECT DISTINCT source FROM segment_density_1mb WHERE method=? ORDER BY source",list(m))$source
    updateSelectInput(session,"overview_ind_source",choices=ss,selected=if(length(ss)) ss[1] else character())
  })
  density_data <- reactive({
    req(input$overview_ind_method,input$overview_ind_source)
    q(con,"SELECT method,source,chr,bin_start,n_calls,n_carriers,total_bp FROM segment_density_1mb WHERE method=? AND source=? ORDER BY CASE WHEN chr GLOB '[0-9]*' THEN CAST(chr AS INT) ELSE 99 END,bin_start",list(input$overview_ind_method,input$overview_ind_source))
  })
  top_bins <- reactive({
    d<-density_data(); if(!nrow(d)) return(d)
    d<-head(d[order(-d$n_carriers,-d$n_calls),,drop=FALSE],25)
    avail<-q(con,"SELECT DISTINCT method FROM segment_density_1mb")$method
    for(m in c("loci_avcf","ibdmix","trace","as3")) {
      nm<-paste0(m,"_N")
      if(m=="loci_avcf") {
        d[[nm]]<-vapply(seq_len(nrow(d)),function(i) {
          z<-q(con,"SELECT MAX(n_carriers) n FROM loci WHERE method='loci_avcf' AND chr=? AND end>? AND start<?",list(as.character(d$chr[i]),as.numeric(d$bin_start[i]),as.numeric(d$bin_start[i])+1000000))
          if(!nrow(z)||is.na(z$n[1])) NA_integer_ else as.integer(z$n[1])
        },integer(1)); next
      }
      if(!m %in% avail) { d[[nm]]<-NA_integer_; next }
      d[[nm]]<-vapply(seq_len(nrow(d)),function(i) {
        z<-q(con,"SELECT COUNT(DISTINCT sample_id) n FROM segments WHERE method=? AND chr=? AND end>? AND start<?",list(m,as.character(d$chr[i]),as.numeric(d$bin_start[i]),as.numeric(d$bin_start[i])+1000000))
        as.integer(z$n[1])
      },integer(1))
    }
    d$region<-paste0("chr",d$chr,":",comma(d$bin_start+1),"-",comma(d$bin_start+1000000))
    d
  })
  output$overview_individual_bins <- renderDT({
    d<-top_bins(); show<-d[,c("region","loci_avcf_N","ibdmix_N","trace_N","as3_N","n_calls","total_bp"),drop=FALSE]
    datatable(show,rownames=FALSE,selection="single",options=list(dom="tip",pageLength=10,scrollX=TRUE),callback=DT::JS(
      "table.on('click','tbody tr',function(){var i=table.row(this).index();Shiny.setInputValue('ind_bin_click',i+1,{priority:'event'});});",
      "table.on('dblclick','tbody tr',function(){var i=table.row(this).index();Shiny.setInputValue('ind_bin_dblclick',i+1,{priority:'event'});});"))
  })
  observeEvent(input$ind_bin_click,{
    i<-as.integer(input$ind_bin_click); d<-top_bins()
    if(is.finite(i)&&i>=1L&&i<=nrow(d)){active_bin(d[i,,drop=FALSE]);browser_target("bin")}
  },ignoreInit=TRUE)
  observeEvent(input$ind_bin_dblclick,{
    i<-as.integer(input$ind_bin_dblclick); d<-top_bins()
    if(is.finite(i)&&i>=1L&&i<=nrow(d)){active_bin(d[i,,drop=FALSE]);browser_target("bin");bslib::nav_select("gu_nav","individual_viewer",session=session)}
  },ignoreInit=TRUE)
  observeEvent(input$hap_locus_click,browser_target("locus"),ignoreInit=TRUE)
  observeEvent(input$overview_ind_method,{active_bin(NULL)},ignoreInit=TRUE)
  observeEvent(input$overview_ind_source,{active_bin(NULL)},ignoreInit=TRUE)
  observeEvent(input$individual_back,bslib::nav_select("gu_nav","overview",session=session),ignoreInit=TRUE)

  output$overview_individual_density <- renderPlotly({
    d<-density_data(); validate(need(nrow(d),"No normalized data for this method/source.")); d$xMb<-d$bin_start/1e6
    b<-active_bin(); if(!is.null(b)) d<-d[d$chr==b$chr[1],,drop=FALSE]
    p<-plot_ly(d,x=~xMb,y=~n_carriers,color=~factor(chr,levels=chr_order),type="scatter",mode="lines",source="density",text=~paste0("chr",chr,":",xMb,"-",xMb+1," Mb<br>carriers=",comma(n_carriers),"<br>calls=",comma(n_calls)),hoverinfo="text")
    if(!is.null(b)) p<-add_markers(p,data=b,x=~(bin_start/1e6+.5),y=~n_carriers,inherit=FALSE,marker=list(color="#e31a1c",size=11),name="selected")
    p %>% layout(xaxis=list(title=if(is.null(b)) "Position within chromosome (Mb)" else paste0("chr",b$chr[1]," position (Mb)")),yaxis=list(title="Carriers per 1 Mb bin"),legend=list(orientation="h",title=list(text="Chromosome")))
  })

  locus_base <- reactive({
    r <- selected_locus(); req(!is.null(r))
    tryCatch(.gu_read_locus_base(r), error=function(e) structure(list(error=conditionMessage(e),row=r),class="gu_locus_error"))
  })

  locus_view <- reactive({
    b <- locus_base()
    if (inherits(b,"gu_locus_error")) return(b)
    n_each <- max(1L,as.integer(input$hap_n_each %||% 7L))

    # Use the loci_avcf-reported source-matching haplotypes. The internal copy ID is
    # retained only for extracting phased alleles from kg.tsv; labels come from the
    # real 1KG sample mapping in haplotype_sample_map.tsv.
    matched <- .gu_pick_matched(b, n_each)
    if (!length(matched)) return(structure(list(error="No source-matching modern haplotypes were found in hap_match.tsv for this locus.",row=b$row),class="gu_locus_error"))

    idx_all <- which(b$polymorphic)
    if (!length(idx_all)) return(structure(list(error="No polymorphic biallelic 1KG SNPs were found in the selected inherited interval.",row=b$row),class="gu_locus_error"))
    priority <- intersect(idx_all,which(b$diagnostic))
    idx <- .gu_cap_idx(idx_all,as.integer(input$hap_max_sites %||% 150L),priority)

    controls <- .gu_pick_controls(b,matched,10L)
    shown_copies <- c(matched,controls)
    M <- .gu_copy_matrix(b,shown_copies,idx)
    A <- b$archaic[.gu_expected_arch,idx,drop=FALSE]
    # Omit references that have no callable base at any displayed site.  A modern
    # MATCH row without an actual reference row is misleading, so fail clearly if
    # the selected interval has no callable archaic sequence at all.
    # `%in%` drops matrix dimensions in R; restore them before rowSums so this
    # also works for loci with only one callable reference/site.
    called_arch <- matrix(A %in% c("A","C","G","T"), nrow=nrow(A), ncol=ncol(A),
                          dimnames=dimnames(A))
    keep_arch <- rowSums(called_arch) > 0L
    A <- A[keep_arch,,drop=FALSE]
    arch_names <- rownames(A)
    if (!nrow(A)) return(structure(list(error="No callable archaic A/C/G/T reference sequence is available for the displayed sites; matched modern rows are therefore hidden.",row=b$row),class="gu_locus_error"))
    rownames(M) <- .gu_copy_label(b,shown_copies)
    rownames(A) <- paste0("[ARCH] ",arch_names)
    B <- rbind(A,M)

    # Technical table and hover summaries use loci_avcf's precomputed *_match values.
    mm <- if (nrow(b$copy_meta)) b$copy_meta[match(matched,copy)] else data.table()
    sim <- data.table(label=.gu_copy_label(b,matched), copy_id=matched)
    if (nrow(mm)) {
      add <- intersect(c("sample_id","sample_haplotype","pop","super_pop","hap_id","hap_n","best_arch","best_match","carry_risk",paste0("match_",.gu_expected_arch)),names(mm))
      sim <- cbind(sim,mm[,..add])
      for (cc in c("best_match",paste0("match_",.gu_expected_arch))) if (cc %in% names(sim))
        sim[[cc]] <- vapply(sim[[cc]], function(x) { z=.gu_pct(x); if (is.finite(z)) round(z,1) else NA_real_ }, numeric(1))
    }
    hover <- setNames(vapply(matched,function(cp) .gu_hover_one(b,cp),character(1)),matched)

    # Pick the panel-wide best reference from the precomputed match percentages.
    # It becomes the green baseline; the remaining reference identities receive
    # stable, distinct colours.
    arch_scores <- setNames(rep(-Inf,length(arch_names)),arch_names)
    if (nrow(mm)) for (an in arch_names) {
      cc<-paste0("match_",an)
      if (cc %in% names(mm)) { z<-vapply(mm[[cc]],.gu_pct,numeric(1)); if (any(is.finite(z))) arch_scores[an]<-mean(z[is.finite(z)]) }
    }
    best_arch <- if (any(is.finite(arch_scores))) names(which.max(arch_scores)) else arch_names[1]

    control_hover <- setNames(paste0(.gu_copy_label(b,controls),"\nRandom non-matched haplotype (negative control)"),controls)
    list(base=b,idx=idx,bases=B,matched=matched,controls=controls,shown_copies=shown_copies,
         similarity=sim,hover=c(hover,control_hover),
         arch_names=arch_names,best_arch=best_arch,match_meta=mm,
         polymorphic_n=sum(b$polymorphic),matching_n=sum(b$matching),display_n=length(idx))
  })

  output$haplotype_title <- renderUI({
    r <- selected_locus()
    if (is.null(r)) return(tags$span(class="text-muted","Select a locus_evidence row on the Overview tab."))
    tags$div(
      tags$h5(sprintf("%s | %s | chr%s:%s-%s",r$trait[1],r$locus_id[1],r$chr[1],format(r$start[1],big.mark=","),format(r$end[1],big.mark=","))),
      tags$p(class="mb-2",sprintf("Reported archaic source: %s",r$source[1]))
    )
  })

  output$genome_browser <- renderUI({
    if (identical(browser_target(),"bin") && !is.null(active_bin())) {
      b<-active_bin(); r<-data.frame(chr=b$chr[1],start=b$bin_start[1]+1,end=b$bin_start[1]+1000000)
    } else r<-selected_locus()
    if (is.null(r)) return(tags$div(class="alert alert-info","Select a locus or an individual-based 1 Mb bin."))
    u <- .gu_browser_urls(r,as.integer(input$browser_flank %||% 250000L))
    tags$div(class="gu-browser-wrap",
      tags$div(class="gu-browser-toolbar",tags$b(paste0(u$build," · ",u$locus)),
        tags$a("Open IGV ↗",href=u$igv,target="_blank",class="btn btn-sm btn-primary"),
        tags$a("Open UCSC ↗",href=u$ucsc,target="_blank",class="btn btn-sm btn-outline-primary")),
      tags$iframe(src=u$igv,class="gu-igv-frame",title="IGV-Web genome browser"))
  })

  individual_region_data <- reactive({
    b<-active_bin(); req(!is.null(b))
    q(con,"SELECT sample_id,method,source_class,chr,start,end,length_bp,haplotype,score,posterior,segment_code FROM segments WHERE chr=? AND end>? AND start<? ORDER BY method,start,end LIMIT 20000",list(as.character(b$chr[1]),as.numeric(b$bin_start[1]),as.numeric(b$bin_start[1])+1000000))
  })
  individual_region_summary <- reactive({
    b<-active_bin(); req(!is.null(b))
    q(con,"SELECT method,source_class,COUNT(*) n_calls,COUNT(DISTINCT sample_id) n_carriers,SUM(length_bp) total_bp,AVG(length_bp) mean_length_bp,AVG(score) mean_score,AVG(posterior) mean_posterior FROM segments WHERE chr=? AND end>? AND start<? GROUP BY method,source_class ORDER BY method,source_class",list(as.character(b$chr[1]),as.numeric(b$bin_start[1]),as.numeric(b$bin_start[1])+1000000))
  })
  output$individual_region_title <- renderUI({b<-active_bin();req(!is.null(b));tags$h5(sprintf("%s / %s — chr%s:%s-%s",b$method[1],b$source[1],b$chr[1],comma(b$bin_start[1]+1),comma(b$bin_start[1]+1000000)))})
  output$individual_region_carriers <- renderText({b<-active_bin();req(!is.null(b));comma(b$n_carriers[1])})
  output$individual_region_calls <- renderText({b<-active_bin();req(!is.null(b));comma(b$n_calls[1])})
  output$individual_region_bp <- renderText({b<-active_bin();req(!is.null(b));fmt_bp(b$total_bp[1])})
  output$individual_region_methods <- renderText({d<-individual_region_summary();length(unique(d$method))})
  output$individual_region_comparison <- renderDT({datatable(individual_region_summary(),rownames=FALSE,options=list(dom="tip",pageLength=10,scrollX=TRUE))})
  output$individual_region_segments <- renderDT({datatable(individual_region_data(),rownames=FALSE,options=list(pageLength=12,scrollX=TRUE))})
  output$individual_region_plot <- renderPlotly({
    d<-individual_region_data();validate(need(nrow(d),"No overlapping segments."));d$track<-paste(d$method,d$source_class,sep=" / ");d$width<-(d$end-d$start)/1e6
    p<-plot_ly();for(tr in unique(d$track)){z<-d[d$track==tr,,drop=FALSE];p<-add_bars(p,data=z,y=~track,x=~width,base=~(start/1e6),orientation="h",name=tr,text=~paste0(sample_id,"<br>",method," / ",source_class,"<br>chr",chr,":",comma(start),"-",comma(end),"<br>",segment_code),hoverinfo="text")}
    p %>% layout(barmode="overlay",xaxis=list(title="Genomic position (Mb)"),yaxis=list(title=""),showlegend=FALSE)
  })

  output$haplotype_matrix <- renderUI({
    v <- locus_view()
    if (inherits(v,"gu_locus_error")) return(tags$div(class="alert alert-warning",v$error))
    b <- v$base; idx <- v$idx; B <- v$bases
    pos <- b$kg$pos[idx]; chr0 <- b$kg$chr[idx][1]
    n_arch <- length(v$arch_names)
    palette <- setNames(c("#16883f","#d62728","#1769d2","#e67e22","#795548")[seq_len(n_arch)],
                        c(v$best_arch,setdiff(v$arch_names,v$best_arch)))
    best_seq <- B[match(paste0("[ARCH] ",v$best_arch),rownames(B)),]

    css <- tags$style(HTML("
      .gu-hap-scroll{overflow-x:auto;max-width:100%;padding-bottom:8px}
      .gu-hap-table{border-collapse:collapse;font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;font-size:13px;line-height:1}
      .gu-hap-table th.gu-rowlab{position:sticky;left:0;z-index:2;background:white;text-align:right;white-space:nowrap;padding:5px 10px 5px 4px;border-right:1px solid #999;font-family:var(--bs-body-font-family);font-weight:500}
      .gu-hap-table td.gu-base{min-width:19px;width:19px;height:22px;text-align:center;vertical-align:middle;padding:0;border-right:1px solid #eee;cursor:default}
      .gu-hap-table td.gu-base.gu-coloured{font-weight:700}
      .gu-hap-table tr.gu-arch-last th,.gu-hap-table tr.gu-arch-last td{border-bottom:2px solid #555}
      .gu-hap-table tr.gu-control-first th,.gu-hap-table tr.gu-control-first td{border-top:2px solid #b42318}
      .gu-hap-table tr.gu-control th.gu-rowlab{color:#c62828;font-weight:700}
      .gu-hap-ruler{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;font-size:12px;margin:0 0 6px 0}
    "))

    row_tags <- lapply(seq_len(nrow(B)), function(i) {
      is_arch <- i <= n_arch
      cp <- if (!is_arch) v$shown_copies[i-n_arch] else NA_character_
      is_control <- !is_arch && cp %in% v$controls
      an <- if (is_arch) v$arch_names[i] else NA_character_
      row_info <- if (is_arch) paste0("Archaic reference: ",an,if (an==v$best_arch) " (panel-wide best match)" else "") else v$hover[[cp]]
      cells <- lapply(seq_along(idx), function(j) {
        aa <- B[i,j]; shown <- if (aa %in% c("A","C","G","T")) aa else ""
        colour <- "#1f1f1f"
        matched_to <- "none"
        if (nzchar(shown) && is_arch) {
          matched_to <- an
          colour <- if (an==v$best_arch || (best_seq[j] %in% c("A","C","G","T") && aa==best_seq[j])) "#16883f" else palette[[an]]
        } else if (nzchar(shown)) {
          mi <- match(cp,v$matched)
          meta <- if (!is.na(mi) && nrow(v$match_meta)) v$match_meta[mi,] else NULL
          pref <- character()
          if (!is.null(meta)) {
            ba <- .gu_canon_arch(meta$best_arch[1] %||% "")
            vals <- setNames(vapply(v$arch_names,function(x) { cc<-paste0("match_",x); if (cc %in% names(meta)) .gu_pct(meta[[cc]][1]) else NA_real_ },numeric(1)),v$arch_names)
            pref <- unique(c(ba,names(sort(vals,decreasing=TRUE,na.last=NA)),v$arch_names))
          } else pref <- v$arch_names
          pref <- pref[pref %in% v$arch_names]
          for (candidate in pref) {
            rr<-match(paste0("[ARCH] ",candidate),rownames(B))
            if (!is.na(rr) && B[rr,j] %in% c("A","C","G","T") && aa==B[rr,j]) { matched_to<-candidate; colour<-palette[[candidate]]; break }
          }
        }
        site_info <- paste0("Site: chr",chr0,":",format(pos[j],big.mark=","),
                            "\nAllele: ",if (nzchar(shown)) shown else "missing/ambiguous",
                            "\nColour match: ",matched_to,
                            "\nREF/ALT: ",b$kg$ref[idx[j]],"/",b$kg$alt[idx[j]])
        tags$td(class=if(nzchar(shown)) "gu-base gu-coloured" else "gu-base",style=paste0("color:",colour),title=paste(row_info,site_info,sep="\n"),shown)
      })
      cls <- if (i == n_arch) "gu-arch-last" else if (is_control && identical(cp,v$controls[1])) "gu-control gu-control-first" else if (is_control) "gu-control" else NULL
      tags$tr(class=cls, tags$th(class="gu-rowlab",title=row_info,rownames(B)[i]), cells)
    })

    tags$div(
      css,
      tags$p(class="gu-hap-ruler",sprintf("chr%s: %s → %s | %d polymorphic SNPs shown (%d in interval). Hover any modern row/base for the precomputed archaic-match statistics.",
        chr0,format(min(pos),big.mark=","),format(max(pos),big.mark=","),length(pos),v$polymorphic_n)),
      tags$div(class="gu-hap-scroll",tags$table(class="gu-hap-table",tags$tbody(row_tags)))
    )
  })

  output$haplotype_similarity <- renderDT({
    v <- locus_view()
    if (inherits(v,"gu_locus_error")) return(datatable(data.frame(message=v$error),options=list(dom='t'),rownames=FALSE))
    datatable(v$similarity,rownames=FALSE,options=list(pageLength=20,scrollX=TRUE,dom='tip'))
  })

  output$haplotype_note <- renderText({
    v <- locus_view()
    if (inherits(v,"gu_locus_error")) return(v$error)
    b <- v$base
    missing_arch <- setdiff(.gu_expected_arch,b$available_arch)
    paste0(
      "Matrix: ",b$paths$matdir,"\n",
      "Sample mapping: ",b$paths$sample_report,"\n",
      "Precomputed haplotype matches: ",b$paths$hap_report,"\n",
      "Polymorphic SNPs in interval: ",v$polymorphic_n,"; displayed: ",v$display_n,"; original AVCF matching sites: ",v$matching_n,".\n",
      "Modern rows: ",length(v$matched)," source-matching + ",length(v$controls)," random non-matched 1KG phased haplotypes (negative controls).\n",
      "Green baseline reference: ",v$best_arch,"; matching alleles inherit their archaic-reference colour; alleles matching none of the displayed references are black.\n",
      if (length(missing_arch)) paste0("Archaic matrix files not available/callable: ",paste(missing_arch,collapse=", ")," (cells left blank).\n") else "",
      "Only A/C/G/T are printed. Missing or ambiguous archaic/modern calls are left blank. Match values come from loci_avcf hap_match.tsv and are scale-normalized to 0-100% for display."
    )
  })

  ind_query <- eventReactive(input$ind_go, {
    req(nzchar(input$ind_sample)); methods0=input$ind_methods; wh="sample_id=? AND chr=?"; pars=list(input$ind_sample,input$ind_chr)
    if(!isTRUE(input$ind_allchr)){wh=paste0(wh," AND end>? AND start<?");pars=c(pars,input$ind_start,input$ind_end)}
    if(length(methods0)){qs=paste(rep('?',length(methods0)),collapse=',');wh=paste0(wh," AND method IN (",qs,")");pars=c(pars,as.list(methods0))}
    q(con,paste0("SELECT * FROM segments WHERE ",wh," ORDER BY start,end"),pars)
  },ignoreInit=TRUE)
  output$individual_plot <- renderPlotly({
    d<-ind_query(); validate(need(nrow(d),"No segments found")); d$track<-paste(d$method,ifelse(is.na(d$haplotype)|d$haplotype=='','person',paste0('hap',d$haplotype)),sep=' / '); d$lenMb=(d$end-d$start)/1e6
    p<-plot_ly()
    for(tr in unique(d$track)){z=d[d$track==tr,]; p<-add_bars(p,data=z,y=~track,x=~lenMb,base=~(start/1e6),orientation='h',name=tr,text=~paste0(method,' / ',source_class,'<br>',chr,':',start,'-',end,'<br>',segment_code,'<br>score=',score,' posterior=',posterior),hoverinfo='text')}
    p %>% layout(barmode='overlay',xaxis=list(title=paste0('chr',input$ind_chr,' position (Mb)')),yaxis=list(title=''),showlegend=FALSE)
  })
  output$individual_codes <- renderDT({ req(input$ind_go); datatable(q(con,"SELECT segment_code,source_class,dosage,n_methods,methods_support,max_score,max_posterior FROM carriers WHERE sample_id=? ORDER BY source_class,segment_code",list(input$ind_sample)),options=list(pageLength=12,scrollX=TRUE)) })
  output$individual_burden <- renderDT({ req(input$ind_go); datatable(q(con,"SELECT * FROM sample_burden WHERE sample_id=? ORDER BY method,source_class",list(input$ind_sample)),options=list(dom='t',scrollX=TRUE)) })

  reg_query <- eventReactive(input$reg_go,{
    wh="chr=? AND end>? AND start<?"; pars=list(input$reg_chr,input$reg_start,input$reg_end)
    if(input$reg_method!='ALL'){wh=paste0(wh,' AND method=?');pars=c(pars,input$reg_method)}
    q(con,paste0("SELECT sample_id,method,source_class,start,end,score,posterior,segment_code FROM segments WHERE ",wh," ORDER BY start LIMIT ",as.integer(input$reg_limit)),pars)
  },ignoreInit=FALSE)
  output$region_plot <- renderPlotly({
    d=reg_query(); validate(need(nrow(d),"No segment calls")); d$mid=(d$start+d$end)/2/1e6; d$length=(d$end-d$start)/1e3
    plot_ly(d,x=~mid,y=~method,color=~source_class,size=~pmax(1,sqrt(length)),type='scatter',mode='markers',text=~paste0(sample_id,'<br>',source_class,'<br>',start,'-',end,'<br>',segment_code),hoverinfo='text') %>% layout(xaxis=list(title=paste0('chr',input$reg_chr,' position (Mb)')),yaxis=list(title='Method'))
  })
  output$region_loci <- renderDT({ d=q(con,"SELECT * FROM loci WHERE chr=? AND end>? AND start<? ORDER BY start",list(input$reg_chr,input$reg_start,input$reg_end)); datatable(d,options=list(pageLength=10,scrollX=TRUE)) })

  locus_query <- eventReactive(input$loc_go,{
    wh=c('1=1'); pars=list()
    if(input$loc_method!='ALL'){wh=c(wh,'method=?');pars=c(pars,input$loc_method)}
    if(nzchar(input$loc_trait)){wh=c(wh,'trait LIKE ?');pars=c(pars,paste0('%',input$loc_trait,'%'))}
    if(input$loc_chr!='ALL'){wh=c(wh,'chr=?');pars=c(pars,input$loc_chr)}
    q(con,paste0('SELECT * FROM loci WHERE ',paste(wh,collapse=' AND '),' ORDER BY chr,start LIMIT 10000'),pars)
  },ignoreInit=FALSE)
  output$locus_plot <- renderPlotly({
    d=locus_query(); validate(need(nrow(d),'No locus evidence for this filter')); d$mid=(d$start+d$end)/2/1e6; d$width=pmax(0.001,(d$end-d$start)/1e6)
    plot_ly(d,x=~mid,y=~method,color=~method,size=~pmax(4,sqrt(width*100)),type='scatter',mode='markers',text=~paste0(method,'<br>',trait,' / ',locus_id,'<br>chr',chr,':',start,'-',end,'<br>source=',source,'<br>status=',status),hoverinfo='text') %>% layout(xaxis=list(title='Position within chromosome (Mb)'),yaxis=list(title=''))
  })
  output$locus_table <- renderDT({ datatable(locus_query(),options=list(pageLength=20,scrollX=TRUE)) })

  catalog_query <- reactive({
    wh=c('1=1'); pars=list()
    if(input$code_source!='ALL'){wh=c(wh,'c.source_class=?');pars=c(pars,input$code_source)}
    if(input$code_chr!='ALL'){wh=c(wh,'c.chr=?');pars=c(pars,input$code_chr)}
    sql=paste0("SELECT c.*,COALESCE(x.n_carriers,0) n_carriers FROM segment_catalog c LEFT JOIN (SELECT segment_code,COUNT(DISTINCT sample_id) n_carriers FROM carriers GROUP BY segment_code) x USING(segment_code) WHERE ",paste(wh,collapse=' AND ')," AND COALESCE(x.n_carriers,0)>=? ORDER BY n_carriers DESC LIMIT 5000")
    q(con,sql,c(pars,list(input$code_min_carriers)))
  })
  output$catalog_table <- renderDT(datatable(catalog_query(),options=list(pageLength=20,scrollX=TRUE)))
  output$code_carriers <- renderDT({ if(!nzchar(input$code_id)) return(datatable(data.frame())); datatable(q(con,"SELECT * FROM carriers WHERE segment_code=? ORDER BY sample_id LIMIT 10000",list(input$code_id)),options=list(pageLength=15,scrollX=TRUE)) })

  cc_data <- eventReactive(input$cc_go,{ req(nzchar(input$cc_sample)); q(con,"SELECT method,start,end FROM segments WHERE sample_id=? AND chr=? AND end>? AND start<? ORDER BY method,start",list(input$cc_sample,input$cc_chr,input$cc_start,input$cc_end)) },ignoreInit=TRUE)
  cc_matrix <- reactive({
    d=as.data.table(cc_data()); if(!nrow(d)) return(data.frame()); ms=sort(unique(d$method)); M=matrix(NA_real_,length(ms),length(ms),dimnames=list(ms,ms))
    for(i in seq_along(ms)) for(j in seq_along(ms)){a=d[method==ms[i],.(start,end)];b=d[method==ms[j],.(start,end)];inter=intersection_bp(a,b); uni=interval_bp(a)+interval_bp(b)-inter;M[i,j]=if(uni>0) inter/uni else NA}
    M
  })
  output$concordance_plot <- renderPlotly({ M=cc_matrix(); validate(need(length(M),"No calls")); plot_ly(x=colnames(M),y=rownames(M),z=M,type='heatmap',zmin=0,zmax=1) %>% layout(xaxis=list(title=''),yaxis=list(title='')) })
  output$concordance_table <- renderDT({ M=cc_matrix(); if(!length(M)) return(datatable(data.frame())); datatable(cbind(method=rownames(M),as.data.frame(round(M,3))),options=list(dom='t',scrollX=TRUE)) })

  data_query <- reactive({
    wh=c('1=1');pars=list(); if(nzchar(input$data_sample)){wh=c(wh,'sample_id LIKE ?');pars=c(pars,paste0('%',input$data_sample,'%'))}; if(input$data_method!='ALL'){wh=c(wh,'method=?');pars=c(pars,input$data_method)}; if(input$data_chr!='ALL'){wh=c(wh,'chr=?');pars=c(pars,input$data_chr)}
    q(con,paste0('SELECT * FROM segments WHERE ',paste(wh,collapse=' AND '),' ORDER BY chr,start LIMIT 20000'),pars)
  })
  output$data_table <- renderDT(datatable(data_query(),options=list(pageLength=20,scrollX=TRUE)))
  output$download_data <- downloadHandler(filename=function() paste0('gu_segments_',Sys.Date(),'.csv'),content=function(file) fwrite(as.data.table(data_query()),file))
}
