#!/usr/bin/env Rscript
# Portable scientific-contract checks. Uses synthetic records only.
args <- commandArgs(trailingOnly=TRUE)
this <- sub('^--file=', '', grep('^--file=', commandArgs(), value=TRUE)[1L])
root <- if(length(args)) normalizePath(args[1L]) else normalizePath(file.path(dirname(this),'../..'))
engine <- file.path(root,'f/dose')
source(file.path(engine,'scoring.R'))
source(file.path(root,'f/common.R'))
expect_error <- function(code) stopifnot(inherits(try(force(code),silent=TRUE),'try-error'))
extract <- function(path,name) {
    matches <- Filter(function(x) is.call(x) && identical(x[[1L]],as.name('<-')) &&
        identical(x[[2L]],as.name(name)),as.list(parse(path)))
    stopifnot(length(matches)==1L);matches[[1L]]
}

# Explicit boundary expectations include the intentionally non-monotonic domains.
x <- data.frame(eid=letters[1:4],matrix(0,4,length(dose_food_rules()$Food)))
names(x) <- c('eid',dose_food_rules()$Food);attr(x,'intake_unit') <- 'servings/day'
x$Whole_grains <- c(1,1+1e-8,3,3+1e-8)
x$Eggs <- c(0,.5,.5+1e-8,NA_real_)
x$Green_leafy_vegetables <- c(0,1.56,1.56+1e-8,NA_real_)
x$Processed_meats <- c(0,1,1+1e-8,NA_real_)
x$Wine <- c(0,2,2+1e-8,NA_real_)
x$Apples_Pears <- c(0,1,NA,NA);x$Dried_fruit <- c(0,0,1,0)
a <- calculate_dose(x)
stopifnot(identical(a$whole,c(0L,1L,1L,0L)),
          identical(a$eggs,c(0L,1L,0L,NA_integer_)),
          identical(a$leafy,c(0L,1L,0L,NA_integer_)),
          identical(a$processed,c(0L,1L,0L,NA_integer_)),
          identical(a$wine,c(0L,1L,0L,NA_integer_)),
          identical(a$fruit_or,c(0L,1L,1L,NA_integer_)),
          identical(a$fruit_half,c(0,.5,NA_real_,NA_real_)))
versions <- c('dose11','dose9','dose9half')
scores <- setNames(lapply(versions,function(v) score_dose(x,v)),versions)
primary <- make_primary_dose(a,scores)
stopifnot(identical(primary$dose,primary$dose9half),
          identical(primary$dose,primary$diet.dose.sum),
          isTRUE(all.equal(primary$dose,rowSums(primary[attr(primary,'component_columns')]),check.attributes=FALSE)))
bad <- x;bad$eid[2] <- bad$eid[1];expect_error(calculate_dose(bad))
bad <- x;bad$Eggs[1] <- -1;expect_error(calculate_dose(bad))
bad <- x;attr(bad,'intake_unit') <- 'grams/day';expect_error(calculate_dose(bad))
cat('PASS frozen recipe, interval boundaries, half-domain missingness and score aliases\n')

# The score stage must not evaluate even the file-input promises of clinical_merge.
e <- new.env(parent=globalenv())
eval(extract(file.path(engine,'0.common.R'),'dose_cached_phase'),e)
old <- Sys.getenv('DOSE_PREPROCESS_PHASES',unset=NA_character_)
Sys.setenv(DOSE_PREPROCESS_PHASES='qc,medito,webq,patterns,diet_merge,upf')
e$dose_cached_phase('clinical_merge',inputs=stop('Clinical read'),outputs=stop('Clinical write'),code=stop('Clinical execution'))
if(is.na(old))Sys.unsetenv('DOSE_PREPROCESS_PHASES') else Sys.setenv(DOSE_PREPROCESS_PHASES=old)
cat('PASS score-only phase never forces clinical inputs or execution\n')

# Canonical profile reconstruction must preserve independent comparator definitions.
tmp <- tempfile('diet-dose-test-',tmpdir=Sys.getenv('DIET_TEST_TMPDIR',tempdir()))
dir.create(tmp,recursive=TRUE)
part <- primary;part$diet.dash.sum <- c(8,12,20,40)
for(v in versions) for(nm in setdiff(names(scores[[v]]),names(part))) part[[nm]] <- scores[[v]][[nm]]
attr(part,'score_columns') <- c('dose','diet.dose.sum','diet.dash.sum')
attr(part,'recipe') <- 'dose9half'
maha <- data.frame(eid=c('c','a','M'),diet.dash.sum=c(25,30,35))
p1 <- file.path(tmp,'maha.rds');p2 <- file.path(tmp,'dose.rds');saveRDS(maha,p1);saveRDS(part,p2)
can <- diet_publish_disco('ukb',tmp,c(maha=p1,dose=p2))
p <- diet_profile(can,'dose','ukb')
j <- match(part$eid,p$eid)
stopifnot(identical(p$diet.dash.sum[j],part$diet.dash.sum),nrow(p)==nrow(part),
          identical(attr(p,'recipe'),'dose9half'),
          identical(can$diet.dash.sum[match('a',can$eid)],30),
          identical(can$dose.diet.dash.sum[match('a',can$eid)],8))
e <- new.env(parent=globalenv());e$.scripts <- engine;e$.diet_root <- tmp
for(nm in c('.require','.canonical','.overlay','.require_current')) eval(extract(file.path(engine,'adapter.R'),nm),e)
old_file <- Sys.getenv('DIET_DISCO_FILE',unset=NA_character_)
Sys.setenv(DIET_DISCO_FILE=file.path(tmp,'disco/ukb.disco.rds'))
p <- e$.canonical('ukb')
for(v in names(scores)) {
 d <- e$.overlay(scores[[v]],p,setdiff(names(scores[[v]]),'eid'))
 stopifnot(all(vapply(names(d)[names(d) %in% names(scores[[v]])],function(nm)
     isTRUE(all.equal(d[[nm]],scores[[v]][[nm]],check.attributes=FALSE)),logical(1))))
}
changed <- p;changed$dose[1] <- 4.5
updated <- e$.overlay(primary,changed,setdiff(names(primary),'eid'))
stopifnot(updated$dose[match(p$eid[1],updated$eid)]==4.5)
phe <- data.frame(eid=p$eid);attr(phe,'canonical_md5') <- attr(p,'canonical_md5')
dir.create(file.path(tmp,'phe'));saveRDS(list(data=phe),file.path(tmp,'phe/ukb.dose.rds'))
e$.require_current('ukb')
attr(can,'changed') <- TRUE;saveRDS(can,Sys.getenv('DIET_DISCO_FILE'))
expect_error(e$.require_current('ukb'))
if(is.na(old_file))Sys.unsetenv('DIET_DISCO_FILE') else Sys.setenv(DIET_DISCO_FILE=old_file)
# Cycle+SEQN is used to align native NHANES caches, without coercing raw SEQN values.
nh <- data.frame(cycle=c('2003-2004','2005-2006'),SEQN=c(7,7),dose=c(0,0))
nhp <- data.frame(eid=c('2005-2006_7','2003-2004_7'),dose_day1=c(2,5))
nhout <- e$.overlay(nh,nhp,c(dose='dose_day1'))
stopifnot(identical(nhout$dose,c(5,2)),identical(nhout$SEQN,nh$SEQN))
unlink(tmp,recursive=TRUE)
cat('PASS canonical profile names, component identity, score readback, stale-artifact rejection and NHANES keys\n')
for(f in list.files(engine,'[.]R$',full.names=TRUE))invisible(parse(f))
cat('PASS DOSE R syntax\n')
