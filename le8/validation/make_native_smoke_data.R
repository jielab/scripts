# Synthetic native C1-C4 input fixture. No UKB participant data are read.
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
out <- args[1]
if (dir.exists(out)) stop('Choose a new fixture directory')
dir.create(file.path(out, 'Rdata'), recursive = TRUE)
dir.create(file.path(out, 'common'))
set.seed(726)
n <- 2400L
eid <- sprintf('synthetic-%07d', seq_len(n))
features <- c('PCSK9','LPA','GDF15','NTPROBNP','MMP12','APOB','IL6R','CRP','AGRN','TNFRSF4','CPTP','IL6')
g <- matrix(rnorm(n * length(features)), n, dimnames = list(NULL, features))
x <- .3 * g + matrix(rnorm(length(g)), n)
age <- runif(n, 40, 70)
sex <- rbinom(n, 1, .5)
latent <- .65 * x[, 'PCSK9'] + .3 * x[, 'IL6R'] + .02 * (age - 55)
failure <- rexp(n, .045 * exp(latent))
prev <- runif(n) < plogis(-2 + .5 * latent)
baseline <- as.Date('2010-01-01')
diagnosis <- as.Date(rep(NA_character_, n))
incident <- failure < 13 & !prev
diagnosis[incident] <- baseline + pmax(1, round(failure[incident] * 365.25))
diagnosis[prev] <- baseline - sample(1:4500, sum(prev), replace = TRUE)
phe <- data.frame(eid, age, sex, tdi = rnorm(n), PC1 = rnorm(n), PC2 = rnorm(n),
                  center = rep(c('synthetic-site-A','synthetic-site-B'), length.out=n),
                  prot.plate = paste0('synthetic-plate-', (seq_len(n)-1L) %% 20L),
                  ethnic.c = 'White', date_attend = baseline,
                  birth_date = baseline - round(age * 365.25),
                  date_lost = as.Date('2023-12-31'), date_death = as.Date(NA),
                  fod_icd10_cvd_cad = diagnosis, family = rep(seq_len(n/2), each = 2),
                  bmi = pmax(18, 26 + 3*x[, 'PCSK9'] + rnorm(n)), bb_TC = 5 + .7*x[, 'APOB'],
                  bb_HDL = 1.3 + .1*rnorm(n), bb_HBA1C = 38 + 5*x[, 'IL6R'] + rnorm(n),
                  sbp = 125 + 12*x[, 'AGRN'] + rnorm(n), dbp = 78 + 8*rnorm(n),
                  p6153_i0 = sample(c('-7','1','2','3'), n, replace=TRUE),
                  cvd_cad.pgs = .4*g[, 'PCSK9'] + rnorm(n), check.names=FALSE)
behavior <- c(diet='GDF15',pa='MMP12',smoke='NTPROBNP',sleep='LPA')
for (v in names(behavior)) phe[[paste0(v,'.pts')]] <-
    (as.integer(cut(x[,behavior[[v]]] + .3*rnorm(n), c(-Inf,-.8,-.2,.2,.8,Inf))) - 1L) * 25
prot <- data.frame(eid, x, check.names = FALSE)
pgs <- data.frame(eid, g, check.names = FALSE)
saveRDS(phe, file.path(out, 'Rdata/all.rds'))
saveRDS(prot, file.path(out, 'Rdata/prot.rds'))
saveRDS(pgs, file.path(out, 'Rdata/prot.pgs.rds'))
dir.create(file.path(out, 'rap/raw'), recursive = TRUE)
raw <- gzfile(file.path(out, 'rap/raw/prot.tab.gz'), 'wt')
write.table(prot, raw, sep='\t', row.names=FALSE, quote=FALSE)
close(raw)
met_features <- c('Ala','Gly','Val','Tyr','L_VLDL_TG.pct','L_VLDL_TG','Total_TG','ApoB','HDL_C','LDL_C','Glucose','Lactate')
met <- prot; names(met)[-1] <- met_features
met_pgs <- pgs; names(met_pgs)[-1] <- met_features
saveRDS(met, file.path(out, 'Rdata/met.rds'))
saveRDS(met_pgs, file.path(out, 'Rdata/met.pgs.rds'))
write.table(data.frame(data_field=paste0('fixture',seq_along(met_features)),met_name=met_features,
                       full_name=met_features,group=c(rep('Amino_acids',4),rep('Lipoprotein_lipids',6),rep('Energy',2)),
                       subgroup='Synthetic annotation',GCST=NA_character_),
            file.path(out,'common/met.lst'),row.names=FALSE,quote=FALSE,sep='\t')
write.csv(data.frame(eid, group = paste0('synthetic-family-', phe$family)), file.path(out,'groups.csv'), row.names=FALSE)
writeLines(c('SYNTHETIC DATA ONLY; no clinical interpretation.',
             paste('N=', n, 'protein_features=', length(features), 'prevalent=',sum(prev),'incident=',sum(incident))),
           file.path(out, 'SIMULATION_NOTICE.txt'))
cat('Created synthetic phenotype, protein and PGS RDS files:', out, '\n')
