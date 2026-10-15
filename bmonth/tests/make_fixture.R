# Small, deterministic input covering all stages without private cohort data.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
dest <- file.path(args[[1]], 'Rdata')
dir.create(dest, recursive = TRUE, showWarnings = FALSE)
set.seed(20261014)
n <- 3000L
month <- rep(1:12, length.out = n)
x <- -cos(2*pi*month/12)
z <- sin(2*pi*month/12)
mediator <- .9*x + .5*z + rnorm(n)
dat <- data.frame(
  eid = as.character(seq_len(n)), birth_month = month,
  birth_north = runif(n, 100000, 1000000),
  birth_east = runif(n, 100000, 600000),
  age = runif(n, 40, 70), sex = sample(c('Female', 'Male'), n, TRUE),
  edu.sco = sample(0:5, n, TRUE), tdi = rnorm(n),
  center = sample(c('A', 'B', 'C'), n, TRUE),
  date_attend = as.Date('2010-01-01'), date_lost = as.Date(NA),
  date_death = as.Date(NA)
)
event_time <- rexp(n, .065*exp(.5*mediator + .15*x))
dat$cvd_cad.t2e <- pmin(event_time, 10)
dat$cvd_cad.Yt2e <- as.integer(event_time <= 10)
# The second outcome exercises date reconstruction and prevalent-case exclusion.
dat$fod_icd10_ra <- as.Date('2010-01-01') + round(rexp(n, .03*exp(.4*mediator))*365.25)
dat$fod_icd10_ra[1:20] <- as.Date('2008-01-01')
dat$date_lost[21:40] <- as.Date('2014-01-01')
dat$date_death[41:60] <- as.Date('2015-01-01')
dat$birth_month[61:66] <- c(-1, -3, 0, 13, NA, 1.5)
dat$birth_north[67:70] <- c(-1, -3, NA, Inf)
dat$birth_east[71:74] <- c(-1, -3, NA, Inf)
dat$tdi[75:78] <- NA
saveRDS(dat, file.path(dest, 'all.rds'))
ph <- data.frame(id = dat$eid, check.names = FALSE)
ph[['250.2']] <- runif(n) < plogis(-1 + .4*x)
ph[['401']] <- runif(n) < plogis(-1 + .3*z)
ph[['714']] <- runif(n) < .2
ph[['493']] <- runif(n) < .15
ph[['008']] <- seq_len(n) <= 5
ph[['008.5']] <- FALSE
ph[['008.51']] <- NA
ph[['401']][79:85] <- NA
saveRDS(ph, file.path(dest, 'PheWAS.rds'))
saveRDS(data.frame(eid=dat$eid, signal=mediator, noise=rnorm(n), constant=1),
        file.path(dest, 'prot.rds'))
saveRDS(data.frame(eid=dat$eid, signal=.8*mediator+rnorm(n, sd=.5), noise=rnorm(n), constant=1),
        file.path(dest, 'met.rds'))
cat('Fixture written:', n, 'participants,', ncol(ph)-1, 'phenotypes\n')
