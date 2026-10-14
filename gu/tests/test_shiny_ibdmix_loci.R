# Run from the repository root with the gu R environment.
suppressPackageStartupMessages({library(shiny); library(DT); library(plotly)})
app_dir <- normalizePath('shiny')
source('shiny/evidence.R')
`%||%` <- function(x,y) if(is.null(x) || length(x)==0 || is.na(x) || !nzchar(x)) y else x
fixture <- data.frame(rank=1:132, segment_code=paste0('test-',1:132), chr=c(rep('1',100),rep('2',32)),
    start=seq(0,by=200000,length.out=132), end=seq(100000,by=200000,length.out=132),
    length_bp=100000, source_class='Neanderthal', genome_build='GRCh37',
    n_tested=1000, peak_fraction=.03, peak_raw_fraction=.03, peak_effective_copies=60,
    mean_fraction=.025, background_fraction=.01, extension_threshold=.01, core_threshold=.02,
    core_windows=2, n_windows=2, excess_support_kb=133-(1:132), n_carriers=65,
    n_references=1, reference_names='Altai', max_lod=10, n_calls=65)
cache <- tempfile('gu-evidence-r-test-'); dir.create(cache)
write.table(data.frame(source_class='Neanderthal',chr='1',start=c(0,50000),
    end=c(50000,100000),n_tested=1000,raw_fraction=.03,smooth_fraction=.03),
    file.path(cache,'windows.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
# Reverse input order to verify ranking is applied before the display cap.
gu_ibdmix_evidence_cached <- function(...) list(loci=fixture[132:1,],metadata=list(warnings=character()),directory=cache)
wrapper <- function(input,output,session) {
    seen <- reactiveValues(selected=NULL,opened=NULL)
    visible <- gu_ibdmix_evidence_server(input,output,session,NULL,
        reactive('test'),reactive('GRCh37'),reactive(1),'unused.sqlite',
        function(r,m)seen$selected<-r,function(r,m)seen$opened<-r)
}
shiny::testServer(wrapper, {
    session$setInputs(density_lineage='all',ibdmix_evidence_chr='ALL',ibdmix_evidence_peak=2)
    stopifnot(nrow(visible())==50L,grepl('132',output$ibdmix_evidence_note))
    stopifnot(grepl('Top 50',output$ibdmix_evidence_note),identical(visible()$rank,1:50))
    stopifnot(length(output$ibdmix_evidence_loci)>0,length(output$ibdmix_evidence_profile)>0)
    session$setInputs(ibdmix_evidence_select='test-50',ibdmix_evidence_open='test-50')
    stopifnot(seen$selected$start==9800000,seen$opened$end==9900000)
    session$setInputs(ibdmix_evidence_chr='2')
    stopifnot(nrow(visible())==32L,identical(visible()$rank,101:132))
    session$setInputs(ibdmix_evidence_chr='3')
    stopifnot(nrow(visible())==0L,length(output$ibdmix_evidence_loci)>0)
    session$setInputs(ibdmix_evidence_chr='ALL',density_lineage='Denisovan')
    stopifnot(nrow(visible())==0L,length(output$ibdmix_evidence_loci)>0)
    session$setInputs(density_lineage='all',ibdmix_evidence_peak=4)
    stopifnot(nrow(visible())==0L)
})
unlink(cache,recursive=TRUE)
cat('PASS: strongest 50 after filtering, full candidate set preserved, empty filters, profile, stable selection/open callbacks\n')
