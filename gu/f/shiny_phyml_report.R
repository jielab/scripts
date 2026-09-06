# Candidate-level report prepared automatically by gu.sh final.
gu_phyml_report_ui <- function() {
  tagList(
    tags$p("逐谱系汇总全部位点。选择一行查看该谱系的候选单倍型、指定 lead SNP 与最佳 tag SNP；未通过候选筛选的高相似序列仅作比较。"),
    tags$p(class="text-muted", "最佳 tag SNP 表示对已发现单倍型的标记能力，不等同于 GWAS 因果位点。坐标区间为 0-based half-open，SNP 位置为 1-based。"),
    downloadButton("report_download_summary", "下载全部谱系汇总"),
    downloadButton("report_download_haplotypes", "下载全部单倍型明细"),
    downloadButton("report_download_validation", "下载逐拷贝验证"),
    card(card_header("PhyML · 位点与谱系"), DTOutput("report_loci")),
    card(card_header("所选谱系 · 单倍型与双 SNP 明细"), DTOutput("report_haplotypes")),
    card(card_header("所选单倍型 · 完整指标"), DTOutput("report_details")),
    card(card_header("候选片段的独立方法证据"),
      tags$p("IBDmix：同一个体、同一古人类谱系，单条片段覆盖候选区间 ≥80%。TRACE：同时列出同个体重叠和同一 haplotype 编号的拷贝重叠；Ghost/Unknown 不能验证 Neanderthal/Denisovan 来源。分母为候选携带者，不是可检出区域总量。"),
      tags$p("not_run = 未运行；partial_overlap = 有重叠但覆盖不足80%；not_detected = 未检出重叠；exploratory = 探索性结果（包括 IBDmix chrX）；no_candidate_to_test = 无候选可比较。"),
      DTOutput("report_validation")),
    card(card_header("所选谱系的候选树"),
      tags$p(class="text-muted", "Bootstrap 是该谱系候选边的支持度。没有合格候选边时留空；高序列相似度本身不是 introgression 结论。"),
      verbatimTextOutput("report_tree_note"), plotOutput("report_tree",height="650px"))
  )
}

gu_phyml_report_server <- function(input,output,session,dataset,build,root) {
  read_report <- function(name) {
    path <- file.path(root,paste0(name,".tsv"))
    reactiveFileReader(3000,session,path,function(p) {
      if (!file.exists(p) || file.info(p)$size < 2) return(data.frame())
      tryCatch(read.delim(p,sep="\t",quote="",comment.char="",check.names=FALSE,
               stringsAsFactors=FALSE,na.strings=c("","NA"),fileEncoding="UTF-8"),error=function(e)data.frame())
    })
  }
  all_summary <- read_report("phyml_locus_report")
  all_haps <- read_report("phyml_haplotype_report")
  all_validation <- read_report("phyml_copy_validation")
  all_trees <- read_report("phyml_lineage_trees")
  scope <- function(d) {
    if (!nrow(d)) return(d)
    d[d$dataset_id==dataset() & d$genome_build==build(),,drop=FALSE]
  }
  summary <- reactive(scope(all_summary()))
  selected <- reactive({
    d <- summary(); if(!nrow(d)) return(d)
    i <- input$report_loci_rows_selected
    if(!length(i) || i<1 || i>nrow(d)) i<-1L
    d[i,,drop=FALSE]
  })
  lineage_rows <- function(d) {
    s<-selected(); if(!nrow(s) || !nrow(d)) return(data.frame())
    d[d$locus_key==s$locus_key[[1]] & d$lineage==s$lineage[[1]],,drop=FALSE]
  }
  haps <- reactive(lineage_rows(scope(all_haps())))
  selected_hap <- reactive({
    d<-haps(); if(!nrow(d)) return(d)
    i<-input$report_haplotypes_rows_selected
    if(!length(i) || i<1 || i>nrow(d)) i<-1L
    d[i,,drop=FALSE]
  })
  # Clear haplotype selection when its lineage changes.
  observeEvent(selected(), { DT::selectRows(DT::dataTableProxy("report_haplotypes",session=session),NULL) },ignoreInit=TRUE)
  labels <- c(chr="Chr",locus_id="位点",lineage="谱系",call="PhyML 结论",input_lead_snp="指定 lead SNP",
    best_tag_snp="最佳 tag SNP",best_tag_allele="Tag allele",best_tag_r2="Tag–haplotype r²",input_tag_r2="指定 SNP–haplotype r²",
    n_candidate_haplotypes="候选类型数",n_candidate_copies="候选拷贝数",tree_bootstrap="Bootstrap",tree_purity="树纯度",
    diagnostic_sites="诊断位点数",diagnostic_matches="匹配诊断位点数",candidate_start="候选 start",candidate_end="候选 end",
    n_copies="拷贝数",n_individuals="人数",hap_id="Haplotype",archaic="最匹配参考",prop_match="核心区匹配比例",
    risk_allele="Risk allele",n_risk_copies="Risk 拷贝数",risk_fraction="候选 risk 比例",role="展示角色",
    ibdmix_status="IBDmix 状态",ibdmix_any_overlap_individuals="IBDmix 任意重叠人数",ibdmix_supported_individuals="IBDmix ≥80%人数",ibdmix_individuals="IBDmix 候选人数",
    trace_status="TRACE 状态",trace_any_overlap_individuals="TRACE 任意重叠人数",trace_supported_individuals="TRACE ≥80%人数",trace_individuals="TRACE 候选人数",
    trace_supported_copies="TRACE 同拷贝支持数",trace_candidate_copies="TRACE 候选拷贝数",ils_probability="ILS P")
  show_table <- function(d,columns,select=FALSE,pages=16) {
    if(!nrow(d)) return(datatable(data.frame(状态="无对应结果；请运行 ./gu.sh final"),rownames=FALSE,options=list(dom="t")))
    d<-d[,intersect(columns,names(d)),drop=FALSE]
    # Coordinates and counts must remain exact; round continuous metrics only.
    for(k in names(d)) if(is.numeric(d[[k]]) && any(d[[k]] %% 1 != 0,na.rm=TRUE)) d[[k]]<-signif(d[[k]],7)
    ix<-match(names(d),names(labels));names(d)[!is.na(ix)]<-unname(labels[ix[!is.na(ix)]])
    datatable(d,rownames=FALSE,selection=if(select)list(mode="single",selected=1)else "none",
      options=list(pageLength=pages,scrollX=TRUE,autoWidth=TRUE))
  }
  validation_columns<-c("ibdmix_status","ibdmix_any_overlap_individuals","ibdmix_supported_individuals","ibdmix_individuals","trace_status",
                       "trace_any_overlap_individuals","trace_supported_individuals","trace_individuals","trace_supported_copies","trace_candidate_copies")
  output$report_loci <- renderDT(show_table(summary(),c("chr","locus_id","lineage","call","n_candidate_haplotypes","n_candidate_copies",
    "tree_bootstrap","tree_purity","diagnostic_sites","risk_allele","n_risk_copies","risk_fraction","input_lead_snp","best_tag_snp","best_tag_allele","best_tag_r2","input_tag_r2",validation_columns),TRUE))
  output$report_haplotypes <- renderDT(show_table(haps(),c("hap_id","lineage","role","call","n_copies","n_individuals","archaic","prop_match",
    "diagnostic_matches","diagnostic_sites","candidate_start","candidate_end","ils_probability","tree_bootstrap","input_lead_snp","input_tag_allele","input_tag_r2",
    "best_tag_snp","best_tag_allele","best_tag_r2","best_tag_ppv","best_tag_sensitivity","n_equivalent_best_tags","risk_allele","n_risk_copies","risk_fraction"),TRUE,12))
  output$report_details <- renderDT({
    d<-selected_hap();if(!nrow(d))return(show_table(d,character()))
    d<-data.frame(指标=names(d),值=vapply(d,function(x)if(is.na(x[[1]]))"—" else as.character(x[[1]]),character(1)))
    datatable(d,rownames=FALSE,options=list(pageLength=15,scrollX=TRUE))
  })
  output$report_validation <- renderDT(show_table(haps(),c("hap_id","call",validation_columns),FALSE,12))
  tree <- reactive(lineage_rows(scope(all_trees())))
  output$report_tree_note <- renderText({
    t<-tree();if(!nrow(t))return("没有该谱系的候选树。")
    paste("Lineage:",t$lineage[[1]]," | ",t$tree_call_reason[[1]]," | bootstrap:",t$candidate_clade_bootstrap[[1]],
          " | purity:",t$candidate_purity[[1]]," | sensitivity:",t$candidate_sensitivity[[1]])
  })
  output$report_tree <- renderPlot({
    t<-tree();validate(need(nrow(t)>0,"没有该谱系的候选树。"))
    phy<-ape::read.tree(text=t$tree_newick[[1]])
    tips<-strsplit(as.character(t$candidate_tips_in_clade[[1]]),",",fixed=TRUE)[[1]]
    colours<-ifelse(phy$tip.label %in% tips,"#13815b",ifelse(grepl("Altai|Vindija|Chagyr|Denisova",phy$tip.label),"#b75217","#4d5964"))
    par(mar=c(1,1,1,1));plot(phy,cex=max(.35,min(.8,24/length(phy$tip.label))),tip.color=colours)
    if(!is.null(phy$node.label))ape::nodelabels(phy$node.label,frame="none",cex=.55,adj=c(1.1,-.4))
  })
  download <- function(name,data) downloadHandler(filename=function()paste0(name,".tsv"),content=function(file)write.table(scope(data()),file,sep="\t",quote=FALSE,row.names=FALSE,na=""))
  output$report_download_summary<-download("phyml_locus_report",all_summary)
  output$report_download_haplotypes<-download("phyml_haplotype_report",all_haps)
  output$report_download_validation<-download("phyml_copy_validation",all_validation)
}
