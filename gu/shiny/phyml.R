# Original COJO lead is the only index SNP. Counts refer to its risk haplotypes.
gu_phyml_overview <- function(summary, haps, validation=data.frame()) {
  if (!nrow(summary)) return(summary)
  summary$record_id <- paste(summary$locus_key,summary$lineage,sep="|")
  summary$archaic_LD_match <- "未评估"
  refs<-c("Vindija","Altai","Chagyr","Denisova","Denisova25")
  {
    summary$archaic_LD_match<-vapply(seq_len(nrow(summary)),function(i)paste(vapply(if(summary$lineage[i] %in% c("Denisovan","Denisova"))refs[4:5]else refs[1:3],function(r) {
      if(!all(paste0(r,c("_LD_matches","_LD_called")) %in% names(summary)))return("—")
      a<-summary[[paste0(r,"_LD_matches")]][i];b<-summary[[paste0(r,"_LD_called")]][i]
      if(is.na(a)||is.na(b))"—"else paste0(a,"/",b)
    },character(1)),collapse=" · "),character(1))
  }
  summary$risk_haplotypes <- ifelse(is.na(summary$n_candidate_haplotypes)|is.na(summary$n_candidate_copies),"未评估",paste0(summary$n_candidate_haplotypes," / ",summary$n_candidate_copies))
  summary$ibdmix_risk_support <- "未评估"
  for(i in seq_len(nrow(summary))) {
    if("ibdmix_status" %in% names(summary) && summary$ibdmix_status[i] %in% "not_run")summary$ibdmix_risk_support[i]<-"未运行"
    if(!nrow(validation))next
    v<-validation[validation$locus_key==summary$locus_key[i] & validation$lineage==summary$lineage[i] & validation$method=="ibdmix",,drop=FALSE]
    if(!nrow(v))next
    summary$ibdmix_risk_support[i]<-if(any(v$method_complete %in% c(1,"TRUE")))"未评估"else "未运行"
    good<-v$method_complete %in% c(1,"TRUE") & v$comparison_available %in% c(1,"TRUE") & v$evidence_eligible %in% c(1,"TRUE")
    ids<-unique(v$sample_id[good]);yes<-unique(v$sample_id[good & v$overlap_pass %in% c(1,"TRUE")])
    if(length(ids))summary$ibdmix_risk_support[i]<-paste0(length(yes)," / ",length(ids),if(length(setdiff(unique(v$sample_id),ids)))"（部分可评估）"else "")
  }
  summary
}

gu_phyml_report_ui <- function() {
  bslib::navset_card_tab(
    bslib::nav_panel("单倍型",
      tags$p(class="gu-overview-note", "每行一种重复出现的核心区间序列。risk 携带原始 GWAS lead 的风险等位基因；nonrisk 为同一棵树中的对照。"),
      DTOutput("report_haplotypes",fill=FALSE)),
    bslib::nav_panel("序列比对",
      layout_columns(
        numericInput("hap_n_match", "风险单倍型", 8, min=1, max=50),
        numericInput("hap_n_control", "对照单倍型", 10, min=0, max=50),
        numericInput("hap_max_sites", "最多显示 SNP", 150, min=20, max=1000, step=20)),
      uiOutput("haplotype_matrix"),
      accordion(open=FALSE,accordion_panel("匹配统计",DTOutput("haplotype_similarity",fill=FALSE)))),
    bslib::nav_panel("方法交叉验证",
      tags$span(class="gu-metric-help",tabindex="0",title="IBDmix：同一个体、同一古人类谱系，单条片段覆盖候选区间 ≥80%。TRACE 同时报告同个体与同拷贝的支持；Ghost/Unknown 不用于确认 Neanderthal/Denisovan 来源。分母为候选携带者。not_run 表示未运行，partial_overlap 表示重叠不足80%，not_detected 表示未检出，exploratory 表示探索性结果。", "验证口径 ⓘ"),
      DTOutput("report_validation",fill=FALSE)),
    bslib::nav_panel("携带者",DTOutput("haplotype_samples",fill=FALSE)),
    bslib::nav_panel("详细指标",DTOutput("report_details",fill=FALSE))
  )
}

gu_phyml_report_server <- function(input,output,session,dataset,build,root,external_record=NULL,external_region=NULL,external_lineage=NULL,external_locus=NULL) {
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
  scope <- function(d) {
    if (!nrow(d)) return(d)
    d[d$dataset_id==dataset() & d$genome_build==build(),,drop=FALSE]
  }
  overview <- reactive(gu_phyml_overview(scope(all_summary()),scope(all_haps()),scope(all_validation())))
  # Review completed tree tests, including those without tree support.
  summary <- reactive({
    d <- overview()
    if(!nrow(d)) return(d)
    lineage <- input$report_lineage
    if(is.null(lineage)) lineage <- "all"
    d[d$call %in% c("tree_supported","tree_not_supported") &
      (lineage=="all" | d$lineage %in% if(lineage=="Denisovan")c("Denisovan","Denisova")else lineage),,drop=FALSE]
  })
  selected_index<-reactiveVal(1L)
  observeEvent(input$report_lineage, {
    d<-summary(); i<-if(!is.null(external_record))match(external_record(),d$record_id)else NA_integer_
    selected_index(if(length(i)==1L && !is.na(i))i else if(nrow(d))1L else NA_integer_)
  },priority=110)
  select_record <- function(value) {
    i<-suppressWarnings(as.integer(value))
    if(length(i)==1L && !is.na(i) && i>=1L && i<=nrow(summary())) selected_index(i)
  }
  observeEvent(input$report_locus_select,select_record(input$report_locus_select),priority=100)
  observeEvent(input$report_locus_open,select_record(input$report_locus_open),priority=100)
  if(!is.null(external_record)) observeEvent(external_record(), {
    d<-summary();i<-match(external_record(),d$record_id)
    if(length(i)==1L && !is.na(i))selected_index(i)
  },priority=100)
  if(!is.null(external_region)) observeEvent(external_region(), {
    r<-external_region();req(!is.null(r))
    d<-summary()
    ix<-which(as.character(d$chr)==as.character(r$chr[[1]]) & d$core_start<r$end[[1]] & d$core_end>r$start[[1]])
    # A viewport changes the highlight, never the overview's full set of rows.
    # Preserve the exact locus/lineage when its locus still overlaps.
    current<-selected_index()
    i<-if(current %in% ix)current else if(length(ix))ix[[1]] else NA_integer_
    selected_index(i)
  },priority=100)
  if(!is.null(external_lineage)) observeEvent(external_lineage(), {
    d<-summary()
    if(!is.null(external_locus)) {
      r<-external_locus();if(is.null(r) || !nrow(r))return()
      ix<-which(d$locus_id==r$locus_id[[1]] & d$lineage==external_lineage())
    } else {
      r<-if(!is.null(external_region))external_region()else NULL
      if(is.null(r))return()
      ix<-which(as.character(d$chr)==as.character(r$chr[[1]]) & d$core_start<r$end[[1]] & d$core_end>r$start[[1]] & d$lineage==external_lineage())
    }
    if(!selected_index() %in% ix)selected_index(if(length(ix))ix[[1]] else NA_integer_)
  },ignoreInit=TRUE,priority=100)
  observeEvent(selected_index(), {
    i<-selected_index()
    DT::selectRows(DT::dataTableProxy("report_loci",session=session),if(is.na(i))NULL else i)
  },ignoreInit=TRUE)
  selected <- reactive({
    d<-summary();if(!nrow(d))return(d)
    i<-selected_index()
    if(is.na(i)||i<1||i>nrow(d))return(d[0,,drop=FALSE])
    d[i,,drop=FALSE]
  })
  lineage_rows <- function(d) {
    if(!is.null(external_locus)) {
      r<-external_locus();lineage<-if(!is.null(external_lineage))external_lineage()else NULL
      if(is.null(r) || !nrow(r) || !nrow(d) || !length(lineage))return(data.frame())
      return(d[d$locus_id==r$locus_id[[1]] & (lineage=="all" | d$lineage==lineage),,drop=FALSE])
    }
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
  observeEvent(list(selected(),haps()), {
    s<-selected();d<-haps();id<-if(nrow(s) && "hap_id" %in% names(s))s$hap_id[[1]]else NA
    if((is.na(id) || !nzchar(id)) && nrow(s) && "best_haplotype" %in% names(s))id<-s$best_haplotype[[1]]
    i<-if(nrow(d))match(id,d$hap_id)else NA
    DT::selectRows(DT::dataTableProxy("report_haplotypes",session=session),if(!nrow(d))NULL else if(is.na(i))1L else i)
  },ignoreInit=TRUE)
  labels <- c(chr="Chr",core_interval="核心区间 (GRCh37)",lineage="检验谱系",index_snp="Index SNP",
    risk_allele="Risk Allele",p_j="COJO P",core_kb="核心 (kb)",n_ld_sites="LD SNP数",n_sites="建树位点",
    risk_haplotypes="单倍型",tree_bootstrap="Bootstrap",call="树检验",risk_frequency_EUR="EUR 风险频率",
    archaic_LD_match="古参考匹配",n_nonrisk_haplotypes="对照类型",n_singleton_copies="单次拷贝",search_edge_warning="边界提示",ils_probability="ILS P（模型）",
    ibdmix_risk_support="IBDmix 支持人数",hap_id="Haplotype",role="类别",n_copies="拷贝数",n_individuals="人数",
    archaic="最相似参考",prop_match="序列相同 (%)",n_compared="可比较位点",n_match="相同位点",reason="原因",
    ibdmix_status="IBDmix 状态",ibdmix_supported_individuals="IBDmix 支持人数",ibdmix_individuals="可评估人数",
    trace_status="TRACE 状态",trace_supported_individuals="TRACE 支持人数",trace_individuals="可评估人数")
  help <- c(
    chr="分析坐标使用 GRCh37。",
    core_interval="与原始 lead 在 1KG EUR 中 phased r² > 0.98 的 SNP 所覆盖的区间，含两端。不是固定 1 Mb；搜索范围默认 lead 两侧各 500 kb。",
    lineage="检验风险单倍型与所选谱系的参考是否共同成支（Neanderthal 三个、Denisovan 两个），并排除另一谱系参考；不表示这一行已经证实古人类来源。",
    index_snp="原始 COJO 文件的 lead，不重新挑选标记。名称内坐标保留输入基因组版本；核心区间和实际分析坐标为 GRCh37。",
    risk_allele="GRCh37 方向的风险等位基因：bJ > 0 时为 refA，bJ < 0 时为另一等位基因。",
    p_j="输入 COJO 的联合分析 P 值 pJ；0 表示原文件数值下溢，不代表概率精确为零。",
    n_ld_sites="核心边界由这些与 lead 的 EUR phased r² > 0.98 的 SNP 定义。不是树中只使用这些 SNP。",
    n_sites="核心区间内可在五个古人类参考中共同比较、且现代样本次要等位基因至少出现两次的 SNP 数。包括不满足高 LD 阈值的区间内位点。",
    risk_haplotypes="重复出现的风险单倍型：序列种类数 / 染色体拷贝数。只按原始 lead 的风险等位基因定义，不按古人类相似度挑选；每种序列至少出现两次。不是人数，也不表示全部获得树支持。",
    tree_bootstrap="预先定义的风险单倍型与所选谱系参考共同成支、且不含另一谱系参考、非风险对照或祖先序列的分支支持率；100 次重采样。空白表示没有该分支或未完成树，不是零，也不是渗入概率。",
    call="树支持表示预定义风险分支 Bootstrap ≥70（程序报告阈值）。不等于 high confidence introgression；仍需结合序列、重组及 IBDmix 等证据。此概览仅显示已完成树检验的 lead，包括树支持和未获树支持；无法评估、未建树及建树失败的记录不显示。",
    risk_frequency_EUR="1KG EUR 中风险等位基因的染色体频率；不是风险单倍型在人群中的疾病效应。",
    core_kb="高 LD 核心区间长度，单位 kb。",
    ils_probability="仅为长度模型敏感性指标：假设重组率 0.53 cM/Mb、世代 29 年、分化 55 万年、古人类年龄 5 万年。未使用该 locus 的局部重组图谱，也未作多重检验校正，不能直接据此宣布高置信渗入。",
    ibdmix_risk_support="本行重复风险单倍型携带者中，同个体、同谱系 IBDmix 单条片段覆盖核心区间 ≥80% 的人数 / 可评估人数，每人只计一次。无论树是否支持均检验；IBDmix 不保证同一染色体拷贝。来自 final/review；IBDmix 更新后运行 gu.sh final。未运行或未评估不等于零。",
    role="risk：风险单倍型；nonrisk：非风险对照；mixed：相同核心序列无法区分 lead 的两种等位基因。",
    n_copies="该序列的染色体拷贝数；一个常染色体个体可贡献两份。",
    n_individuals="携带该序列的人数，每人只计一次。",
    prop_match="可比较的核心 SNP 中，与该参考相同的比例；不是全基因组一致率，也不是渗入概率。",
    archaic_LD_match="Neanderthal 依次为 Vindija · Altai · Chagyr；Denisovan 依次为 Denisova · Denisova25：高 LD SNP 中与风险相关等位基因相同的位点数 / 该参考可比较的位点数。0/0 表示无可比较位点；相同不等于渗入。",
    n_nonrisk_haplotypes="同一棵树内重复出现的非风险单倍型种类数；不按相似度抽样。",
    n_singleton_copies="完整核心序列仅出现一次，按论文规则未纳入树的拷贝数。",
    search_edge_warning="1 表示高 LD 标记接近搜索边界，核心可能不完整，需扩大 --phyml-window-bp 后复核。0 仅表示未接近边界，不保证不存在更远的 LD。")
  show_table <- function(d,columns,select=FALSE,pages=16,selected_row=1L,full=FALSE,length_change=TRUE) {
    if(!nrow(d)) return(datatable(data.frame(状态="暂无该范围的结果"),rownames=FALSE,options=list(dom="t")))
    row_colours<-if("lineage" %in% names(d))ifelse(d$lineage=="Neanderthal","#2878B5",ifelse(d$lineage %in% c("Denisovan","Denisova"),"#8B5A2B","inherit"))else rep("inherit",nrow(d))
    d<-d[,intersect(columns,names(d)),drop=FALSE]
    if("call" %in% names(d)) {
      calls<-c(tree_supported="树支持",tree_not_supported="未获树支持",tree_not_requested="未建树",tree_failed="建树失败",risk_haplotype="风险单倍型",nonrisk_control="非风险对照",not_evaluable="未评估",insufficient_high_LD_markers="高 LD 位点不足",lead_absent_or_ambiguous="lead 缺失或不唯一",lead_allele_mismatch="等位基因不一致",ancestral_sequence_unavailable="祖先序列缺失",insufficient_recurrent_haplotypes="重复单倍型不足",insufficient_tree_sites="建树位点不足",risk_nonrisk_sequence_unresolved="风险序列无法区分")
      ix<-match(d$call,names(calls));d$call[!is.na(ix)]<-unname(calls[ix[!is.na(ix)]])
    }
    if("role" %in% names(d)) {
      d$role<-gsub("input_marker_representative","指定 SNP 代表",d$role,fixed=TRUE)
      d$role<-gsub("raw_similarity_best","最高原始相似度",d$role,fixed=TRUE)
      d$role<-gsub("sequence_candidate","序列候选",d$role,fixed=TRUE)
      d$role<-gsub(";","；",d$role,fixed=TRUE)
    }
    for(k in intersect(c("input_lead_snp","best_tag_snp"),names(d)))d[[k]][is.na(d[[k]])|!nzchar(d[[k]])]<-"未找到"
    for(k in intersect(c("tree_purity","tree_sensitivity"),names(d)))d[[k]]<-round(100*as.numeric(d[[k]]),1)
    if("prop_match" %in% names(d))d$prop_match<-round(100*as.numeric(d$prop_match),2)
    for(k in intersect(c("input_vs_best_snp_r2"),names(d)))d[[k]]<-round(as.numeric(d[[k]]),3)
    # Coordinates and counts must remain exact; round continuous metrics only.
    for(k in names(d)) if(is.numeric(d[[k]]) && any(d[[k]] %% 1 != 0,na.rm=TRUE)) d[[k]]<-signif(d[[k]],7)
    # Format after rounding so JSON/browser conversion cannot expose float tails.
    if("p_j" %in% names(d)) d$p_j <- vapply(as.numeric(d$p_j),function(p) {
      if(is.na(p)) NA_character_ else as.character(signif(p,2))
    },character(1))
    tips<-vapply(names(d),function(k)if(k %in% names(help))help[[k]]else "",character(1))
    header_callback<-DT::JS(paste0("function(thead){var tips=",jsonlite::toJSON(unname(tips)),";$(thead).find('th').each(function(i){if(tips[i]){this.removeAttribute('title');if(!this.querySelector('.gu-note-mark')){var mark=document.createElement('span');mark.className='gu-note-mark';mark.tabIndex=0;mark.setAttribute('role','button');mark.setAttribute('aria-label','查看列说明');mark.textContent='▤';this.appendChild(mark);}this.querySelector('.gu-note-mark').setAttribute('data-gu-note',tips[i]);}});}"))
    ix<-match(names(d),names(labels));names(d)[!is.na(ix)]<-unname(labels[ix[!is.na(ix)]])
    core_col<-match("核心区间 (GRCh37)",names(d))-1L
    row_callback<-if(!is.na(core_col))DT::JS(paste0("function(row,data,displayNum,displayIndex,index){var colours=",jsonlite::toJSON(unname(row_colours)),";$('td',row).eq(",core_col,").css({'color':colours[index],'font-weight':'600'});}"))else NULL
    p_col<-match("COJO P",names(d))-1L
    column_defs<-if(!is.na(p_col))list(list(targets=p_col,width="72px",className="gu-cojo-p",type="num"))else list()
    datatable(d,rownames=FALSE,class="compact stripe hover nowrap",selection=if(select)list(mode="single",selected=if(length(selected_row) && !is.na(selected_row))selected_row else NULL)else "none",
      options=list(pageLength=pages,paging=!full,lengthChange=length_change,dom=if(full)"t"else if(length_change)"lftip"else "ftip",displayStart=if(select && length(selected_row) && !is.na(selected_row))floor((selected_row-1L)/pages)*pages else 0L,scrollX=TRUE,autoWidth=FALSE,columnDefs=column_defs,rowCallback=row_callback,headerCallback=header_callback),
      callback=DT::JS("table.on('click','tbody tr',function(){if(table.table().node().closest('#report_loci')){var i=table.row(this).index();if(i!==undefined)Shiny.setInputValue('report_locus_select',i+1,{priority:'event'});}});table.on('dblclick','tbody tr',function(){if(table.table().node().closest('#report_loci')){var i=table.row(this).index();if(i!==undefined)Shiny.setInputValue('report_locus_open',i+1,{priority:'event'});}});"))
  }
  validation_columns<-c("ibdmix_status","ibdmix_any_overlap_individuals","ibdmix_supported_individuals","ibdmix_individuals","trace_status",
                       "trace_any_overlap_individuals","trace_supported_individuals","trace_individuals","trace_supported_copies","trace_candidate_copies")
  output$report_loci <- renderDT(show_table(summary(),c("chr","core_interval","lineage","index_snp","risk_allele","p_j","core_kb","n_ld_sites","n_sites","archaic_LD_match","risk_haplotypes",
    "tree_bootstrap","call","reason","search_edge_warning","ils_probability","ibdmix_risk_support"),TRUE,pages=10,selected_row=isolate(selected_index()),length_change=FALSE),server=FALSE)
  output$report_haplotypes <- renderDT(show_table(haps(),c("hap_id","role","call","n_copies","n_individuals","archaic","prop_match","n_compared","n_match"),TRUE,12))
  output$report_details <- renderDT({
    d<-selected_hap();if(!nrow(d))return(show_table(d,character()))
    d<-data.frame(指标=names(d),值=vapply(d,function(x)if(is.na(x[[1]]))"—" else as.character(x[[1]]),character(1)))
    datatable(d,rownames=FALSE,options=list(pageLength=15,scrollX=TRUE))
  })
  output$report_validation <- renderDT(show_table(haps(),c("hap_id","call",validation_columns),FALSE,12))
  download <- function(name,data) downloadHandler(filename=function()paste0(name,".tsv"),content=function(file)write.table(scope(data()),file,sep="\t",quote=FALSE,row.names=FALSE,na=""))
  output$report_download_summary<-download("phyml_locus_overview",overview)
  output$report_download_haplotypes<-download("phyml_haplotype_report",all_haps)
  output$report_download_validation<-download("phyml_copy_validation",all_validation)
  selected
}
