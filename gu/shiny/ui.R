# Application UI
ui <- page_navbar(
  title = "GU — Archaic Introgression Browser", id = "nav", fillable = FALSE,
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  header = tagList(
    tags$script(src="gu_assets/vendor/igv-3.0.0.min.js"),
    tags$script(src="gu_assets/overview.js"),
    tags$script(src="gu_assets/report-notes.js"),
    tags$link(rel="stylesheet",href="gu_assets/report-notes.css"),
    tags$script(HTML("$(document).on('shiny:connected',function(){Shiny.addCustomMessageHandler('gu-region-selected',function(r){document.querySelectorAll('[id^=overview_][id$=curve]').forEach(function(el){if(el.data&&el.layout){Plotly.relayout(el,{shapes:[{type:'rect',xref:'x',yref:'paper',x0:r.start,x1:r.end,y0:0,y1:1,fillcolor:'rgba(240,160,30,.2)',line:{color:'#d99419',width:1}}]});}});});});")),
    tags$head(tags$style(HTML("
      .gu-browser-toolbar{display:flex;align-items:center;gap:.5rem;flex-wrap:wrap;margin-bottom:.6rem}
      .gu-browser-toolbar .gu-locus-label{font-weight:700;margin-right:auto}
      .gu-igv-frame{width:100%;height:330px;border:1px solid #ccd3da;border-radius:6px;background:#fff}
      .gu-phyml-toolbar{display:flex;gap:1rem;align-items:end;flex-wrap:wrap}
      .gu-phyml-toolbar .form-group{margin-bottom:0}
      .gu-metric-help{cursor:help;border-bottom:1px dotted #8a99a8}
      .gu-summary-panel{border-top:1px solid #e1e6e9;padding-top:1.4rem;margin-top:1rem}
      .gu-summary-layout{display:grid;grid-template-columns:minmax(0,2.2fr) minmax(330px,1fr);gap:1.7rem;align-items:start}
      .gu-summary-map-scroll{overflow-x:auto;min-width:0}
      .gu-summary-map{min-width:820px}
      .gu-summary-text{border-left:1px solid #e1e6e9;padding-left:1.5rem;font-size:.88rem}
      .gu-summary-text p{margin-bottom:.75rem}
      .gu-summary-kicker{font-size:.72rem;letter-spacing:.1em;color:#687883;font-weight:700}
      .gu-summary-value{font-size:2.8rem;line-height:1.2;font-weight:700;color:#b3222b;margin-top:.4rem}
      .gu-summary-subvalue{font-size:1.1rem;color:#52616d}
      .gu-summary-scope{background:#f3f6f8;border-radius:5px;padding:.55rem .7rem;color:#435563}
      .gu-summary-table{font-size:.76rem;font-variant-numeric:tabular-nums}
      .gu-summary-table-panel{margin-top:1.2rem;min-width:0}
      .gu-summary-table-scroll{overflow-x:auto}
      .gu-summary-cell-scope{display:block;font-size:.7rem;color:#687883;font-weight:400}
      .gu-summary-table th,.gu-summary-table td{text-align:right;white-space:nowrap}
      .gu-summary-table th:first-child,.gu-summary-table td:first-child{text-align:left}
      .gu-summary-caption{font-size:.8rem;color:#687883;padding:.5rem .4rem}
      @media(max-width:1200px){.gu-summary-layout{grid-template-columns:1fr}.gu-summary-text{border-left:0;border-top:1px solid #e1e6e9;padding:1rem 0 0}}
      .gu-header-filters{display:flex;gap:1rem;align-items:end;flex-wrap:wrap}
      .gu-method-note{border-left:4px solid #18bc9c;padding:.65rem .9rem;background:#f4fbf9;margin:.4rem 0 .8rem}
      .gu-hap-scroll{overflow-x:auto;max-width:100%;padding-bottom:.5rem}
      .gu-hap-table{border-collapse:collapse;font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;font-size:13px;line-height:1}
      .gu-hap-table th.gu-rowlab{position:sticky;left:0;z-index:3;background:#fff;text-align:right;white-space:nowrap;padding:5px 10px 5px 4px;border-right:1px solid #8d99a6;font-family:var(--bs-body-font-family);font-weight:500}
      .gu-hap-table th.gu-site{height:94px;min-width:20px;max-width:20px;padding:0;vertical-align:bottom;color:#607080;font-size:10px;font-weight:400}
      .gu-hap-table th.gu-site span{display:inline-block;transform:rotate(-65deg);transform-origin:bottom left;white-space:nowrap;margin-left:12px}
      .gu-hap-table td.gu-base{min-width:20px;width:20px;height:23px;text-align:center;vertical-align:middle;padding:0;border-right:1px solid #edf0f2;font-weight:800}
      .gu-hap-table td.gu-base-match{background:#dcf4e4;box-shadow:inset 0 -3px #198754}
      .gu-hap-table tr.gu-arch-last th,.gu-hap-table tr.gu-arch-last td{border-bottom:2px solid #34495e}
      .gu-hap-table tr.gu-control-first th,.gu-hap-table tr.gu-control-first td{border-top:2px solid #c0392b}
      .gu-hap-table tr.gu-control th.gu-rowlab{color:#a93226}
      .gu-hap-legend{display:flex;gap:.8rem;align-items:center;flex-wrap:wrap;margin:.2rem 0 .7rem}
      .gu-hap-legend span{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;font-weight:800}
      table.dataTable tbody tr{cursor:pointer}
      .card{height:auto;max-height:none}
      .card > .card-body{flex:0 0 auto;overflow:visible;min-height:0}
      .card .html-widget-output{flex-shrink:0}
    "))),
    div(class = "container-fluid py-2 gu-header-filters",
        selectInput("genome_build", NULL, choices = character(0), width = "220px"),
        selectInput("target_dataset", NULL, choices = character(0), width = "220px")),
    div(class="container-fluid gu-shared-browser",
    gu_card(
      card_header(tags$b("Compact genomic context — IGV-Web / UCSC")),
      layout_columns(
        selectInput("browser_locus", "PhyML locus", choices = character(0)),
        numericInput("browser_flank", "Context flank (bp)", 250000, min = 0, max = 5000000, step = 50000),
        col_widths = c(8, 4)
      ),
      uiOutput("genome_browser"),
      tags$div(id="gu_igv_status",class="text-muted","正在加载 IGV 参考序列与注释…"),
      tags$div(id="gu_igv",style="min-height:300px")
    )
    )
  ),
  nav_panel("Overview", value = "overview",
    downloadButton("report_download_summary", "下载全部 lead 汇总"),
    downloadButton("report_download_haplotypes", "下载全部单倍型明细"),
    downloadButton("report_download_validation", "下载逐拷贝验证"),
    gu_card(card_header("Locus evidence overview — double click to display its phylogeny"),
         selectInput("report_lineage",NULL,choices=c("All"="all","Neanderthal"="Neanderthal","Denisova"="Denisovan"),selected="all",width="180px"),
         DTOutput("report_loci",fill=FALSE), min_height = "540px"),
    gu_card(
      card_header("Haplotype matching methods"),
      div(style="display:flex;align-items:center;gap:12px;flex-wrap:wrap",
        conditionalPanel("input.density_method == 'ibdmix'",
          selectInput("density_lineage",NULL,choices=c("All"="all","Neanderthal"="Neanderthal","Denisova"="Denisovan"),selected="all",width="180px")),
        selectInput("density_method",NULL,choices=c("IBDmix"="ibdmix","Trace"="trace","AS3"="as3"),selected="ibdmix",width="260px"),
        uiOutput("density_chromosomes",inline=TRUE)),
      conditionalPanel("input.density_method == 'ibdmix'",
      tags$h5("IBDmix · Genome-wide overview"),
      plotlyOutput("introgression_density", height="860px")),
      conditionalPanel("input.density_method == 'trace'",gu_matching_ui("overview_trace","Trace · Selected region")),
      conditionalPanel("input.density_method == 'as3'",gu_matching_ui("overview_as3","AS3 · Selected region")),
      conditionalPanel("input.density_method == 'ibdmix'",
        div(class="gu-summary-panel",
          tags$h5("IBDmix · Genome-wide summary"),
          tags$p(class="text-muted","当前数据集的静态人群汇总，不随上方区间缩放改变。"),
          div(class="gu-summary-layout",
            div(class="gu-summary-map-scroll",
              div(class="gu-summary-map",plotOutput("ibdmix_summary_map",height="620px")),
              tags$p(class="gu-summary-caption","圆饼红色扇区为人群平均 Altai Neanderthal 覆盖比例（0–100%）；旁注给出百分比数值。地点为采样地或祖籍的近似位置，南亚侨居人群按祖籍定位；引线仅用于避让标签。底图：Natural Earth。")),
            div(class="gu-summary-text",uiOutput("ibdmix_summary_text"))),
          div(class="gu-summary-table-panel",uiOutput("ibdmix_summary_table"))))
    )
  ),
  nav_panel("PhyML", value = "phyml",
    gu_card(
      card_header(uiOutput("haplotype_title")),
      div(class="gu-phyml-toolbar",
        selectInput("phyml_locus", "位点", choices=character(0), width="280px"),
        selectInput("phyml_tree_choice", "谱系 / 系统树", choices=character(0), width="360px"),
        selectInput("tree_display_min_copies", "显示单倍型", choices=c("拷贝数 > 10"="11","拷贝数 ≥ 2"="2"),selected="11",width="180px"),
        downloadButton("download_panel_b_pdf","下载树图 PDF")),
      plotOutput("phyml_tree", height="800px")
    ),
    gu_phyml_report_ui(),
    accordion(id="phyml_details",open=FALSE,
      accordion_panel("全部位点",
        layout_columns(
          selectInput("locus_chr", "Chromosome", c("ALL", chr_order), "ALL"),
          selectInput("locus_status", "Run status", c("ALL", "pass"), "ALL")),
        DTOutput("loci_table",fill=FALSE)),
      accordion_panel("运行与数据明细",
        verbatimTextOutput("tree_summary"),
        verbatimTextOutput("haplotype_note"),
        DTOutput("haplotype_table",fill=FALSE)))
  ),
  nav_panel("IBDmix", value="ibdmix", gu_method_ui("ibdmix")),
  nav_panel("Trace", value="trace", gu_method_ui("trace")),
  nav_panel("AS3", value="as3", gu_method_ui("as3"))
)
