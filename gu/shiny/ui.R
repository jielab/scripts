page_navbar(
  title="GU — Archaic Introgression Browser", id="gu_nav", fillable=FALSE,
  theme=bs_theme(version=5, bootswatch="flatly"),
  header=tags$head(tags$style(HTML(".gu-section-title{font-size:1.22rem;font-weight:700;margin:1.1rem 0 .55rem}.gu-browser-toolbar{display:flex;align-items:center;gap:8px;margin-bottom:8px}.gu-browser-toolbar b{margin-right:auto}.gu-igv-frame{width:100%;height:560px;border:1px solid #ccd3da;border-radius:6px;background:white}table.dataTable tbody tr{cursor:pointer}"))),
  nav_panel("Overview", value="overview",
    card(card_header(tags$b("Interactive genomic context — IGV-Web")),
      layout_columns(numericInput("browser_flank","Context flank (bp)",250000,min=0,max=5000000,step=50000),tags$p(class="text-muted mt-4","Single-click either results table to update this genomic context. Double-click a row to open its detailed viewer."),col_widths=c(3,9)),uiOutput("genome_browser")),
    tags$div(class="gu-section-title","Loci based results"),
    card(tags$p(class="text-muted","Results from loci_avcf. Single-click updates IGV; double-click opens Locus viewer."),DTOutput("overview_loci_detail")),
    tags$div(class="gu-section-title",textOutput("individual_section_title",inline=TRUE)),
    card(layout_columns(selectInput("overview_ind_method","Method",c("ibdmix","trace","as3"),"ibdmix"),selectInput("overview_ind_source","Archaic source",NULL),tags$p(class="text-muted mt-4","Carrier columns compare all available methods in the same 1 Mb region; unavailable analyses are shown as NA."),col_widths=c(2,3,7)),
      card_header(tags$b("Highest-density bins")),DTOutput("overview_individual_bins"),tags$hr(),card_header(tags$b("Genome-wide density and peaks (1 Mb bins)")),plotlyOutput("overview_individual_density",height="420px"))
  ),
  nav_panel("Locus viewer", value="locus_viewer",
    card(card_header(div(style="display:flex;justify-content:space-between",span("Selected locus"),actionButton("locus_back","← Overview",class="btn-sm btn-outline-secondary"))),uiOutput("haplotype_title")),
    card(card_header("Archaic references, matched haplotypes, and negative controls (A/C/G/T)"),layout_columns(numericInput("hap_n_each","Matched modern haplotypes",7,min=1,max=50),numericInput("hap_max_sites","Max polymorphic SNPs",150,min=20,max=1000,step=20),col_widths=c(3,3)),tags$p(class="text-muted","Callable archaic genomes are shown first; source-matching phased haplotypes and ten random non-matched negative controls follow."),uiOutput("haplotype_matrix")),
    accordion(id="hap_details",open=FALSE,accordion_panel("Similarity table / technical details",DTOutput("haplotype_similarity"),tags$hr(),verbatimTextOutput("haplotype_note")))
  ),
  nav_panel("Individual viewer", value="individual_viewer",
    card(card_header(div(style="display:flex;justify-content:space-between",span("Selected individual-based region"),actionButton("individual_back","← Overview",class="btn-sm btn-outline-secondary"))),uiOutput("individual_region_title")),
    layout_columns(value_box(title="Selected-method carriers",value=textOutput("individual_region_carriers",inline=TRUE)),value_box(title="Segment calls",value=textOutput("individual_region_calls",inline=TRUE)),value_box(title="Called sequence",value=textOutput("individual_region_bp",inline=TRUE)),value_box(title="Methods with data",value=textOutput("individual_region_methods",inline=TRUE))),
    card(card_header("Cross-method carrier comparison"),DTOutput("individual_region_comparison")),
    card(card_header("Segments overlapping the selected 1 Mb region"),plotlyOutput("individual_region_plot",height="480px"),DTOutput("individual_region_segments"))
  ),
  nav_panel("Explorer", value="explorer",layout_sidebar(sidebar=sidebar(selectInput("reg_chr","Chromosome",chr_order,"3"),numericInput("reg_start","Start",45000000,min=1),numericInput("reg_end","End",47000000,min=2),selectInput("reg_method","Method","ALL"),sliderInput("reg_limit","Max segments",100,5000,1000,100),actionButton("reg_go","Query")),card(card_header("Regional segment landscape"),plotlyOutput("region_plot",height="540px")),card(card_header("Overlapping loci"),DTOutput("region_loci")))),
  nav_panel("Data", value="data",card(card_header("Database"),verbatimTextOutput("db_info")),card(card_header("Filtered segments"),layout_columns(textInput("data_sample","Sample contains"),selectInput("data_method","Method","ALL"),selectInput("data_chr","Chromosome",c("ALL",chr_order)),col_widths=c(4,4,4)),DTOutput("data_table"),downloadButton("download_data","Download current rows")))
)
