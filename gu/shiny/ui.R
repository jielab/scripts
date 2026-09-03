page_navbar(
  title = "GU — Archaic Introgression Browser", id = "nav", fillable = FALSE,
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  header = tagList(
    tags$head(tags$style(HTML("
      .gu-browser-toolbar{display:flex;align-items:center;gap:.5rem;flex-wrap:wrap;margin-bottom:.6rem}
      .gu-browser-toolbar .gu-locus-label{font-weight:700;margin-right:auto}
      .gu-igv-frame{width:100%;height:330px;border:1px solid #ccd3da;border-radius:6px;background:#fff}
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
    "))),
    div(class = "container-fluid py-2 gu-header-filters",
        selectInput("genome_build", "Genome build", choices = character(0), width = "220px"),
        selectInput("target_dataset", "Target dataset", choices = character(0), width = "220px"))
  ),
  nav_panel("Overview", value = "overview",
    card(
      card_header(tags$b("Compact genomic context — IGV-Web / UCSC")),
      layout_columns(
        selectInput("browser_locus", "PhyML locus", choices = character(0)),
        numericInput("browser_flank", "Context flank (bp)", 250000, min = 0, max = 5000000, step = 50000),
        selectInput("support_population", "Population", "ALL"),
        selectInput("support_denominator", "Curve denominator",
                    choices = c("PhyML-positive"="phyml_positive", "All PhyML-tested"="tested"),
                    selected = "phyml_positive"),
        tags$p(class = "text-muted mt-4",
               "IGV supplies gene context. GU evidence is shown in the aligned prevalence curves below; download the bedGraph to add those tracks to IGV."),
        col_widths = c(3, 2, 2, 2, 3)
      ),
      uiOutput("genome_browser")
    ),
    layout_columns(
      value_box("Samples tested", textOutput("n_samples", inline = TRUE)),
      value_box("Tree-clade carriers", textOutput("n_phyml_positive", inline = TRUE)),
      value_box("Independent overlap", textOutput("n_other_supported", inline = TRUE)),
      value_box("Selected-locus state", textOutput("locus_evidence_state", inline = TRUE)),
      value_box("Methods completed here", textOutput("n_methods", inline = TRUE))
    ),
    card(card_header("Five-locus validation — click to update browser; double-click for sequence view"),
         DTOutput("overview_loci_table")),
    card(
      card_header("Selected-locus interpretation"),
      tags$div(class = "gu-method-note",
               "PhyML carriers are defined by a bootstrap-supported archaic/modern tree edge. Pairwise sequence identity is descriptive and is not used as the carrier rule. IBDmix, TRACE, and AS3 remain separate evidence streams; not run, unsupported, and exploratory are never converted to zero evidence."),
      uiOutput("locus_interpretation"),
      downloadButton("download_locus_evidence", "Download locus evidence")
    ),
    card(card_header("Independent method support"), DTOutput("locus_support_table")),
    accordion(
      id = "overview_details", open = FALSE,
      accordion_panel("Regional prevalence curves", plotlyOutput("support_trajectory", height = "430px")),
      accordion_panel("Individual support rows",
        layout_columns(
          checkboxInput("support_phyml_only", "Show tree-clade carriers only", TRUE),
          downloadButton("download_locus_tracks", "Download IGV bedGraph"),
          downloadButton("download_locus_support", "Download individual support"),
          col_widths = c(5, 3, 4)
        ),
        tags$p(class = "text-muted", "The interactive table is capped at 2,000 rows; the download contains every matching individual."),
        DTOutput("locus_sample_support"))
    )
  ),
  nav_panel("PhyML loci", value = "phyml",
    layout_sidebar(
      sidebar(
        selectInput("locus_chr", "Chromosome", c("ALL", chr_order), "ALL"),
        selectInput("locus_status", "Status", c("ALL", "pass"), "ALL"),
        numericInput("hap_n_match", "Matched haplotypes", 8, min = 1, max = 50),
        numericInput("hap_n_control", "Non-matched controls", 10, min = 0, max = 50),
        numericInput("hap_max_sites", "Maximum SNPs", 150, min = 20, max = 1000, step = 20),
        open = "open"
      ),
      card(card_header("Tree-defined modern–archaic haplotype evidence"), DTOutput("loci_table")),
      card(card_header("Selected locus"), uiOutput("haplotype_title")),
      card(
        card_header("Archaic references, candidate-clade modern haplotypes, and controls (A/C/G/T)"),
        tags$p(class = "text-muted",
               "Archaic references are shown first. A green cell background marks a modern base matching that row's best archaic reference; the letter colour identifies A/C/G/T."),
        uiOutput("haplotype_matrix")
      ),
      accordion(
        id = "phyml_details", open = FALSE,
        accordion_panel("Match statistics", DTOutput("haplotype_similarity")),
        accordion_panel("Raw haplotype groups", DTOutput("haplotype_table")),
        accordion_panel("Display provenance / technical details", verbatimTextOutput("haplotype_note"))
      ),
      card(card_header("PhyML candidate-clade tree for selected locus"),
           tags$p(class = "text-muted", "The shaded edge-side is selected directly from the tree: one coherent archaic lineage, at least two modern haplotypes, and bootstrap ≥70 for a passing clade. Without an explicit ancestral outgroup this is an unrooted ML topology; display rooting does not imply ancestry direction."),
           verbatimTextOutput("tree_summary"), plotOutput("phyml_tree", height = "520px")),
      card(card_header("Samples carrying selected haplotypes"), DTOutput("haplotype_samples"))
    )
  ),
  nav_panel("Segments", value = "segments",
    layout_sidebar(
      sidebar(selectInput("seg_method", "Method", "ALL"), selectInput("seg_chr", "Chromosome", c("ALL", chr_order), "ALL"),
              numericInput("seg_start", "Start (0-based)", 0, min = 0), numericInput("seg_end", "End", 200000000, min = 1),
              checkboxInput("show_reference", "Overlay published callset", TRUE),
              selectInput("reference_population", "Published population", "ALL"),
              sliderInput("seg_limit", "Maximum rows", 100, 10000, 2000, 100), actionButton("seg_go", "Query"), open = "open"),
      card(card_header("Regional segment landscape"), plotlyOutput("segment_plot", height = "470px")),
      card(card_header("Segment calls — click a row to update the Overview browser"), DTOutput("segment_table")),
      card(card_header("Published callset overlay — external_reference, never model input"), DTOutput("published_callset_table")),
      card(card_header("GU-AS3 overlap / reciprocal-overlap by 1KG population"), DTOutput("reference_overlap_table"))
    )
  ),
  nav_panel("Individuals", value = "individuals",
    card(card_header("Per-sample burden"),
         layout_columns(textInput("sample_search", "Sample contains", ""), selectInput("sample_method", "Method", "ALL"),
                        selectInput("burden_type", "Burden definition", c("raw_call", "nonredundant_union", "consensus_catalog"), "nonredundant_union"), col_widths = c(4, 4, 4)),
         DTOutput("sample_table"))
  ),
  nav_panel("Data", value = "data",
    card(card_header("Database"), verbatimTextOutput("db_info")),
    card(card_header("Exports"), tags$p("Published callsets are marked external_reference and are used only for overlays/validation; they are never AS3 model inputs or individual predictions. The SQLite database and normalized TSV files are computation artifacts."), downloadButton("download_segments", "Download filtered segments"))
  )
)
