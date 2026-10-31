ui<-navbarPage('GRID · Twin evidence',id='nav',
  header=tags$head(tags$style(HTML('
    body{background:#f6f7f9} .container-fluid{max-width:1500px} .well{background:white;border:1px solid #e2e5e9;border-radius:8px}
    .grid-card{background:white;padding:18px 22px;margin:10px 0;border:1px solid #e2e5e9;border-radius:8px}
    .grid-source{font-size:12px;overflow-wrap:anywhere;color:#596574} h3{margin-top:5px;font-size:20px}
    .grid-controls{display:flex;gap:14px;flex-wrap:wrap;align-items:end} .grid-controls .form-group{min-width:160px}
  '))),
  tabPanel('Results',value='results',sidebarLayout(
    sidebarPanel(width=3,
      selectInput('trait','Trait',choices='all'),
      selectInput('module','Module',choices='all'),
      selectizeInput('table','Table / run',choices=NULL,options=list(placeholder='Choose a result table')),
      selectInput('ancestry','Ancestry',choices='ALL'),
      selectInput('subset','Evaluation subset',choices='ALL'),
      selectInput('source_metric','Recorded metric (long tables)',choices='ALL'),
      selectizeInput('methods','Methods',choices=NULL,multiple=TRUE),
      selectInput('metric','Plot metric (wide tables)',choices=NULL),
      tags$hr(),downloadButton('download_data','Filtered table'),tags$br(),tags$br(),downloadButton('download_plot','Current figure'),
      tags$hr(),tags$small('Confidence denotes reference support, not the probability that one individual’s outcome will be correct.')),
    mainPanel(width=9,
      div(class='grid-card',h3(textOutput('title')),textOutput('scope_note'),plotOutput('plot',height='480px')),
      div(class='grid-card',DTOutput('table_view'),textOutput('row_note'),tags$details(tags$summary('Source and interpretation'),div(class='grid-source',textOutput('source'))))
    ))),
  tabPanel('ABM',value='abm',
    div(class='grid-card',h3('ABM · count-based assessment'),
      tags$p('Fixed-baseline expected loss and absolute candidate benefit use count-only heads. Selection is frozen before assessment outcomes.'),
      div(class='grid-controls',selectInput('abm_trait','Trait',choices=NULL),
          selectInput('abm_run','Run / candidate',choices=NULL),
          selectInput('abm_view','View',choices=c('Baseline selected / control loss'='predictability_curve',
            'Gain enrichment'='benefit_curve','Full-query policy comparisons'='policy_performance',
            'Count additional information'='count_incremental_contrasts'))),
      plotOutput('abm_plot',height='500px'),
      tags$p('Tables retain nominal and actual coverage, evaluable N / groups / events and unavailable reasons. Count contrasts measure assessment-target MSE reduction, in squared assessment-loss units.'),
      DTOutput('abm_table'),downloadButton('abm_download_data','Filtered table'),
      downloadButton('abm_download_plot','Current figure'),textOutput('abm_source'))),
  tabPanel('Sources / audit',value='audit',
    div(class='grid-card',h3('Available source tables'),tags$p('Choose a row to open that table. Historical metrics retain their original names.'),DTOutput('catalogue')),
    div(class='grid-card',h3('Import audit'),DTOutput('audit_table')),
    div(class='grid-card',h3('Snapshot'),verbatimTextOutput('metadata')))
)
