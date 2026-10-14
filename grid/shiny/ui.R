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
  tabPanel('Sources / audit',value='audit',
    div(class='grid-card',h3('Available source tables'),tags$p('Choose a row to open that table. Historical metrics retain their original names.'),DTOutput('catalogue')),
    div(class='grid-card',h3('Import audit'),DTOutput('audit_table')),
    div(class='grid-card',h3('Snapshot'),verbatimTextOutput('metadata')))
)
