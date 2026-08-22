# ui.R
source("global.R", local = TRUE, encoding = "UTF-8")

## -------- I/O panel (minimal) --------
panel_setIO <- tabPanel(
  "I/O",
  fluidRow(
    column(
      width = 6,
      helpText("GWAS files are read from: ", strong(PATHS$GWAS_DIR)),
      helpText("You can also upload additional GWAS-VCF/txt files if you like (optional)."),
      fileInput("gwasFilecs", "Upload GWAS files (optional)", multiple = TRUE,
                buttonLabel = "Browse...", placeholder = PATHS$GWAS_DIR)
    )
  ),
  hr()
)

## -------- Run panel --------
panel_run <- tabPanel(
  fluidRow(
    column(
      width = 2,
      actionButton("analyze", "Run All", icon = icon("play")),
      align = "left"
    )
  )
)

## -------- Correlation panel --------
panel_Correlation <- tabPanel(
  "Correlation",
  sidebarLayout(
    sidebarPanel(
      h5("Pick GWAS for X vs Y, then click Run."),
      fluidRow(
        column(width = 6, uiOutput("files_C21")),
        column(width = 6, uiOutput("files_C22"))
      ),
      actionButton("runC2", "Run")
    ),
    mainPanel(uiOutput("report_C2"))
  )
)

## -------- Causation panel --------
panel_Causation <- tabPanel(
  "Causation",
  sidebarLayout(
    sidebarPanel(
      h5("Select exposure and outcome GWAS, then click Run."),
      fluidRow(
        column(width = 6, uiOutput("files_C31")),
        column(width = 6, uiOutput("files_C32"))
      ),
      textInput("n1", "Number of SNPs in exposure", "2000"),
      textInput("n2", "Number of SNPs in outcome",  "2000"),
      textInput("pvalue_cutoff", "p-value cutoff for instruments", "5e-10"),
      fluidRow(
        column(4, actionButton("defaultC3", "Default")),
        column(4, actionButton("resetC3",   "Reset")),
        column(4, actionButton("runC3",     "Run"))
      )
    ),
    mainPanel(uiOutput("report_C3"))
  )
)

## -------- Colocalization panel --------
panel_Colocalization <- tabPanel(
  "Colocalization",
  sidebarLayout(
    sidebarPanel(
      uiOutput("files_C41"),
      selectInput("gbuild", "Genome build", choices = c("hg19", "hg38"), selected = "hg19"),
      textAreaInput("gene", "Input: Gene list"),
      checkboxInput("GeneList", "multi_genes", TRUE),
      fluidRow(
        column(6, textInput("sigpvalue_eQTL", "sigpvalue_eQTL", "0.05")),
        column(6, textInput("sigpvalue_GWAS", "sigpvalue_GWAS", "5e-8"))
      ),
      actionButton("defaultC4", "Default"),
      actionButton("resetC4",   "Reset"),
      actionButton("runC4",     "Run")
    ),
    mainPanel(
      navbarPage(
        "",
        tabPanel("hyprcoloc", tableOutput("hyprcoloc")),
        tabPanel("heatmap", plotOutput("heatmap")),
        tabPanel("eQTLpLot", uiOutput("report_C41")),
        tabPanel("locuscomparer", uiOutput("report_C42"))
      )
    )
  )
)

## -------- Completion panel --------
panel_Completion <- tabPanel(
  "Completion",
  sidebarLayout(
    sidebarPanel(
      textAreaInput("genelist", "Input: Gene list"),
      selectInput("keyType", "ID Type", choices = c("SYMBOL", "ENSEMBL", "ENTREZID"), selected = "ENTREZID"),
      numericInput("pvalueCutoff", "p adjusted value cutoff", value = 0.05),
      numericInput("qvalueCutoff", "qvalue cutoff", value = 0.05),
      selectInput("pAdjustMethod", "p Adjust Method:",
                  choices = c("BH", "holm", "hochberg", "hommel", "bonferroni", "BY", "fdr", "none"),
                  selected = "BH"),
      checkboxInput("readable", "Convert gene id to gene symbol", TRUE),
      checkboxGroupInput("database5", "Annotation database",
                         choices = c("BP", "MF", "CC", "DO"), selected = "BP"),
      fluidRow(
        column(4, actionButton("defaultC5", "Default")),
        column(4, actionButton("resetC5",   "Reset")),
        column(4, actionButton("runC5",     "Run"))
      )
    ),
    mainPanel(uiOutput("report_C5"))
  )
)

## -------- Coregulation panel --------
panel_Coevolution <- tabPanel(
  "Coregulation",
  sidebarLayout(
    sidebarPanel(
      uiOutput("files_C6"),
      actionButton("runC6", "Run")
    ),
    mainPanel(uiOutput("report_C6"))
  )
)

## -------- Assemble UI --------
shinyUI(
  fluidPage(
    titlePanel("Progressively Establish Novel Genetics Underlying Interrelation"),
    navbarPage(
      title = NULL,
      panel_setIO,
      panel_Correlation,
      panel_Causation,
      panel_Colocalization,
      panel_Coevolution,
      panel_Completion,
      panel_run
    )
  )
)
