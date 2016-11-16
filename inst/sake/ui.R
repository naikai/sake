library(shiny)
library(shinyBS)
library(plotly)
library(d3heatmap)
library(htmltools)
library(shinythemes)
library(shinydashboard)
library(DT)

tg <- tags$div()
tg <- attachDependencies(
  tg,
  htmlDependency(
    name    = "plotlyjs",
    version = "1.13.0",
    src     = c(href = "https://cdn.plot.ly"),
    script  = "plotly-latest.min.js"
  )
)

header <- dashboardHeader(
  title = "SAKE"
)

sidebar <- dashboardSidebar(
  sidebarMenu(id = "sidebarmenu",
    menuItem("File", tabName="File", icon=icon('file-text-o')),
    menuItem("Quality control", tabName="quality", icon = icon('bar-chart-o')),
    menuItem("Filtering", tabName="filter", icon = icon('scissors'),
             menuSubItem("Identify Variable Genes", tabName="featureSelection"),
             menuSubItem("Sample Scatter Plot", tabName="scatter")
    ),
    menuItem("NMF", tabName="NMF", icon = icon('sitemap'),
             menuSubItem("NMF Run", tabName="nmfRun"),
             menuSubItem("Identify Groups", tabName="nmfGroups"),
             menuSubItem("Extract features", tabName="nmfFeatures")
    ),
    menuItem("Visualization", tabName="Visualization", icon = icon('calendar', lib = "glyphicon")
    ),
             conditionalPanel("input.sidebarmenu === 'Visualization'",
                              selectInput("VisType",
                                          label = "What kind of plot?",
                                          choices = c("Heatmap", "PCA", "t-SNE", "Cluster Distribution"),
                                          selected = "Cluster Distribution")
             ),
    menuItem("Differential analysis", tabName="DE", icon = icon('tasks', lib="glyphicon"),
             menuSubItem("DESeq2", tabName="DESeq2")
             # menuSubItem("SCDE", tabName="SCDE")
    ),
    menuItem("Enrichment analysis", tabName="enrichment", icon = icon('cogs')
    ),
             conditionalPanel("input.sidebarmenu === 'enrichment'",
                              selectInput("EnrichType",
                                          label = "Which Category?",
                                          choices = c("GO", "KEGG"),
                                          selected = "KEGG")
             )
  )
)

body <- dashboardBody(
  tags$head(tags$script(src = "google-analytics.js")),
  tags$head(tags$style(".dwnld{background-color:#8E24BF;} .dwnld{color: white;} .dwnld{border-color:$9932CC;}")),
  tags$head(tags$style(".act{background-color:#337AB7;} .act{color: #FFFFFF;} .act{border-color:$2E6DA4;}")),
  tags$head(tags$style(".stop{background-color:#4F1E0A;} .stop{color: #FFFFFF;} .stop{border-color:#612C0E;}")),
  tags$head(tg),

  tabItems(
    tabItem("File",
            fluidRow(
              tabBox( width=5, id = "tabset1", height = "575px",
                      tabPanel("Upload file",
                               selectInput("selectfile",
                                           label = "How do you want to upload the data?",
                                           choices = list("Select from preloaded data"='preload',
                                                          "Upload rawdata"='upload',
                                                          "Upload saved NMF run"='saved'),
                                           selected = "saved"),
                               conditionalPanel(
                                 condition = "input.selectfile == 'upload'",
                                 fileInput('file1', 'Choose text File',
                                           accept=c('text/csv',
                                                    'text/comma-separated-values',
                                                    'text/tab-separated-values',
                                                    'text/plain',
                                                    '.csv',
                                                    '.tsv',
                                                    '.out'))
                               ),
                               conditionalPanel(
                                 condition = "input.selectfile == 'preload'",
                                 selectizeInput(inputId = 'inputdata',
                                                label = 'Choose DataSet',
                                                choices = filenames,
                                                options = list(
                                                  placeholder = 'Please select a data set',
                                                  onInitialize = I('function() { this.setValue(""); }')
                                                )
                                 )
                               ),
                               conditionalPanel(
                                 condition = "input.selectfile == 'upload' || input.selectfile == 'preload'",
                                 tags$hr(),
                                 checkboxInput('header', 'Header', TRUE),
                                 radioButtons('sep', 'Separator',
                                              c(Comma=',',
                                                Semicolon=';',
                                                Tab='\t'),
                                              selected='\t')
                               ),
                               conditionalPanel(
                                 condition = "input.selectfile == 'saved'",
                                 fileInput('rda', 'Choose saved RData',
                                           accept=c('.rda',
                                                    '.RData'))
                               )
                      ),
                      tabPanel("System setup",
                               selectInput("ncores", label = "Specify number of cores to run analysis",
                                           choices = c(1,2,4,8,16),
                                           selected = 16)
                      )
              ),
              column(width = 7,
                     box(width = NULL, title="Preview Raw Data", status = "warning",
                         DT::dataTableOutput('rawdatatbl', height = "439px"),
                         tags$hr(),
                         p(
                           class = "text-muted",
                           paste("Note: You can preview the loaded data, ",
                                 "make sure the format is what you want to then proceed."
                           )
                         )
                     )
              )
            ),
            fluidRow(
              infoBoxOutput("nameBox"),
              infoBoxOutput("sampleBox"),
              infoBoxOutput("featureBox")
            )
    ),
    tabItem("quality",
            fluidRow(
              box(title="Datatable", width=12, solidHeader=TRUE, status="success",
                  fluidRow(
                    column(width=3,
                           selectInput("normdata", label = "Normalization",
                                         choices = list("RPM normalization" = "takeRPM",
                                                      "DESeq size factor normalization" = "takesizeNorm",
                                                      "Upper Quartile normalization" = "takeuq",
                                                      "None" = "none"),
                                       selected = "none")
                    ),
                    column(width=3,
                           selectInput("transformdata", label = "Transformation",
                                       choices = list("Variance Stablizing Transformation" = "takevst",
                                                      "Log transformation" = "takelog",
                                                      "None" = "none"),
                                       selected = "none")
                    ),
                    column(width=2, br(),
                           actionButton("run_transf", "Run Transform", icon("play-circle"), class = 'act')
                    ),
                    column(width=2, br(),
                           downloadButton("dl_transformdata", "Download Data", class="dwnld")
                    )
                  )
              )
            ),
            fluidRow(
              column(width=6,
                     box(title="Read Distribution", width=NULL, solidHeader=TRUE, status="info", height="570px",
                         plotlyOutput("readDistrib", height=500)
                     )
              ),
              column(width=6,
                box(title="Samples", width=NULL, solidHeader=TRUE, status="info", height="350px",
                    DT::dataTableOutput("selsamp_tb")
                ),
                box(title="Filter", width=NULL, solidHeader=TRUE, status="info", height="200px",
                    verbatimTextOutput("selsamp_txt"),
                    actionButton("run_filtsamp", "Submit", icon("play-circle"), class = 'act')
                )
              )
            )
    ),
    tabItem("featureSelection",
            fluidRow(
              box(width=12, collapsible = TRUE,
                  fluidRow(
                    column(width=2,
                           selectInput("list_type",
                                       label = "How to extract genes",
                                       choices = c("Whole transcriptome",
                                                   "Top Ranks",
                                                   "Upload gene list"),
                                       selected = "Top Ranks"
                           )
                    ),
                    conditionalPanel(
                      condition = "input.list_type == 'Top Ranks'",
                      column(width=2, selectInput("orders",
                                                  label = "Rank order by",
                                                  choices = c("top", "bottom"),
                                                  selected = "top")),
                      column(width=2, selectInput("math",
                                                  label = "Arithmetic",
                                                  choices = c("mean", "median", "mad", "var"),
                                                  selected = "mad")),
                      column(width=4, sliderInput("top.num",
                                                  label = "How many genes to select", ticks = FALSE,
                                                  min=500, max=20000, value=5000, step = 500)
                      )
                    ),
                    conditionalPanel(
                      condition = "input.list_type == 'Upload gene list'",
                      column(width=3,
                             fileInput('genefile', 'Choose gene list file',
                                       accept=c('text/csv',
                                                'text/comma-separated-values,text/plain',
                                                '.csv'))
                      )
                    )
                  )
              )
            ),
            fluidRow(
              column(width=6,
                     featureUI("sample", title = "Sample Correlation")
              ),
              column(width=6,
                     networkUI("gene", title = "Gene Network")
              )
            )
    ),
    tabItem("scatter",
            fluidRow(
              box(title="Parameters", width=12, solidHeader = TRUE, status="success",
                  fluidRow(
                    column(width=2, selectizeInput("cor_sample1", label = "Sample on X-axis", choices=NULL, multiple=FALSE)),
                    column(width=2, selectizeInput("cor_sample2", label = "Sample on Y-axis", choices=NULL, multiple=FALSE)),
                    column(width=2, checkboxInput("scatter_log", label = "Log transformation", TRUE)),
                    column(width=2, checkboxInput("show_r2", label = "Show R2", FALSE)),
                    column(width=2, checkboxInput("show_fc", label = "Show Fold-Change line", TRUE))
                  )
              ),
              box(width=6, height="520px", plotlyOutput('scatterPlot')),
              box(width=6, height="520px", DT::dataTableOutput('scatterdatatbl'))
            )
    ),
    tabItem("nmfRun",
            fluidRow(
              box(title="Estimate Number of Clusters", width=12, solidHeader=TRUE, status="success",
                  fluidRow(
                    column(width=2,
                           selectInput("mode",
                                       label = "Which module?",
                                       choices = c("Estimate number of k" = "estim",
                                                   "Real run" = "real"),
                                       selected = "real"),
                           bsTooltip("mode", "We recommend you run estimate number of cluster (k) first <br> Then use the best estimate cluster (k) for real run",
                                     "right", options = list(container = "body"))
                    ),
                    column(width=2,
                           numericInput("num_cluster",
                                        label = "Estimate clusters",
                                        min=2, max=10, value=3, step=1),
                           bsTooltip("num_cluster", "For Estimate number of k, it will run NMF from 2 to k clusters <br><br> For Real run, it will run NMF for that specific k cluster",
                                     "right", options = list(container = "body"))
                    ),
                    column(width=2,
                           numericInput("nrun",
                                        label = "Num of runs",
                                        min=10, max=200, value=7, step=10),
                           bsTooltip("nrun", "Recommend 20-30 for estimate k <br> 50-100 for real run",
                                     "right", options = list(container = "body"))
                    ),
                    column(width=2,
                           actionButton("runNMF", "Run NMF", icon("play-circle"), class = 'act'),
                           checkboxInput('nmf_runmoreopt', 'More Options', FALSE)
                    ),
                    uiOutput("dlNMF_UI"),
                    column(width=6, bsAlert("NMFAlert")),
                    column(width=6, bsAlert("YabiAlert")),
                    bsModal("nmfModal1", "Warning!\nYour gene size is above 1000, it will take longer than usual.\nAre you sure you want to continue?",
                            trigger = "", size = "small",
                            textInput("mailto", label = "To:", placeholder = "To:"),
                            actionButton("nmf_yabi_yes", 'Yes'),
                            actionButton("nmf_yabi_no", 'No')
                    )
                  ),
                  fluidRow(
                    conditionalPanel(
                      condition = "input.nmf_runmoreopt == true",
                      column(width=2,
                             selectInput("algorithm",
                                         label = "NMF algorithm",
                                         choices = c("brunet", "lee", "nsNMF", "KL", "Frobenius"),
                                         selected = "brunet"),
                             bsTooltip("algorithm", "Different NMF algorithms",
                                       "right", options = list(container = "body"))
                      ),
                      column(width=2,
                             selectInput("nmf_seed",
                                         label = "Random seed",
                                         choices = c("Yes"=0, "No"=123211),
                                         selected = 123211),
                             bsTooltip("nmf_seed", "Whether to use random starting condition for NMF run <br> Default: fixed seed",
                                       "right", options = list(container = "body"))
                      )
                    )
                  )
              )
            ),
            conditionalPanel(
              condition = "input.mode == 'estim'",
              fluidRow(
                box(title="NMF Plot", width=8, solidHeader = TRUE, status="info", height="1000px",
                    collapsible = TRUE,
                    fluidRow(
                      column(width=4, selectInput("estimtype",
                                                  label = "What kind of plot?",
                                                  choices = c("Quality measures" = "quality",
                                                              "Consensus matrics" = "consensus"),
                                                  selected = "consensus")
                      ),
                      column(width=3, br(),
                             downloadButton("dl_nmf_estimplot", "Save plot", class="dwnld")
                      ),
                      fluidRow(
                        column(width=4,
                               br(),
                               checkboxInput('dl_nmf_estimmoreopt', 'More Options', FALSE),
                               conditionalPanel(
                                 condition = "input.dl_nmf_estimmoreopt == true",
                                 strong("PDF size (inches):"),
                                 sliderInput(inputId="estim_pdf_w", label = "width:", min=3, max=20, value=14, ticks=F),
                                 sliderInput(inputId="estim_pdf_h", label = "height:", min=3, max=20, value=14, ticks=F)
                               )
                        )
                      )
                    ),
                    fluidRow(
                      column(width=12,
                             plotOutput('estimPlot')
                      )
                    )
                ),
                box(title="Stats", width=4, solidHeader = TRUE, status="warning", height="1000px",
                    DT::dataTableOutput('estimSummary')
                )
              )
            ),
            conditionalPanel(
              condition = "input.mode == 'real'",
              fluidRow(
                box(title="NMF Plot", width=8, solidHeader = TRUE, status="info", height = "1000px",
                    fluidRow(
                      column(width=4, selectInput("plottype",
                                                  label = "What kind of plot?",
                                                  choices = c("samples", "features", "consensus"),
                                                  selected = "samples")
                      ),
                      column(width=3, br(),
                             downloadButton("dl_nmf_realplot", "Save plot", class="dwnld")
                      ),
                      fluidRow(
                        column(width=4,
                               br(),
                               checkboxInput('dl_nmf_realmoreopt', 'More Options', FALSE),
                               conditionalPanel(
                                 condition = "input.dl_nmf_realmoreopt == true",
                                 checkboxInput('nmfplot_silhouette', 'Match silhouette order', TRUE),
                                 strong("PDF size (inches):"),
                                 sliderInput(inputId="real_pdf_w", label = "width:", min=3, max=20, value=14, ticks=F),
                                 sliderInput(inputId="real_pdf_h", label = "height:", min=3, max=20, value=14, ticks=F)
                               )
                        )
                      )
                    ),
                    fluidRow(
                      column(width=12, plotOutput('nmfplot'))
                    )
                ),
                box(title="Silhouette Plot", width=4, solidHeader=TRUE, status="info", height="1000px",
                    plotOutput('silhouetteplot')
                )
              )
            )
    ),
    tabItem("nmfGroups",
            fluidRow(
              box(title="Identify Groups", width=4, solidHeader=TRUE, status="success", height = "740px",
                  fluidRow(
                    column(width=12,
                           selectInput("predict",
                                label = "How to assign NMF groups?",
                                choices = c("consensus", "samples"),
                                selected = "samples"),
                           checkboxInput('nmf_prob', 'Change dot size based on probability?', FALSE),
                           bsTooltip("nmf_prob", "Samples with higher uncertainty (lower probability of being correctly assigned) will become larger dots",
                                     "right", options = list(container = "body")),
                           DT::dataTableOutput('nmfGroups')
                    )
                  )
              ),
              tabBox(title=tagList(shiny::icon("th-large"), "t-SNE Plot"),
                     width=8, id = "tabset2", side = "right", height = "740px",
                  tabPanel("2D",
                           column(width=2, numericInput("nmftsne_perplexity", label = "Perplexity", value=10)),
                           column(width=2, br(),
                                  actionButton("run_nmftSNE", "Run t-SNE!", icon("play-circle"), class = "act")
                           ),
                           uiOutput("dlperplexity_UI"),
                           column(width=6, bsAlert("perplexityAlert")),
                           column(width=12,
                                  plotlyOutput('nmftsneplot', height=600, width=600)
                           )
                  )
              )
            )
    ),
    tabItem("nmfFeatures",
            fluidRow(
              box(title="Extract Features", width=12, solidHeader=TRUE, status="success",
                  fluidRow(
                    column(width=3,
                           selectInput("select_method",
                                       label = "Feature selection method",
                                       choices = c("Total"="total", "By default"="default", "rank"="rank"),
                                       selected = "default")
                    ),
                    conditionalPanel(
                      condition = "input.select_method == 'default' | input.select_method == 'rank'",
                      column(width=3,
                             selectInput("select_feature_num",
                                         label = "Number of features from each group",
                                         choices = c("By default"=0, "10"=10, "20"=20, "30"=30, "50"=50, "100"=100, "200"=200, "500"=500, "1000"=1000, "1500"=1500, "2000"=2000, "3000"=3000, "5000"=5000),
                                         selected = 0)
                      ),
                      conditionalPanel(
                        condition = "input.select_method == 'rank'",
                        column(width=3,
                               selectInput("select_math",
                                           label = "Rank based on:",
                                           choices = c("mean", "median", "mad", "var"),
                                           selected = "mad")
                        ),
                        column(width=3,
                               numericInput("select_FScutoff",
                                            label = "FeatureScore cutoff",
                                            min=0.5, max=1, value=0.86, step = 0.01)
                        )
                      )
                    )
                  ),
                  fluidRow(
                    column(width=8, DT::dataTableOutput('nmfFeatures')),
                    column(width=4, plotlyOutput('nmf_boxplot', height = "500px"))
                  )
              )
            )
    ),
    tabItem("Visualization",
            fluidRow(
              box(width=12,
                  fluidRow(
                    column(width=3, selectInput("heatmapGene",
                                                label = "How to select genes?",
                                                choices = c("From NMF features" = "nmf",
                                                            "Ranks from data" = "rank",
                                                            "Upload gene list" = "upload",
                                                            "Manually select genes" = "manual",
                                                            "Preloaded gene list" = "preload"),
                                                selected = "rank")
                    ),
                    conditionalPanel(
                      condition = "input.heatmapGene == 'upload'",
                      column(width=3, fileInput('heatmapGeneFile', 'Choose text File',
                                                accept=c('text/csv',
                                                         'text/comma-separated-values,text/plain',
                                                         '.csv',
                                                         '.out'))
                      )
                    ),
                    conditionalPanel(
                      condition = "input.heatmapGene == 'rank'",
                      column(width=3, selectInput("heat.math",
                                                  label = "Rank by which metrics?",
                                                  choices = c("Mean Expression" = "mean",
                                                              "Median Expression" = "median",
                                                              "Median Absolute Deviation" = "mad",
                                                              "Variance" = "var"),
                                                  selected = "mad")
                      ),
                      column(width=2, selectInput("heat.top.num",
                                                  label = "How many top genes:",
                                                  choices = c(seq(10,100,10), seq(200, 1000, 100), seq(2000, 15000, 1000)),
                                                  selected = 5000)
                      )
                    ),
                    conditionalPanel(
                      condition = "input.heatmapGene == 'preload'",
                      column(width=3,
                             selectizeInput("predefined_list",
                                            label = "Select from predefined lists",
                                            choices = list(
                                              Melanoma = c(`Keith Hoek (105 genes, 2006, Pigment Cell & Melanoma)`='Hoek',
                                                           `Widmer (Hoek's group) (97 genes, 2012, Pigment Cell & Melanoma)`='Widmer',
                                                           `Stein Aert (core 59 genes, 2015, Nature Communication)`='Stein_core',
                                                           `Stein Aert (whole 1412 genes, 2015, Nature Communication)`='Stein',
                                                           `SKCM NMF-K3 (53 genes, 2015, TCGA)`='SKCM-NMF-K3',
                                                           `SKCM NMF-K3-VST (101 genes, 2015, TCGA)`='SKCM-NMF-K3-VST-1000nrun',
                                                           `SKCM NMF-K4 (61 genes, 2015, TCGA)`='SKCM-NMF-K4',
                                                           `Stephanie Kreis (20 miRNAs, 2015, Oncotarget)`='Melanoma-20miRNA'
                                              ),
                                              Hallmark = c(`Cell Cycle phases`='CellCycle',
                                                           `Regulation of Cell cycle`='RegCC',
                                                           `EMT`='EMT',
                                                           `KRAS`='KRAS',
                                                           `WNT`='WNT'
                                              ),
                                              Others = c(`Resistance markers (combined from literature search)`='PubResistList',
                                                         `Resistance markers (Only select core)`='PubResistList_core',
                                                         `My 4 groups Gene List (May modify later)`='MyFourGroups',
                                                         `Aviv Regev cell cycle genes`='Aviv'
                                              ),
                                              Mouse =  c(`Cell Cycle phases`='Mouse_CellCycle',
                                                         `Regulation of Cell cycle`='Mouse_RegCellCycle'
                                              ),
                                              Breast = c(`PAM50`='PAM50',
                                                         `BRCA-NMF-K5`="BRCA-NMF-K5"
                                              )
                                            ),
                                            selected = 'EMT', multiple=TRUE
                             )
                      )
                    ),
                    conditionalPanel(
                      condition = "input.heatmapGene == 'manual'",
                      column(width=3, selectizeInput("manual_genes",
                                                     label = "Select your favorite genes",
                                                     choices=NULL, multiple=TRUE,
                                                     options = list(placeholder='Type in gene name, e.g. MYC, ERBB2')
                      )
                      )
                    ),

                    conditionalPanel(
                      condition = "input.VisType == 'Heatmap' ",
                      column(width = 3, br(),
                             downloadButton("dl_heatmap", "Download", class = "butt"),
                             tags$head(tags$style(".butt{background-color:#8E24BF;} .butt{color: white;} .butt{border-color:$9932CC;}"))
                      )
                    ),

                    column(width=2, checkboxInput('heatmoreopt', 'More Options', FALSE)),
                    conditionalPanel(
                      condition = "input.VisType == 't-SNE' ",
                      column(width=2, actionButton("runtSNE", "Run t-SNE!", icon("play-circle"), width="100px", class = 'act'))
                    )
                  ),
                  conditionalPanel(
                    condition = "input.heatmoreopt == true & input.VisType == 'Heatmap' ",
                    fluidRow(
                      column(width=6,
                             selectInput("heat_opttype",
                                         label = "Which paramter type",
                                         choices = list('Figure adjustment' = c(`Column Side Colors` = 'ColSideColors',
                                                                                `Row Side Colors` = 'RowSideColors',
                                                                                `Heatmap Color Scheme` = 'heatcolor',
                                                                                `PDF Size` = 'pdf',
                                                                                `Label Size and Apperance` = 'labels'),
                                                        'Clustering' = c(`Cluster Row or Column` = 'cluster',
                                                                         `How to Scale Data` = 'scale',
                                                                         `Distance and Linkage` = 'distnlink')),
                                         selected = "cluster")
                      )
                    ),
                    conditionalPanel(
                      condition = "input.heat_opttype == 'cluster'",
                      fluidRow(
                        column(width=3,
                               selectizeInput("OrdCol",
                                              label = "Cluster column?",
                                              choices=NULL, multiple=FALSE,
                                              options = list(placeholder='Hierarchical')),
                               conditionalPanel(
                                 condition = "input.OrdCol == 'Filename'",
                                 numericInput("sortcolumn_num", label = "Sort by which group:", value = 1, min = 1)
                               ),
                               conditionalPanel(
                                 condition = "input.OrdCol == 'GeneExpr'",
                                 selectizeInput("sortgene_by",
                                                label = "Select a gene to sort on",
                                                choices=NULL, multiple=FALSE,
                                                options = list(placeholder='Type in gene name, e.g. MYC, ERBB2'))
                               )
                        ),
                        column(width=3, selectInput("OrdRow", "Cluster row?",
                                                    choices = list("Default ordering" = "default",
                                                                   "Hierarchical" = "hier"),
                                                    selected = "hier"))
                      )


                    ),
                    conditionalPanel(
                      condition = "input.heat_opttype == 'scale'",
                      fluidRow(
                        column(width=3, selectInput("scale_method", label = "How to scale data:", choices = c("mean", "median"), selected = "mean")),
                        column(width=4, checkboxInput("scale_first", "Scale data before dendrogram?", FALSE))
                      )
                    ),
                    conditionalPanel(
                      condition = "input.heat_opttype == 'distnlink'",
                      fluidRow(
                        column(width=3, selectInput("linkage",
                                                    label = "Which linkage methods?",
                                                    choices = c("average", "centroid", "complete","mcquitty", "median", "ward.D", "ward.D2", "single"),
                                                    selected = "ward.D2")),
                        column(width=3, selectInput("distance",
                                                    label = "Which distance methods?",
                                                    choices = c("binary", "canberra","euclidean", "manhattan", "maximum", "minkowski"),
                                                    selected = "euclidean"))
                      )
                    ),
                    conditionalPanel(
                      condition = "input.heat_opttype == 'labels'",
                      fluidRow(
                        column(width=2, textInput("title", label= "Figure title", value=NULL)),
                        column(width=2, numericInput("cexRow", label = "Row LabelSize:", value=0.5, step=0.1)),
                        column(width=2, numericInput("cexCol", label = "Col LabelSize:", value=0.4, step=0.1)),
                        column(width=3, selectInput('row_legend', label='Show legend for RowSide Colors?',
                                                    choices = c("Yes" = TRUE, "No" = FALSE))),
                        column(width=3, selectInput('col_legend', label='Show legend for ColumnSide Colors?',
                                                      choices = c("Yes" = TRUE, "No" = FALSE)))
                      )
                    ),
                    conditionalPanel(
                      condition = "input.heat_opttype == 'heatcolor'",
                      fluidRow(
                        column(width=3, selectInput("heatColorScheme",
                                                    label = "Heat map color scheme",
                                                    choices = c("RdBu_Brewer", "RdYlBu_Brewer", "Red_Green", "Red_Blue","Rainbow", "Black_White"),
                                                    selected = "RdBu_Brewer")),
                        column(width=3, numericInput("color_num", label = "Num", value = 9, min=5))
                      )
                    ),
                    conditionalPanel(
                      condition = "input.heat_opttype == 'pdf'",
                      fluidRow(
                        column(width=3, sliderInput("heat_pdf_wd", label = "Width (inch)", value=14, min=5, max=20, ticks=F)),
                        column(width=3, sliderInput("heat_pdf_ht", label = "Height (inch)", value=12, min=5, max=20, ticks=F))
                      )
                    ),
                    conditionalPanel(
                      condition = "input.heat_opttype == 'RowSideColors'",
                      fluidRow(
                        column(width=3, selectInput("RowColScheme",
                                                    label = "Row Color Scheme:",
                                                    choices = c("naikai", "naikai2", "Set1", "Set2", "Set3","YlGn", "YlGnBu", "YlOrRd", "Greys", "Pastel1", "Pastel2", "Paired", "Dark2"),
                                                    selected = "naikai2")),
                        column(width=3, numericInput("RowColScheme.num", label = "Num of colors:", value = 7))
                      )
                    ),
                    conditionalPanel(
                      condition = "input.heat_opttype == 'ColSideColors'",
                      fluidRow(
                        column(width=2, selectizeInput("ColClrBy",
                                       label = "Color sample by?",
                                       choices=NULL, multiple=FALSE,
                                       options = list(placeholder='Hierarchical'))
                        ),
                        column(width=2, numericInput("ColSideColorsSize", label = "Column Color Size:", value=2.0, step=0.1)),
                        conditionalPanel(
                          condition = "input.ColClrBy == 'Filename'",
                          column(width=2, selectInput("SelColSideColors",
                                                      label = "Show columns in order:",
                                                      multiple = TRUE,
                                                      choices = c(1,2,3,4,5),
                                                      selected = 1)),
                          column(width=2, sliderInput("ColSideColorsNum", label = "Modify Which Column:", min=1, max=5, value=1, ticks = F)),
                          # Need to fix this later
                          conditionalPanel(
                            condition = "input.ColSideColorsNum == 1",
                            column(width=2, selectInput("ColScheme1",
                                                        label = "Column Color 1:",
                                                        choices = c("naikai", "naikai2", "Set1", "Set2", "Set3","YlGn", "YlGnBu", "YlOrRd", "Greys", "Pastel1", "Pastel2", "Paired", "Dark2"),
                                                        selected = "naikai")),
                            column(width=2, numericInput("ColScheme1.num", label = "Num of colors:", value = 5))
                          ),
                          conditionalPanel(
                            condition = "input.ColSideColorsNum == 2",
                            column(width=2, selectInput("ColScheme2",
                                                        label = "Column Color 2:",
                                                        choices = c("naikai", "naikai2", "Set1", "Set2", "Set3","YlGn", "YlGnBu", "YlOrRd", "Greys", "Pastel1", "Pastel2", "Paired", "Dark2"),
                                                        selected = "Paired")),
                            column(width=2, numericInput("ColScheme2.num", label = "Num of colors:", value = 10))
                          ),
                          conditionalPanel(
                            condition = "input.ColSideColorsNum == 3",
                            column(width=2, selectInput("ColScheme3",
                                                        label = "Column Color 3:",
                                                        choices = c("naikai", "naikai2", "Set1", "Set2", "Set3","YlGn", "YlGnBu", "YlOrRd", "Greys", "Pastel1", "Pastel2", "Paired", "Dark2"),
                                                        selected = "Greys")),
                            column(width=2, numericInput("ColScheme3.num", label = "Num of colors:", value = 10))
                          ),
                          conditionalPanel(
                            condition = "input.ColSideColorsNum == 4",
                            column(width=2, selectInput("ColScheme4",
                                                        label = "Column Color 4:",
                                                        choices = c("naikai", "naikai2", "Set1", "Set2", "Set3","YlGn", "YlGnBu", "YlOrRd", "Greys", "Pastel1", "Pastel2", "Paired", "Dark2"),
                                                        selected = "YlOrRd")),
                            column(width=2, numericInput("ColScheme4.num", label = "Num of colors:", value = 10))
                          ),
                          conditionalPanel(
                            condition = "input.ColSideColorsNum == 5",
                            column(width=2, selectInput("ColScheme5",
                                                        label = "Column Color 5:",
                                                        choices = c("naikai", "naikai2", "Set1", "Set2", "Set3","YlGn", "YlGnBu", "YlOrRd", "Greys", "Pastel1", "Pastel2", "Paired", "Dark2"),
                                                        selected = "Pastel1")),
                            column(width=2, numericInput("ColScheme5.num", label = "Num of colors:", value = 10))
                          )
                        )
                      )
                    )
                  ),
                  conditionalPanel(
                    condition = "input.heatmoreopt == true & (input.VisType=='PCA' || input.VisType=='t-SNE')",
                    fluidRow(
                      column(width=2, selectizeInput("pt_col",
                                                     label = "Color sample by?",
                                                     choices=NULL, multiple=FALSE,
                                                     options = list(
                                                       onInitialize = I('function() { this.setValue(""); }')
                                                     ))
                      ),
                      conditionalPanel(
                        condition = "input.pt_col == 'GeneExpr'",
                        column(width=2, selectizeInput("pt_allgene",
                                                       label = "Select gene",
                                                       choices=NULL, multiple=FALSE,
                                                       options = list(
                                                         onInitialize = I('function() { this.setValue(""); }')
                                                       ))
                        )
                      ),
                      conditionalPanel(
                        condition = "input.pt_col == 'NMF Feature'",
                        column(width=2, selectizeInput("pt_nmfgene",
                                                       label = "Select gene",
                                                       choices=NULL, multiple=FALSE,
                                                       options = list(
                                                         onInitialize = I('function() { this.setValue(""); }')
                                                       ))
                        )
                      ),
                      conditionalPanel(
                        condition = "input.pt_col == 'Filename'",
                        column(width=2, selectizeInput("pt_file_grp",
                                                       label="Which Group",
                                                       choices=NULL, multiple=TRUE,
                                                       options = list(
                                                         onInitialize = I('function() { this.setValue(""); }')
                                                       ))
                        )
                      ),
                      conditionalPanel(
                        condition = "input.VisType == 'PCA'",
                        column(width=2, selectInput("pca_x", label="PC on X-axis", choices=seq(1,10), selected = 1)),
                        column(width=2, selectInput("pca_y", label="PC on Y-axis", choices=seq(1,10), selected = 2))
                      ),
                      conditionalPanel(
                        condition = "input.VisType == 't-SNE'",
                        column(width=1, numericInput("tsne_perplexity", label = "Perplexity", value=10)),
                        column(width=1, numericInput("tsne_iter", label = "Iterations", value=20))
                      ),
                      column(width=1, numericInput("plot_point_size", label = "DotSize", value=7)),
                      column(width=1, numericInput("plot_point_alpha", label = "Alpha", value=0.75, step=0.05)),
                      column(width=1, numericInput("plot_label_size", label = "LabelSize", value=9)),
                      column(width=1, checkboxInput('plot_label', 'Add label', FALSE)),
                      column(width=1, checkboxInput('plot_legend', 'Add legend', TRUE)),
                      column(width=2, selectInput("pt_sel_grp_order", label="Group Order", choices=c("Default", "Manual"), selected="Default")),
                      conditionalPanel(
                        condition = "input.pt_sel_grp_order == 'Manual' & (input.pt_col == 'Filename' || input.pt_col == 'NMF Group')",
                        column(width=2, selectizeInput("pt_grp_order",
                                                       label="Select group order",
                                                       choices=NULL, multiple=TRUE,
                                                       options = list(
                                                         onInitialize = I('function() { this.setValue(""); }')
                                                       ))
                        )
                      ),
                      uiOutput("visualperplexity_UI"),
                      column(width=6, bsAlert("visualperplexityAlert"))
                    )
                  ),
                  conditionalPanel(
                    condition = "input.heatmoreopt == true & input.VisType=='Cluster Distribution'",
                    fluidRow(
                      column(width=2, numericInput("cl_alpha", label = "Alpha", value=0.4, step=0.05, max=1, min=0.05)),
                      column(width=2, numericInput("cl_lgsize", label = "Legend size", value=10, step=1, max=14, min=6)),
                      column(width=2, numericInput("cl_tickfont", label = "Tick font size", value=10, step=1, max=20, min=4)),
                      column(width=2, numericInput("cl_axisfont", label = "Axis font size", value=14, step=1, max=20, min=4))
                    )
                  )
              ),
              conditionalPanel(
                condition = "input.VisType == 'Heatmap'",
                box(width=12, d3heatmapOutput("d3heatmap", height=700))
              ),
              conditionalPanel(
                condition = "input.VisType == 'PCA'",
                box(width=6, title="2D PCA Plot", solidHeader=TRUE, status="info", collapsible = TRUE, plotlyOutput("pcaPlot_2D", height=500)),
                box(width=6, title="3D PCA Plot", solidHeader=TRUE, status="info", collapsible = TRUE, plotlyOutput("pcaPlot_3D", height=500))
              ),
              conditionalPanel(
                condition = "input.VisType == 't-SNE'",
                box(width=6, title="2D t-SNE Plot", solidHeader=TRUE, status="info", collapsible = TRUE, plotlyOutput('tsneplot_2d', height=500)),
                box(width=6, title="3D t-SNE Plot", solidHeader=TRUE, status="info", collapsible = TRUE, plotlyOutput("tsneplot_3d", height=500))
              ),
              conditionalPanel(
                condition = "input.VisType == 'Cluster Distribution'",
                box(width=4, title="Transcriptome variance", solidHeader=TRUE, status="info", collapsible = TRUE, plotlyOutput("nmfgrp_var", height=350)),
                box(width=4, title="No. of expressed genes", solidHeader=TRUE, status="info", collapsible = TRUE, plotlyOutput('nmfgrp_expgene', height=350)),
                box(width=4, title="Intra-group correlation coefficient", solidHeader=TRUE, status="info", collapsible = TRUE, plotlyOutput('nmfgrp_coef', height=350))
              )
            )
    ),
    tabItem("DESeq2",
            fluidRow(
              box(title="Differential Analysis Parameters", width=12, solidHeader=TRUE, status="success",
                  fluidRow(
                    column(width=2, selectizeInput("de_group1",
                                                   label = "Reference group",
                                                   choices=NULL, multiple=FALSE,
                                                   options = list(placeholder='Choose a NMF group'))
                    ),
                    column(width=2, selectizeInput("de_group2",
                                                   label = "Treatment group",
                                                   choices=NULL, multiple=FALSE,
                                                   options = list(placeholder='Choose a NMF group'))
                    ),
                    column(width=2, numericInput("de_alpha", label = "alpha cutoff", min=0.01, max=0.1, value=0.05, step = 0.01)),
                    column(width=2, selectInput("de_conv2int", label = "Force convert to integer",
                                                choices = c("Yes" = TRUE, "No" = FALSE),
                                                selected = FALSE)
                    ),
                    column(width=2, br(),
                           actionButton("runDESeq", "Run DESeq!", icon("play-circle"), class = 'act')
                    )
                  ),
                  fluidRow(
                    column(width=8, DT::dataTableOutput('deseq_table')),
                    column(width=4, plotlyOutput('deseq_boxplot', height = "500px"))
                  ),
                  fluidRow(
                    column(width=12, bsAlert("alert"))
                  )
              )
            )
    ),
    tabItem("enrichment",
            fluidRow(
              box(title="Enrichment Analysis Parameters", width=12, solidHeader=TRUE, status="success",
                  fluidRow(
                    column(width=2,
                           selectizeInput("pathway_group",
                                          label = "Pick one NMF group",
                                          choices=NULL, multiple=FALSE,
                                          options = list(placeholder='Choose a NMF group')),
                           bsTooltip("pathway_group", "Choose one NMF group to perform pathway enrichment test",
                                     "right", options = list(container = "body"))
                    ),
                    column(width=2,
                           selectInput("species",
                                       label = "Species",
                                       choices = c("Human", "Mouse"),
                                       selected = "Human"),
                           bsTooltip("species", "Default set to Human <br>Please select the corresponding species for your data, the program will run the conversion if necessary",
                                     "right", options = list(container = "body"))
                    ),
                    conditionalPanel(
                      condition = "input.EnrichType == 'GO'",
                      column(width=2, selectInput("goTerm",
                                                  label = "GO category",
                                                  choices = c("BP", "MF", "CC"),
                                                  selected = "MF")
                      )
                    ),
                    column(width=2, br(),
                           actionButton("run_enrich", "Run Analysis", icon("play-circle"), class = 'act')
                    ),
                    column(width=2, checkboxInput('enrichmoreopt', 'More Options', FALSE))
                  ),
                  conditionalPanel(
                    condition = "input.enrichmoreopt == true",
                    fluidRow(
                      column(width=2, selectInput("rank_method",
                                                  label = "Rank genes by?",
                                                  choices = c("NMF FeatureScore"="featureScore",
                                                              "Log Fold Change" = "logFC"),
                                                  selected = "logFC")
                      ),
                      column(width=2, selectInput("samedir",
                                                  label = "Need logFC in same direction",
                                                  choices = c("Yes"=TRUE, "NO"=FALSE),
                                                  selected = TRUE)
                      ),
                      column(width=2, numericInput("min_rowMean",
                                                  label = "Expression cutoff:",
                                                  min=0, max=10, value=0.5, step=0.05)
                      ),
                      column(width=2, numericInput("qval_cutoff",
                                                  label = "q_val cutoff:",
                                                  min=0.05, max=1, value=0.4, step=0.05)
                      ),
                      column(width=2, selectInput("go_logscale",
                                                  label = "Use Log scale",
                                                  choices = c("Yes"=TRUE, "No"=FALSE),
                                                  selected = TRUE)
                      )
                    ),
                    fluidRow(
                      column(width=2, selectInput("go_colorscale",
                                                  label = "Color Scheme:",
                                                  # choices = c("Greys", "YlGnBu", "Greens", "YlOrRd", "Bluered", "RdBu", "Reds", "Blues", "Picnic", "Rainbow", "Portland", "Jet", "Hot", "Blackbody", "Earth", "Electric", "Viridis"),
                                                  choices = c("YlOrRd", "YlOrBr", "YlGnBu", "PuBuGn", "Spectral", "RdYlGn", "RdYlBu", "RdGy", "RdBu", "PRGn", "BuPu", "BuGn", "Blues", "Greys", "PiYG", "BrBG"),
                                                  selected = "Spectral")
                      ),
                      column(width=2, numericInput("go_xfontsize",
                                                   label = "x-axis font size:",
                                                   min=2, max=30, value=10, step=1)
                      ),
                      column(width=2, numericInput("go_yfontsize",
                                                   label = "y-axis font size:",
                                                   min=2, max=30, value=8, step=1)
                      ),
                      column(width=2, numericInput("go_leftmargin",
                                                   label = "left margin:",
                                                   min=50, max=500, value=200, step=10)
                      ),
                      column(width=2, numericInput("go_barwidth",
                                                   label = "colorbar width:",
                                                   min=10, max=40, value=14, step=2)
                      ),
                      column(width=2, br(),
                             downloadButton("downloadPathData", "Download Result", class="dwnld")
                      )
                    )
                  )
              )
            ),
            fluidRow(
              column(width=4,
                     box(width=NULL, title="GO term", height = "500px",
                         DT::dataTableOutput('go_summary'),
                         tags$hr(),
                         p(
                           class = "text-muted",
                           paste("Note: Make sure you select the correct species.")
                         )
                     )
              ),
              column(width=4,
                     box(width=NULL, title="Summary Plot - Greater", solidHeader=TRUE, status="info", height="500px",
                         plotlyOutput('goplot_hi')
                     )
              ),
              column(width=4,
                     box(width=NULL, title="Summary Plot - Less", solidHeader=TRUE, status="info", height="500px",
                         plotlyOutput('goplot_lo')
                     )
              ),
              conditionalPanel(
                condition = "input.EnrichType == 'KEGG'",
                column(width=12,
                       box(width=NULL, height = "850px",
                           imageOutput("pathviewImg")
                       )
                )
              )
            )
    )
  )
)

dashboardPage(
  header,
  sidebar,
  body
)






#       conditionalPanel(
#         condition = "input.Modules == 'Pathway'",
#         selectizeInput("pathway_group",
#                        label = "Display results for which NMF groups?",
#                        choices=NULL, multiple=FALSE,
#                        options = list(placeholder='Choose a NMF group')
#         ),
#         tags$hr(),
#         selectInput("gmt_type",
#                     label = "Which gene sets file to use?",
#                     choices = c("Preloaded gene sets" = "preload",
#                                 "Upload your own gene set (.gmt) file" = "upload"),
#                     selected = "preload"
#         ),
#         conditionalPanel(
#           condition = "input.gmt_type == 'preload'",
#           selectizeInput("predefined_list",
#                          label = "Select from predefined gene sets",
#                          choices = list(
#                            GSEA = c(`Hallmark`='h.all',
#                                     `c1 - Positional`='c1.all',
#                                     `c2 - Curated`='c2.all',
#                                     `c3 - Motif` = "c3.all",
#                                     `c4 - Computational module` = "c4.all",
#                                     `c5 - Gene Ontology (GO)` = "c5.all",
#                                     `c6 - Oncogenic` = "c6.all",
#                                     `c7 - Immunological` = "c7.all"
#                            )
#                          ),
#                          selected = 'c5.all', multiple=F
#           )
#         ),
#         conditionalPanel(
#           condition = "input.gmt_type == 'upload'",
#           fileInput('gmtfile1', 'Uploda your own gene set file',
#                     accept=c('text/csv',
#                              'text/comma-separated-values,text/plain',
#                              '.csv'))
#         ),
#         downloadButton('downloadPathData', 'Download Enrichment Analysis Result'),
#
#         checkboxInput('moreopt2', 'More Options', FALSE),
#         conditionalPanel(
#           condition = "input.moreopt2 == true",
#           selectInput("piano_algorithm",
#                       label = "What kind of algorithm (for piano)?",
#                       choices = c("fisher", "stouffer", "reporter", "tailStrength", "wilcoxon", "mean", "median", "sum", "maxmean", "gsea", "page"),
#                       selected = "mean"),
#           sliderInput("nPerm",
#                       label = "Number of permutation:",
#                       min=100, max=1000, value=200, step=50),
#           selectInput("pathway_mode",
#                       label = "Enrichment Analysis or Raw LogFC?",
#                       choices = c("Enrichment Analysis" = "enrich",
#                                   "Expression log Fold Change across groups" = "logFC"),
#                       selected = "enrich"
#           )
#         )
#       ) # end if Modules = "Pathway"
#     ),
#     mainPanel(
