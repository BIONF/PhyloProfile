if (!require("shiny")) {install.packages("shiny")}
if (!require("shinyBS")) {install.packages("shinyBS")}
if (!require("DT")) {install.packages("DT")}
if (!require("colourpicker")) {install.packages("colourpicker")}
if (!require("shinyjs")) {install.packages("shinyjs")}

shinyUI(fluidPage(
  # Application title
  titlePanel("Phylogenetic profile app"),
  
  ################### wellpanel for INPUTS ##########################
  wellPanel(
    fluidRow(
      column(3, offset = 0,
             fileInput("file1","Presence/absence file: "),
             bsButton("addTaxa","Add new taxa",disabled=TRUE),
             bsButton("parse","Get info from input",disabled=TRUE),
             h5(""),
             bsButton("AddFile","Upload additional file(s)",disabled = TRUE)
      ),
      column(3,
             shinyjs::useShinyjs(),
             uiOutput("rankSelect"),
             uiOutput("select"),
             bsButton("getConfig","FASTA config",style="info")
      ),
      column(1,
             numericInput("number","# rows ",min=1,max=1600,step=10,value=30,width=100),
             numericInput("stIndex","start at:",min=1,max=1600,value=1,width=100),
             bsButton("do", "PLOT",type="action",style="danger",size = "large",disabled = TRUE)
      ),
      column(1,
             numericInput("width","Width(px)",min=600,max=3200,step=50,value=600,width=100),
             numericInput("height","Height(px)",min=600,max=1600,step=50,value=600,width=100),
             actionButton("setSize","Set label size",style='padding:4px; font-size:100%')
      ),
      column(1,
             # radioButtons(
             #   inputId="xAxis",
             #   label="x-Axis:",
             #   choices=list("taxa","genes"),
             #   selected="taxa"),
             selectInput("xAxis", label = "x-Axis:",
                         choices = list("Taxa"="taxa", "Genes"="genes"), 
                         selected = "taxa",
                         width = 80),
             selectInput("legendPos", label = "Legend:",
                         choices = list("Right"="right", "Left"="left","Top"="top","Bottom"="bottom", "Hide"="none"), 
                         selected = "right",
                         width = 80),
             actionButton("setColor","Set colors",style='padding:4px; font-size:100%')
      ),
      column(2,
             sliderInput("fas",
                         "FAS cutoff: ",
                         min = 0,
                         max = 1,
                         step = 0.05,
                         value = 0.0,
                         width = 200),
             sliderInput("percent",
                         "% of present species:",
                         min = 0,
                         max = 1,
                         step = 0.05,
                         value = 0.0,
                         width = 200)
        )
      )
  ),
  
  ####### popup to confirm parsing data from input file
  bsModal("addTaxaWindows", "Add new taxa", "addTaxa", size = "medium",
          helpText(em("Use this form to add taxon that does not exist in NCBI taxonomy database")),
          textInput("newID","ID (must be a number and greater than 1835343, e.g. 2000001)",2000001,width=500),
          textInput("newName","Name (e.g. Saccharomyces cerevisiae strain ABC)","",width=500),
          textInput("newRank","Rank (e.g. \"norank\" (for strain), species, order, etc.)","norank",width=500),
          textInput("newParent","Parent ID (NCBI taxonomy ID of the next higher rank, e.g. 4932 (S.cerevisiae species))",4932,width=500),
          actionButton("newAdd","Add"),
          actionButton("newDone","Done")
  ),
  
  ####### popup to confirm parsing data from input file
  bsModal("parseConfirm", "Get info from input", "parse", size = "medium",
          HTML("Parsing taxonomy information from input file?"),
          actionButton("BUTyes", "Yes"),
          actionButton("BUTno", "No"),
          helpText(em("***Note: Please run this step whenever you have a new taxa set. For instance, if you have a new matrix file but the taxa remain the same, then DO NOT re-run this step!***"))
  ),

  ####### popup windows for uploading additional files (FAS & Traceability)
  bsModal("modalnew", "Upload additional file(s)", "AddFile", size = "medium",
          fileInput("file3","Feature architectures:"),
          fileInput("file2","Traceability matrix:")
  ),
  
  ####### popup windows for upload list of genes of interest
  bsModal("geneListBs", "Gene list", "geneList", size = "small",
          radioButtons(
            inputId="geneList_selected",
            label="Select list of genes of interest:",
            choices=list(
              "all",
              "from file"
            ),
            selected="all"),
          conditionalPanel(
            condition = "input.geneList_selected == 'from file'",
            fileInput("list","List of genes of interest:")
          )
  ),
  
  ####### popup windows for setting plot colors
  bsModal("color", "Set colors for profile", "setColor", size = "small",
          shinyjs::useShinyjs(),
          colourpicker::colourInput("lowColor_fas", "Low FAS", value = "darkorange"),
          colourpicker::colourInput("highColor_fas", "High FAS", value = "steelblue"),
          actionButton("defaultColorFas","Default",style='padding:4px; font-size:100%'),
          hr(),
          colourpicker::colourInput("lowColor_trace", "Low traceability", value = "grey95"),
          colourpicker::colourInput("highColor_trace", "High traceability", value = "khaki"),
          actionButton("defaultColorTrace","Default",style='padding:4px; font-size:100%')
  ),
  
  ####### popup windows for setting axis label size
  bsModal("axisSize", "Set size for", "setSize", size = "small",
          numericInput("xSize","X-axis label (px)",min=8,max=99,step=1,value=8,width=200),
          numericInput("ySize","Y-axis label (px)",min=8,max=99,step=1,value=8,width=200),
          actionButton("defaultSize","Default",style='padding:4px; font-size:100%')
  ),
  
  ####### popup windows for FASTA configurations
  bsModal("config", "FASTA config", "getConfig", size = "small",
          selectInput("input_type", "Choose location for:",
                      c("oneSeq.extended.fa", "Fasta folder")
          ),
          hr(),
          conditionalPanel(
            condition = "input.input_type == 'oneSeq.extended.fa'",
#            textInput("oneseq.file","Path:",""),
            fileInput("oneSeqFasta",""),
            uiOutput("oneSeq.existCheck")

          ),
          conditionalPanel(
            condition = "input.input_type == 'Fasta folder'",
            textInput("path","Main path:","")
            ,selectInput("dir_format","Directory format:",choices=list("path/speciesID.fa*"=1,"path/speciesID/speciesID.fa*"=2),selected="Path/speciesID.fasta")
            ,selectInput("file_ext","File extension:",choices=list("fa"="fa","fasta"="fasta","fas"="fas","txt"="txt"),selected="fa")
            ,selectInput("id_format","ID format:",choices=list(">speciesID:seqID"=1,">seqID"=2),selected=2)
          )
  ),
  
  ################### SIDEBAR PANEL AND MAIN PANEL ################### 
  sidebarLayout(
    fluid = FALSE,
    
    sidebarPanel(
      textOutput("testOutput"),    ### use for testing output ###
      uiOutput("highlight"),
      uiOutput("geneIn"),
      uiOutput("taxaIn"),
      h5(""),
      bsButton("geneList","Upload sequence list",style="warning"),
      bsButton("do2", "Plot selected sequence(s)",disabled=TRUE)
    ),
    
    ################ Main page
    mainPanel(
      tabsetPanel(
#        tabPanel ("Distribution",plotOutput("plot1")),
        tabPanel ("Presence/absence profile",
                  uiOutput("plot.ui"),
                  downloadButton('plotDownload','Download plot'),
                  bsModal("modalBS", "Detailed plot", "go", size = "large",
                          uiOutput("detailPlot.ui"),
                          numericInput("detailedHeight","plot_height(px)",min=100,max=1600,step=50,value=100,width=100)
                          ,verbatimTextOutput("detailClick")
                          ,actionButton("do3", "Show domain architecture")
                          ,br()
                          ,h4("Sequence:")
                          ,verbatimTextOutput("fasta")
                  ),
                  bsModal("helpBS", "Help", "help", size = "large",
                          uiOutput("help.ui")
                  ),
                  bsModal("plotSeq","Plot selected sequence","do2", size = "large",
                          uiOutput("selectedPlot.ui"),
                          fluidRow(
                            column(2,
                                 br(),
                                 numericInput("selectedHeight","Plot_height(px)",min=100,max=1600,step=50,value=400,width=100),
                                 numericInput("selectedWidth","Plot_Width(px)",min=100,max=1000,step=50,value=800,width=100)
                            ),
                            column(2,
                                   br(),
                                   radioButtons(
                                     inputId="xAxis_selected",
                                     label="x-Axis:",
                                     choices=list(
                                       "taxa",
                                       "genes"
                                     ),
                                     selected="taxa"),
                                   radioButtons(
                                     inputId="legend",
                                     label="Legend pos:",
                                     choices=list(
                                       "top",
                                       "right"
                                     ),
                                     selected="top")
                            ),
                            column(5,
                                   br(),
                                   HTML("<strong>Point's info:</strong>"),
                                   verbatimTextOutput("selectedClick")
                            ),
                            column(2,
                                   br(),
                                   br(),
                                   actionButton("selectedDownload","Download plot",disabled=TRUE)
                            )
                          ),
                          verbatimTextOutput("fasta_selected")
                  ),
                  bsModal("plotArchi","Domain architecture","do3", size = "large",
                          uiOutput("archiPlot.ui"),
                          numericInput("archiHeight","plot_height(px)",min=100,max=1600,step=50,value=400,width=100),
                          numericInput("archiWidth","plot_width(px)",min=100,max=1600,step=50,value=800,width=100)
                  )
        ),
        tabPanel ("Data",dataTableOutput("dis"),
                  downloadButton('downloadData', 'Download filtered data'))
        )
    )
  ),
  
  ############# HELP button
  absolutePanel(
    bottom = 5, right = 30,
    fixed = TRUE,
    actionButton("help", "HELP",style='padding:4px; font-size:80%'),
    style = "opacity: 0.80"
  ),

  ############# PONIT's INFO BOX
  absolutePanel(
    bottom = 5, left = 30,
    fixed = TRUE,
    h5("Point's info:"),
    verbatimTextOutput("pointInfo"),
    bsButton("go", "Detailed plot", style="success", disabled = FALSE),
    style = "opacity: 0.80"
  )
))
