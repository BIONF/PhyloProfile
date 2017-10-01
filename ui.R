######## pacman NOT YET WORK WITH shinyapp.io #########
# if(!("pacman" %in% installed.packages())) install.packages("pacman")
# library(pacman)
# p_load(shiny,shinyBS,shinyjs,DT,colourpicker,install=T)
#######################################################
# packages <- c("shiny","shinyBS","shinyjs","DT","colourpicker")
# sapply(packages, require, character.only = TRUE)
#######################################################

if (!require("shiny")) {install.packages("shiny")}
if (!require("shinyBS")) {install.packages("shinyBS")}
if (!require("DT")) {install.packages("DT")}
if (!require("colourpicker")) {install.packages("colourpicker")}
if (!require("shinyjs")) {install.packages("shinyjs")}

### showing spinner while waiting for the profile plot
mycss <- "
#plot-container {
position: relative;
}
#loading-spinner {
position: absolute;
left: 10%;
top: 10%;
z-index: -1;
margin-top: -33px;  /* half of the spinner's height */
margin-left: -33px; /* half of the spinner's width */
}
#plot.recalculating {
z-index: -2;
}
"

############## MAIN UI #####################

shinyUI(fluidPage(
  tags$style(type="text/css", "body {padding-top: 80px;}"),
  tags$head(tags$style(HTML(mycss))),

  # Application title
  titlePanel(""),

  useShinyjs(),
  ################### TOP wellpanel for plot configuration ##########################
  conditionalPanel(
    condition = "input.tabs=='Main profile'",
    wellPanel(
      fluidRow(
        column(2,
               radioButtons(
                 inputId="xAxis",
                 label="Choose type of x-axis:",
                 choices=list("taxa","genes"),
                 selected="taxa",
                 inline=T
               ),
               hr(),
               checkboxInput("autoUpdate",strong(em("Auto update plot")), value = TRUE, width = NULL)
        ),
        column(1,
               numericInput("width","Width (px)",min=600,max=3200,step=50,value=600,width=100),
               actionButton("mainPlotConfig","Appearance")
        ),
        column(1,
               numericInput("height","Height (px)",min=600,max=1600,step=50,value=600,width=100)
        ),
        column(2,
               uiOutput("var1_cutoff.ui")
        ),column(2,
               uiOutput("var2_cutoff.ui")
        ),
        column(2,
               sliderInput("percent",
                           "% of present taxa:", min = 0, max = 1, step = 0.025, value = c(0.0,1.0), width = 200)
        ),
        column(2,
               shinyBS::bsButton("resetMain","Reset cutoffs",style="danger"),
               hr(),
               downloadButton('plotDownload','Download profile'),
               tags$head(
                 tags$style(HTML('#plotDownload{background-color:#A9E2F3}'))
               )
        )
      )
    )
  ),

  conditionalPanel(
    condition = "input.tabs=='Customized profile'",
    wellPanel(
      fluidRow(
        column(2,
               radioButtons(
                 inputId="xAxis_selected",
                 label="Choose type of x-axis:",
                 choices=list("taxa","genes"),
                 selected="taxa",
                 inline=T
               ),
               hr(),
               checkboxInput("autoUpdateSelected",strong(em("Auto update plot")), value = TRUE, width = NULL)
        ),
        column(1,
               numericInput("selectedWidth","Width (px)",min=100,max=1000,step=50,value=600,width=100),
               actionButton("selectedPlotConfig","Appearance")
        ),
        column(1,
               numericInput("selectedHeight","Height (px)",min=100,max=1600,step=50,value=600,width=100)
        ),
        column(2,
               uiOutput("var1Filter.ui")
        ),
        column(2,
               uiOutput("var2Filter.ui")
        ),
        column(2,
               uiOutput("percentFilter.ui")
        ),
        column(2,
               shinyBS::bsButton("resetSelected","Reset cutoffs",style="danger"),
               hr(),
               downloadButton('selectedDownload', 'Download profile'),
               tags$head(
                 tags$style(HTML('#selectedDownload{background-color:#A9E2F3}'))
               )
        )
      )
    )
  ),

  ################### main narvarpage tabs ##########################
  navbarPage(
    em(strong("PhyloProfile Tool")),
    id ="tabs",
    collapsible = TRUE,
    inverse = TRUE,
    fluid = TRUE,
    position = "fixed-top",

    ########## INPUT TAB ###########
    tabPanel(
      "Input & settings",
      column(4,
             strong(h4("Main input:")),
             conditionalPanel(
               condition = "input.do",
               em(strong("RELOAD THIS TOOL TO UPLOAD A NEW INPUT FILE!!!",style = "color:red"))
             ),

             checkboxInput("demo",em(strong("Use demo files"),style = "color:darkblue")),
             uiOutput("mainInputFile.ui"),
             uiOutput("taxaInfoCheck.ui"),
             fluidRow(
               column(5,
                      uiOutput("var1_id.ui")
               ),
               column(5,
                      selectInput("var1_aggregateBy", label = h5("Aggregate by:"),
                                  choices = list("Max"="max", "Min"="min","Mean"="mean","Median"="median"),
                                  selected = "max",
                                  width = 130)
               )
             ),
             fluidRow(
               column(5,
                      uiOutput("var2_id.ui")
               ),
               column(5,
                      selectInput("var2_aggregateBy", label = h5("Aggregate by:"),
                                  choices = list("Max"="max", "Min"="min","Mean"="mean","Median"="median"),
                                  selected = "max",
                                  width = 130)
               )
             ),

             hr(),
             strong(h4("Additional annotation input:")),
             radioButtons(inputId="annoChoose", label="", choices=list("from file","from folder"), selected="from file", inline=T),
             uiOutput("domainInputFile.ui"),

             hr(),
             em("Click here to download demo files:"),
             em(a("(1) Main inputs,", href="https://github.com/BIONF/phyloprofile-data/tree/data/demo", target="_blank")),
             em(a("(2) Domain annotations (optional),", href="https://github.com/BIONF/phyloprofile-data/tree/data/demo/domain_files", target="_blank")),
             em(a("(3) FASTA sequence files (optional)", href="https://github.com/BIONF/phyloprofile-data/tree/data/demo/microsporidia_fasta", target="_blank"))
      ),
      column(3,
             conditionalPanel(
               condition = 'output.unkTaxaStatus == 1',
               strong(h4("New taxa were found:")),
               dataTableOutput("unkTaxaFull")
             ),

             conditionalPanel(
               condition = 'output.unkTaxaStatus == 0',
               strong(h4("Choose genes of interest:")),
               radioButtons(inputId="geneList_selected", label="", choices=list("all","from file"), selected="all", inline=T),
               conditionalPanel(
                 condition = "input.geneList_selected == 'from file'",
                 fileInput("list","")
               ),
               hr(),

               checkboxInput('ordering',strong('Order sequence IDs'),value = TRUE),
               hr(),

               HTML("<b>Order taxa</b>"),
               radioButtons(inputId="order_taxa", label="", choices=list("automatically","by user defined tree"), selected="automatically", inline=T),
               conditionalPanel(
                 condition = "input.order_taxa == 'by user defined tree'",
                 fileInput("inputTree","")
               ),
               hr(),

               shinyBS::bsButton("getConfig","FASTA config"),
               h5(""),
               actionButton("setColor","COLORS config",style='padding:4px; font-size:100%'),
               hr()
             )
      ),
      column(4,
             conditionalPanel(
               condition = 'output.unkTaxaStatus',
               strong(h4("PLEASE CHECK:")),

               em("Do you have any taxon, which doesn't exist in the NCBI taxonomy database?"),
               radioButtons("newTaxaAsk","", c("Yes" = "Yes", "No" = "No"), inline=T, selected = "No"),
               conditionalPanel(condition="input.newTaxaAsk == 'Yes'",
                                shinyBS::bsButton("addTaxa","Add info for new taxa",disabled=FALSE,style="warning")
               ),
               conditionalPanel(condition="input.newTaxaAsk == 'No'",
                                shinyBS::bsButton("BUTparse","Get taxonomy info from NCBI *",disabled=FALSE,style="warning"),
                                helpText(em("(*) Taxonomy information for a given taxa list contains all taxonomy ranks and their correspoding NCBI IDs"))
               ),
               # conditionalPanel(
               #   condition = 'output.invalidIDstatus == 1',
               #   strong(h4("Invalid taxa were found:")),
               #   
               # ),
               
               hr(),
               uiOutput("endParsingMsg"),
               tableOutput("invalidID.output")
               # strong(h4("PLEASE RELOAD THIS TOOL AFTER ADDING NEW TAXA!!!"),style = "color:red")
             ),

             conditionalPanel(
               condition = 'output.unkTaxaStatus == 0',
               strong(h4("Seed (super)taxon:")),
               br(),

               strong(h5("Select taxonomy rank:")),
               uiOutput("rankSelect"),
               br(),
               strong(h5("Choose (super)taxon of interest:")),
               uiOutput("select"),
               br(),
               shinyBS::bsButton("do", "PLOT",type="action",style="danger",size = "large",disabled = TRUE),
               h5("")
             )
      )
    ),

    ########## MAIN PROFILE TAB ###########
    tabPanel(
      "Main profile",
      sidebarLayout(
        sidebarPanel(
          textOutput("testOutput"),    ### use for testing output ###
          column(4,offset = 0,numericInput("number","# of genes:",min=1,max=1600,step=10,value=30,width=150),style='padding:0px;'),
          column(4,numericInput("stIndex","1st index:",min=1,max=1600,value=1,width=200)),

          column(4,uiOutput("highlightGeneUI")),
          #bsPopover("highlightGeneUI","","OR double click on heatmap","right"),

          uiOutput("highlightTaxonUI"),
          bsPopover("highlightTaxonUI","","OR double click on heatmap","right"),

          conditionalPanel(
            condition = "input.autoUpdate == false",
            shinyBS::bsButton("updateBtn","Update plot",style="warning")
          )
        ),

        mainPanel(
          uiOutput("plot.ui"),

          conditionalPanel(
            condition = "input.mainXAxisGuide == true | input.mainYAxisGuide == true",
            absolutePanel(
              id="absAxis",
              bottom = 0, left = 0,
              heigh = NULL, width = NULL,
              fixed = TRUE,
              draggable = TRUE,
              style = "opacity: 0.80",

              uiOutput("mainAxisRender")
            )
          )
        )
      )
    ),

    ########## CUSTOMIZED PROFILE TAB ###########
    tabPanel(
      "Customized profile",
      sidebarLayout(
        sidebarPanel(
          width=4,
          column(12,
                 style='padding:0px;',
                 strong("Select sequence(s) of interest:")
          ),

          column(12,
                 fluidRow(
                   column(8,
                          style='padding:0px;',
                          uiOutput("geneIn")
                   ),
                   column(4,
                          fileInput("customFile","",width='100%')
                   )
                 )
          ),
          column(12,
                 style='padding:0px;',
                 strong("Select (super)taxon/(super)taxa of interest:")
          ),
          column(12,
                 fluidRow(
                   column(8,
                          style='padding:0px;',
                          uiOutput("taxaIn")
                   ),
                   column(4,
                          h3(""),
                          shinyBS::bsButton("cusTaxa","Browse...")
                   )
                 )
          ),

          h5(""),
          uiOutput("plotCustomBtn")
        ),
        mainPanel(
          uiOutput("selectedPlot.ui")
        )
      )
    ),

    ########## FUNCTION TAB ###########
    navbarMenu(
      "Function",
      tabPanel(
        "Profiles clustering",
        h4(strong("Profiles clustering")),

        wellPanel(
          fluidRow(
            column(3,
                   selectInput("distMethod", label = h5("Distance measure method:"),
                               choices = list("euclidean" = "euclidean", "maximum" = "maximum", "manhattan" = "manhattan",
                                              "canberra" = "canberra", "binary" = "binary"),
                               selected = "euclidean")
            ),
            column(3,
                   selectInput("clusterMethod", label = h5("Cluster method:"),
                               choices = list("single" = "single", "complete" = "complete",
                                              "average (UPGMA)" = "average", "mcquitty (WPGMA)" = "mcquitty", "median (WPGMC)" = "median","centroid (UPGMC)"="centroid"),
                               selected = "complete")
            ),
            column(1,
                   numericInput("clusterPlot.width",h5("Width (px)"),min=200,max=3200,step=50,value=600,width=100)
            ),
            column(2,
                   numericInput("clusterPlot.height",h5("Height (px)"),min=200,max=3200,step=50,value=400,width=100)
            ),
            column(3,
                   checkboxInput("applyCluster",em(strong("Apply clustering to profile plot", style="color:red")),value = FALSE),
                   uiOutput("applyClusterCheck.ui"),

                   tags$head(
                     tags$style(HTML('#downloadCluster{background-color:#A9E2F3}'))
                   ),
                   downloadButton('downloadCluster', 'Download plot')
            )
          )
        ),

        column(8,
               uiOutput("cluster.ui")
        ),
        column(4,
               downloadButton('downloadClusterGenes', 'Download gene list'),
               checkboxInput('addClusterCustomProfile',strong(em("Add to Customized profile")), value = FALSE, width = NULL),
               uiOutput("addClusterCustomProfileCheck.ui"),
               tableOutput("brushedCluster.table")
        )
      ),

      tabPanel(
        "Distribution analysis",
        h4(strong("Distribution analysis")),

        wellPanel(
          fluidRow(
            column(2,
                   uiOutput("selected.distribution")
            ),
            column(2,
                   uiOutput("var1_dist.ui")
            ),
            column(2,
                   uiOutput("var2_dist.ui")
            ),
            column(2,
                   uiOutput("percent_dist.ui")
            ),
            column(1,
                   numericInput("dist_textSize","Label size",min=2,max=99,step=1,value=12,width=100)
            ),
            column(2,
                   strong("Download"),
                   tags$head(
                     tags$style(HTML('#plotDownload_dist{background-color:#A9E2F3}'))
                   ),
                   downloadButton('plotDownload_dist','Download plot')
            )
          )
        ),

        uiOutput("dist_plot.ui")
      ),

      tabPanel(
        "Gene age estimation",
        h4(strong("Gene age estimation")),

        wellPanel(
          fluidRow(
            column(2,
                   uiOutput("var1_age.ui")
            ),
            column(2,
                   uiOutput("var2_age.ui")
            ),
            column(2,
                   uiOutput("percent_age.ui")
            ),
            column(2,
                   strong("Download"),
                   tags$head(
                     tags$style(HTML('#geneAgePlotDownload{background-color:#A9E2F3}'))
                   ),
                   downloadButton("geneAgePlotDownload","Download plot")
            )
          )
        ),

        uiOutput("geneAge.ui"),
        conditionalPanel(
          condition = "input.do",
          em(h6("01_Species; 02_Family; 03_Class; 04_Phylum;
                05_Kingdom; 06_Superkingdom; 07_Last universal common ancestor;
                Undef_Genes have been filtered out"))
        ),
        hr(),
        column(4,
               downloadButton("geneAgeTableDownload","Download gene list"),
               checkboxInput("addCustomProfile",strong(em("Add to Customized profile")), value = FALSE, width = NULL),
               uiOutput('addCustomProfileCheck.ui')
        ),
        tableOutput("geneAge.table")
      ),

      tabPanel(
        "Core gene identification",
        h4(strong("Core gene identification")),

        wellPanel(
          fluidRow(
            column(2,
                   uiOutput("var1_cons.ui")
            ),
            column(2,
                   uiOutput("var2_cons.ui")
            ),
            column(2,
                   uiOutput("percent_cons.ui")
            ),
            column(6,
                   uiOutput("taxaList_cons.ui"),
                   shinyBS::bsButton("browseTaxaCons","Browse")
            )
          )
        ),
        hr(),
        column(4,
               downloadButton("consGeneTableDownload","Download gene list"),
               checkboxInput("addConsGeneCustomProfile",strong(em("Add to Customized profile")), value = FALSE, width = NULL),
               uiOutput('addConsGeneCustomProfileCheck.ui')
        ),
        dataTableOutput("consGene.table")
      ),

      tabPanel("Search for NCBI taxonomy IDs",
               column(3,
                      fileInput("taxaList",h4("Upload taxa list")),
                      shinyBS::bsButton("idSearch","Search")
               ),
               column(9,
                      h4("Mismatch(es):"),
                      dataTableOutput("notfoundTaxa"),
                      downloadButton("downloadNotFoundTaxa","Download"),

                      hr(),
                      h4("Retrieved taxonomy ID(s):"),
                      dataTableOutput("taxaID"),
                      downloadButton("downloadTaxaID","Download")
               )
        )
      ),

      ########## DATA TAB ###########
      navbarMenu("Download filtered data",
                 tabPanel("Main data",
                          column(4,
                                 checkboxInput("getRepresentativeMain",strong(em("Download representative sequences")), value = FALSE, width = NULL)
                          ),
                          column(4,
                                 conditionalPanel(
                                   condition = "input.getRepresentativeMain == true",
                                   uiOutput("refVarMain.ui")
                                 )
                          ),
                          column(4,
                                 conditionalPanel(
                                   condition = "input.getRepresentativeMain == true",
                                   radioButtons(inputId="refTypeMain", label="Select representative by", choices=list("max","min"), selected="max", inline=T)
                                 )
                          ),
                          column(12,
                                 dataTableOutput("filteredMainData")
                          ),
                          column(3,
                                 downloadButton('downloadData', 'Download filtered data')
                          ),
                          column(3,
                                 downloadButton('downloadFasta', 'Download FASTA sequences')
                          )
                 ),
                 tabPanel("Customized data",
                          column(12,
                                 conditionalPanel(
                                   condition = "input.getRepresentativeMain == true",
                                   uiOutput("representativeInfo.ui")
                                 )
                          ),
                          # hr(),
                          # dataTableOutput("filteredCustomData"),
                          # downloadButton('downloadCustomData', 'Download customized data')
                          
                          column(12,
                                 dataTableOutput("filteredCustomData")
                          ),
                          column(3,
                                 downloadButton('downloadCustomData', 'Download customized data')
                          ),
                          column(3,
                                 downloadButton('downloadCustomFasta', 'Download FASTA sequences')
                          )
                  )
      ),

      ########## OTHERS TAB ###########
      navbarMenu("More",
                 tabPanel("Help",
                          HTML('<iframe src="https://player.vimeo.com/video/225373912" width="640" height="360" frameborder="0" webkitallowfullscreen mozallowfullscreen allowfullscreen></iframe>'),
                          br(),
                          h3(a("Click here for a detail manual", href="https://trvinh.github.io/phyloprofile_slides/", target="_blank"))
                 ),
                 tabPanel("Q&A",
                          uiOutput("help.ui")
                 ),
                 tabPanel(a("Readme", href="https://BIONF.github.io/PhyloProfile/", target="_blank")
                 )
      )
    ),

    ################### LIST OF POP-UP WINDOWS ##########################

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
    bsModal("parseConfirm", "Get info from input", "BUTparse", size = "small",
            HTML("Processing...<br><br>"),
            strong("PLEASE RELOAD THIS TOOL WHEN FINISHED!!!",style = "color:red"),
            # conditionalPanel(
              # condition = 'output.taxonomyParseStatus == 1',
              # strong(h4("Invalid NCBI IDs:")),
              dataTableOutput("invalidIDout")
            # )
    ),

    ####### popup windows for setting plot colors
    bsModal("color", "Set colors for profile", "setColor", size = "small",
            colourpicker::colourInput("lowColor_var1", "Low variable 1", value = "darkorange"),
            colourpicker::colourInput("highColor_var1", "High variable 1", value = "steelblue"),
            actionButton("defaultColorVar1","Default",style='padding:4px; font-size:100%'),
            hr(),
            colourpicker::colourInput("lowColor_var2", "Low variable 2", value = "grey95"),
            colourpicker::colourInput("highColor_var2", "High variable 2", value = "khaki"),
            actionButton("defaultColorVar2","Default",style='padding:4px; font-size:100%'),
            hr(),
            colourpicker::colourInput("paraColor", "Color for inparalogs", value = "#07d000"),
            actionButton("defaultColorPara","Default",style='padding:4px; font-size:100%')
    ),

    ####### popup windows for FASTA configurations
    bsModal("config", "FASTA config", "getConfig", size = "small",
            selectInput("input_type", "Choose location for:",
                        c("Fasta folder","oneSeq.extended.fa")
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
              ,selectInput("id_format","ID format:",choices=list(">speciesID:seqID"=1,">speciesID@seqID"=2,">speciesID|seqID"=3),selected=1)
            )
    ),

    ####### popup windows for setting main plot configurations
    bsModal("mainPlotConfigBs", "Plot properties configuration", "mainPlotConfig", size = "small",
            column(6,
                   numericInput("xSize","X-axis label size (px)",min=8,max=99,step=1,value=8,width=100)
            ),
            column(6,

                   numericInput("ySize","Y-axis label size (px)",min=8,max=99,step=1,value=8,width=100)
            ),

            column(6,
                   numericInput("legendSize","Legend label size (px)",min=8,max=99,step=1,value=8,width=150)
            ),
            column(6,
                   selectInput("mainLegend", label = "Legend position:",
                               choices = list("Right"="right", "Left"="left","Top"="top","Bottom"="bottom", "Hide"="none"),
                               selected = "right",
                               width = 150)
            ),

            column(12,
                   HTML("<strong>Zooming factor (α) for dots on profile</strong>:<br>"),
                   sliderInput("dotZoom","", min = -1, max = 3, step = 0.1, value = 0, width = 250),
                   HTML("<em>size = (1+α)*default_size<br>default_size=[0:5]</em>"),
                   uiOutput("dotSizeInfo"),
                   br()
            ),

            br(),
            hr(),
            shinyBS::bsButton("resetMainConfig","Reset",style="danger"),
            shinyBS::bsButton("applyMainConfig","Done",style="warning")
    ),

    ####### popup windows for setting main plot configurations
    bsModal("selectedPlotConfigBs", "Plot properties configuration", "selectedPlotConfig", size = "small",
            column(6,
                   numericInput("xSizeSelect","X-axis label size(px)",min=8,max=99,step=1,value=8,width=150)
            ),
            column(6,

                   numericInput("ySizeSelect","Y-axis label size (px)",min=8,max=99,step=1,value=8,width=100)
            ),

            column(6,
                   numericInput("legendSizeSelect","Legend label size (px)",min=8,max=99,step=1,value=8,width=150)
            ),
            column(6,
                   selectInput("selectedLegend", label = "Legend position:",
                               choices = list("Right"="right", "Left"="left","Top"="top","Bottom"="bottom", "Hide"="none"),
                               selected = "right",
                               width = 150)
            ),

            column(12,
                   HTML("<strong>Zooming factor (α) for dots on profile</strong>:<br>"),
                   sliderInput("dotZoomSelect","", min = -1, max = 3, step = 0.1, value = 0, width = 250),
                   HTML("<em>size = (1+α)*default_size<br>default_size=[0:5]</em>"),
                   uiOutput("dotSizeInfoSelect"),
                   br()
            ),

            br(),
            hr(),
            shinyBS::bsButton("resetSelectedConfig","Reset",style="danger"),
            shinyBS::bsButton("applySelectedConfig","Done",style="warning")
    ),

    ####### popup windows for select taxa on Customized Profile
    bsModal("cusTaxaBS", "Select taxon/taxa of interest", "cusTaxa", size = "small",
            uiOutput("rankSelectCus"),
            uiOutput("taxaSelectCus"),
            checkboxInput("applyCusTaxa",strong("Apply to customized profile", style="color:red"),value = FALSE)
    ),

    ####### popup windows for select taxa on Consensus gene finding
    bsModal("browseTaxaConsBS", "Select taxon/taxa of interest", "browseTaxaCons", size = "small",
            uiOutput("rankSelectCons"),
            uiOutput("taxaSelectCons"),
            checkboxInput("applyConsTaxa",strong("Apply", style="color:red"),value = FALSE)
    ),

    ####### popup windows for detailed plot
    bsModal("modalBS", "Detailed plot", "go", size = "large",
            uiOutput("detailPlot.ui"),
            numericInput("detailedHeight","plot_height(px)",min=100,max=1600,step=50,value=100,width=100)
            ,verbatimTextOutput("detailClick")
            ,shinyBS::bsButton("doDomainPlot", "Show domain architecture",disabled = TRUE)
            ,uiOutput("checkDomainFiles")
            ,br()
            ,h4("Sequence:")
            ,verbatimTextOutput("fasta")
    ),

    ####### popup windows for domain architecture plot
    bsModal("plotArchi","Domain architecture","doDomainPlot", size = "large",
            fluidRow(
              column(2,
                     numericInput("archiHeight","plot_height(px)",min=100,max=1600,step=50,value=400,width=100)
              ),
              column(2,
                     numericInput("archiWidth","plot_width(px)",min=100,max=1600,step=50,value=800,width=100)
              ),
              column(2,
                     numericInput("titleArchiSize","Title size(px)",min=8,max=99,step=1,value=11,width=150)
              ),
              column(2,
                     numericInput("labelArchiSize","SeqID size(px)",min=8,max=99,step=1,value=11,width=150)
              ),
              column(2,
                     numericInput("labelDescSize","Text size(px)",min=0,max=99,step=1,value=3,width=150)
              )
            ),
            uiOutput("archiPlot.ui"),
            downloadButton("archiDownload","Download plot")
    ),

    ################### POINT INFO BOX ##########################
    conditionalPanel(
      condition = "input.tabs=='Main profile' || input.tabs=='Customized profile'",
      ############# PONIT's INFO BOX
      absolutePanel(
        bottom = 5, left = 30,
        fixed = TRUE,
        draggable = TRUE,
        h5("Point's info:"),
        verbatimTextOutput("pointInfo"),
        shinyBS::bsButton("go", "Detailed plot", style="success", disabled = FALSE),
        style = "opacity: 0.80"
      )
    )
  )
)
