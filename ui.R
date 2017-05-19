if (!require("shiny")) {install.packages("shiny")}
if (!require("shinyBS")) {install.packages("shinyBS")}
if (!require("DT")) {install.packages("DT")}
if (!require("colourpicker")) {install.packages("colourpicker")}
if (!require("shinyjs")) {install.packages("shinyjs")}

### showing waiting spinner while plotting
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
               checkboxInput("autoUpdate",strong(em("Auto update plot")), value = TRUE, width = NULL),
               bsButton("resetMain","Reset",style="danger")
               # ,
               # checkboxInput("mainXAxisGuide","X-axis guide", value = FALSE, width = NULL),
               # checkboxInput("mainYAxisGuide","Y-axis guide", value = FALSE, width = NULL)
        ),
        column(1,
               numericInput("width","Width (px)",min=600,max=3200,step=50,value=600,width=100),
               numericInput("height","Height (px)",min=600,max=1600,step=50,value=600,width=100)
        ),
        column(2,
               numericInput("xSize","X-axis label size(px)",min=8,max=99,step=1,value=8,width=150),
               numericInput("ySize","Y-axis label size(px)",min=8,max=99,step=1,value=8,width=150)
        ),
        column(2,
               numericInput("legendSize","Legend label size(px)",min=8,max=99,step=1,value=8,width=150),
               selectInput("mainLegend", label = "Legend position:",
                           choices = list("Right"="right", "Left"="left","Top"="top","Bottom"="bottom", "Hide"="none"), 
                           selected = "right",
                           width = 150)
        ),
        column(2,
               uiOutput("var1_cutoff"),
               uiOutput("var2_cutoff")
        ),
        column(2,
               sliderInput("percent",
                           "% of present taxa:", min = 0, max = 1, step = 0.025, value = c(0.0,1.0), width = 200),
               tags$head(
                 tags$style(HTML('#plotDownload{background-color:#A9E2F3}'))
               ),
               downloadButton('plotDownload','Download heatmap')
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
               checkboxInput("autoUpdateSelected",strong(em("Auto update plot")), value = TRUE, width = NULL),
               bsButton("resetSelected","Reset",style="danger")
        ),
        column(1,
               numericInput("selectedWidth","Width (px)",min=100,max=1000,step=50,value=600,width=100),
               numericInput("selectedHeight","Height (px)",min=100,max=1600,step=50,value=600,width=100)
        ),
        column(2,
               numericInput("xSizeSelect","X-axis label size(px)",min=8,max=99,step=1,value=8,width=150),
               numericInput("ySizeSelect","Y-axis label size(px)",min=8,max=99,step=1,value=8,width=150)     
        ),
        column(2,
               numericInput("legendSizeSelect","Legend label size(px)",min=8,max=99,step=1,value=8,width=150),
               selectInput(
                 "selectedLegend", label = "Legend position:",
                 choices = list("Right"="right", "Left"="left","Top"="top","Bottom"="bottom", "Hide"="none"), 
                 selected = "right",
                 width = 150)
        ),
        column(2,
               uiOutput("var1Filter.ui"),
               uiOutput("var2Filter.ui")
        ),
        column(2,
               uiOutput("percentFilter.ui"),
               tags$head(
                 tags$style(HTML('#selectedDownload{background-color:#A9E2F3}'))
               ),
               downloadButton('selectedDownload', 'Download heatmap')
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
             
             fileInput("mainInput",h5("Upload input file:")),
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
                                  width = 150)
               )
             ),
             
             hr(),
             strong(h4("Additional annotation file:")),
             fileInput("fileDomain",""),
             hr(),
             em(a("Click here to download demo data", href="https://github.com/trvinh/phyloprofile/tree/master/data/demo", target="_blank"))
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
               
               strong(h4("Ordering sequence IDs by:")),
               radioButtons("ordering", "", choices = c("alphabetical","hierarchical cluster","none"), selected = "alphabetical",
                            inline = F),
               hr(),
               
               bsButton("getConfig","FASTA config"),
               h5(""),
               actionButton("setColor","COLORS config",style='padding:4px; font-size:100%'),
               hr()
             )
      ),
      column(4,
             conditionalPanel(
               condition = 'output.unkTaxaStatus',
               strong(h4("PLEASE CHECK:")),
               em("Does the taxa in your recently uploaded presence/absence file change? (Note: for 'first time users', please choose 'YES')"),
               radioButtons("parseAsk","", c("Yes" = "Yes", "No" = "No"), inline=T, selected = "Yes"),
               conditionalPanel(condition="input.parseAsk == 'Yes'",
                                bsButton("parse","Get info from input",disabled=TRUE)
               ),
               h5(""),
               
               em("Do you have any taxon, which doesn't exist in the NCBI taxonomy database?"),
               radioButtons("newTaxaAsk","", c("Yes" = "Yes", "No" = "No"), inline=T, selected = "No"),
               conditionalPanel(condition="input.newTaxaAsk == 'Yes'",
                                bsButton("addTaxa","Add new taxa",disabled=TRUE)
               ),
               hr(),
               
               strong(h4("PLEASE RELOAD THIS TOOL AFTER ADDING NEW TAXA!!!"),style = "color:red")
             ),
             
             conditionalPanel(
               condition = 'output.unkTaxaStatus == 0',
               strong(h4("Seed (super)taxon:")),
               uiOutput("rankSelect"),
               uiOutput("select"),
               h5(""),
               bsButton("do", "PLOT",type="action",style="danger",size = "large",disabled = TRUE),
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
            bsButton("updateBtn","Update plot",style="warning")
          )
        ),
        
        mainPanel(
          tabsetPanel(
            tabPanel("Main plot",
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
            ),
            
            tabPanel("Distribution plot",
                     uiOutput("selected.distribution"),
                     conditionalPanel(
                       condition = "input.selected_dist == input.var1_id",
                       downloadButton("var1Download","Download"),
                       uiOutput("var1Dist.ui")
                     ),
                     
                     conditionalPanel(
                       condition = "input.selected_dist == input.var2_id",
                       downloadButton("var2Download","Download"),
                       uiOutput("var2Dist.ui")
                     ),
                     
                     conditionalPanel(
                       condition = "input.selected_dist == '% present taxa'",
                       downloadButton("presSpecDownload","Download"),
                       uiOutput("presSpec.ui")
                     )
            ),
            
            tabPanel("Gene age estimation",
                     downloadButton("geneAgePlotDownload","Download plot"),
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
                            checkboxInput("addCustomProfile",strong(em("Add to Customized profile")), value = FALSE, width = NULL)
                     ),
                     tableOutput("geneAge.table"),
                     hr()   
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
          strong("Select sequence(s) of interest:"),
          column(6,
                 style='padding:0px;',
                 uiOutput("geneIn")
          ),
          column(3,
                 fileInput("customFile","",width='100%')
          ),
          uiOutput("taxaIn"),
          
          h5(""),
          uiOutput("plotCustomBtn")
        ),
        mainPanel(
          uiOutput("selectedPlot.ui")
        )
      )
    ),
    
    ########## DATA TAB ###########
    navbarMenu("Download filtered data",
               tabPanel("Main data",
                        dataTableOutput("filteredMainData"),
                        downloadButton('downloadData', 'Download filtered data')
               ),
               tabPanel("Customized data",
                        dataTableOutput("filteredCustomData"),
                        downloadButton('downloadCustomData', 'Download customized data')
               )
    ),
    
    ########## OTHERS TAB ###########
    navbarMenu("More",
               tabPanel("Description",
                        HTML('<iframe width="560" height="315" src="https://www.youtube.com/embed/Udt316KoM6Y" frameborder="0" allowfullscreen></iframe>')
               ),
               tabPanel("Q&A",
                        uiOutput("help.ui")
               ),
               tabPanel("Search for NCBI taxonomy IDs",
                        column(3,
                               fileInput("taxaList",h4("Upload taxa list")),
                               bsButton("idSearch","Search")
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
               ),
               tabPanel(a("About", href="https://trvinh.github.io/phyloprofile/", target="_blank")
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
  bsModal("parseConfirm", "Get info from input", "parse", size = "medium",
          HTML("Parsing taxonomy information from input file?"),
          actionButton("BUTyes", "Yes"),
          actionButton("BUTno", "No"),
          helpText(em("***Note: Please run this step whenever you have a new taxa set. For instance, if you have a new matrix file but the taxa remain the same, then DO NOT re-run this step!***"))
  ),
  
  ####### popup windows for setting plot colors
  bsModal("color", "Set colors for profile", "setColor", size = "small",
          colourpicker::colourInput("lowColor_var1", paste("Low variable 1"), value = "darkorange"),
          colourpicker::colourInput("highColor_var1", "High variable 1", value = "steelblue"),
          actionButton("defaultColorVar1","Default",style='padding:4px; font-size:100%'),
          hr(),
          colourpicker::colourInput("lowColor_var2", "Low variable 2", value = "grey95"),
          colourpicker::colourInput("highColor_var2", "High variable 2", value = "khaki"),
          actionButton("defaultColorTrace","Default",style='padding:4px; font-size:100%')
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
  
  ####### popup windows for detailed plot
  bsModal("modalBS", "Detailed plot", "go", size = "large",
          uiOutput("detailPlot.ui"),
          numericInput("detailedHeight","plot_height(px)",min=100,max=1600,step=50,value=100,width=100)
          ,verbatimTextOutput("detailClick")
          ,actionButton("do3", "Show domain architecture")
          ,br()
          ,h4("Sequence:")
          ,verbatimTextOutput("fasta")
  ),
  
  ####### popup windows for domain architecture plot
  bsModal("plotArchi","Domain architecture","do3", size = "large",
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
      h5("Point's info:"),
      verbatimTextOutput("pointInfo"),
      bsButton("go", "Detailed plot", style="success", disabled = FALSE),
      style = "opacity: 0.80"
    )
  )
  
  ############# Axis guide BOX
  # conditionalPanel(
  #   condition = "input.tabs=='Main profile'",
  #   
  #   absolutePanel(
  #     bottom = 0, left = 130,
  #     fixed = TRUE,
  #     draggable = TRUE,
  #     column(3,
  #            checkboxInput("mainXAxisGuide","X-axis guide", value = FALSE, width = NULL)
  #     ),
  #     column(3,style='padding:0px;',
  #            checkboxInput("mainYAxisGuide","Y-axis guide", value = FALSE, width = NULL)
  #     ),
  #     style = "opacity: 0.80"
  #   )
  # )
  
  
  )
)


