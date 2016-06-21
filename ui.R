library(shiny)
library(shinyBS)
library(DT)

shinyUI(fluidPage(
  # Application title
  titlePanel("Phylogenetic profile App"),
  
  # Sidebar for input file
  wellPanel(
  fluidRow(
    column(2, offset = 0,
           fileInput("file1","Upload input matrix file: "),
           actionButton("parse","Get info from input",style='padding:4px; font-size:85%')
    ),
    column(2,
           uiOutput("rankSelect")
    ),
    column(2,
           uiOutput("select"),
           actionButton("do", "Plot")
    ),
    column(2,
           numericInput("number","# rows for profile",min=1,max=1600,step=10,value=30,width=150),
           textOutput("totalRows"),
           numericInput("stIndex","start at:",min=1,max=1600,value=1,width=100)
    ),
    column(1,
           numericInput("width","Width(px)",min=600,max=3200,step=50,value=600,width=100),
           numericInput("height","Height(px)",min=600,max=1600,step=50,value=600,width=100)
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
  
  sidebarLayout(
    fluid = FALSE,
    
    sidebarPanel(
      textOutput("testOutput"),    ### use for testing output ###
      uiOutput("highlight"),
      uiOutput("geneIn"),
      actionButton("do2", "Plot this selected sequence")
    ),
    
    # Main page
    mainPanel(
      tabsetPanel(
#        tabPanel ("Distribution",plotOutput("plot1")),
        tabPanel ("Present/absent profile",uiOutput("plot.ui"),
                  downloadButton('plotDownload','Download plot'),
                  bsModal("modalBS", "Detailed plot", "go", size = "large",
                          uiOutput("detailPlot.ui"),
                          numericInput("detailedHeight","plot_height(px)",min=100,max=1600,step=50,value=100,width=100)
                          ),
                  bsModal("helpBS", "Help", "help", size = "large",
                          uiOutput("help.ui")
                  ),
                  bsModal("plotSeq","Plot selected sequence","do2", size = "large",
                          uiOutput("selectedPlot.ui"),
                          numericInput("selectedHeight","Plot_height(px)",min=100,max=1600,step=50,value=400,width=100),
                          numericInput("selectedWidth","Plot_Width(px)",min=100,max=1600,step=50,value=800,width=100)
                  )
                ),
        tabPanel ("Data",dataTableOutput("dis"),
                  downloadButton('downloadData', 'Download filtered data'))
      )
    )
  ),
  
  ### help button
  absolutePanel(
    bottom = 5, right = 30,
    fixed = TRUE,
    actionButton("help", "HELP",style='padding:4px; font-size:80%'),
    style = "opacity: 0.80"
  ),

  ### show click info
  absolutePanel(
    bottom = 5, left = 30,
    fixed = TRUE,
#    wellPanel(
      h5("Point's info:"),
      verbatimTextOutput("pointInfo"),
      actionButton("go", "Detailed plot"),
#    ),
    style = "opacity: 0.80"
  )
))
