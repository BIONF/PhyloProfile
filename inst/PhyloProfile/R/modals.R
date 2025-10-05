
# * popup for plotting detailed plot -------------------------------------------
detailedPlotModal <- function() {
    modalDialog(
        title = "Detailed plot",
        size = "l",
        easyClose = TRUE,
        footer = modalButton("Close"),
        
        fluidRow(
            column(
                2, createPlotSize("detailedHeight", "Height (px)", 100)
            ),
            column(
                3, createTextSize("detailedText", "Text size (px)", 12, 150)
            ),
            column(
                7,
                checkboxInput(
                    "detailedRemoveNA",
                    strong("Hide taxa that have no ortholog (NAs)",
                           style = "color:red"),
                    value = FALSE
                ),
                checkboxInput(
                    "detailedFilter",
                    strong("Apply filters",
                           style = "color:red"),
                    value = FALSE
                )
            )
        ),
        hr(),
        createDetailedPlotUI("detailedPlot"),
        actionButton("doDomainPlot", "Show domain architecture"),
        uiOutput("checkDomainFiles"),
        br(),
        h4("Sequence:"),
        verbatimTextOutput("fasta"),
        br(),
        h4("Links:"),
        uiOutput("dbLink.ui")
    )
}

# * popup for plotting domain architecture plot --------------------------------
archiPlotModal <- function(id, title = "Domain architecture") {
    modalDialog(
        title = title,
        size = "l",
        easyClose = TRUE,
        footer = modalButton("Close"),
        createArchitecturePlotUI("archiPlot")
    )
}

# * popup for confirming parsing taxa from input file ------------------
parseConfirmModal <- function() {
    modalDialog(
        title = "Get taxonomy info",
        size = "s",
        easyClose = FALSE,
        footer = NULL,
        HTML(
            '<p>Fetching Missing Taxonomy Information and
      Post-processing.</p><p><em>This window will close
      automatically when everything is done. Please wait...</em></p>
      <p><strong><span style="color: #ff0000;">PLEASE RELOAD THIS
      TOOL WHEN FINISHED!!!</span></strong></p>'
        )
    )
}

# * popup for upload gene names ----------------------------------------
uploadGeneNameModal <- function() {
    modalDialog(
        title = "Upload gene names",
        size = "m",
        em(paste(
            "This option is used to map gene IDs from the main input with",
            "their gene names. Please upload the mapping file for gene IDs - ",
            "gene names in tab-delimited format! Check",
            "https://github.com/BIONF/PhyloProfile/wiki/Input-Data#gene-names",
            "for more info."
        )),
        br(),
        fileInput("geneName", "")
    )   
}

# * popup for upload gene category -------------------------------------
uploadGeneNameModal <- function() {
    modalDialog(
        title = "Upload gene categories",
        size = "s",
        em(paste(
            "This option is used to map gene IDs from the main input with",
            "a category (e.g. metabolic pathway, phenotype, lifestyle, etc.). ",
            "Please upload the mapping file for gene ID - ",
            "category in tab-delimited format! Check",
            "https://github.com/BIONF/PhyloProfile/wiki/Input-Data#gene-categories",
            "for more info."
        )),
        br(),
        fileInput("geneCategory", "")
    )   
}

# * popup for FASTA upload ---------------------------------------------
uploadFastaModal <- function() {
    modalDialog(
        title = "FASTA upload",
        size = "s",
        selectInput(
            "inputType", "Choose location for:",
            c("Concatenated fasta file", "Fasta folder")
        ),
        hr(),
        conditionalPanel(
            condition = "input.inputType == 'Concatenated fasta file'",
            fileInput("concatFasta", ""),
            uiOutput("concatFasta.existCheck")
        ),
        conditionalPanel(
            condition = "input.inputType == 'Fasta folder'",
            textInput("path", "Main FULL path:", ""),
            selectInput(
                "dirFormat", "Directory format:",
                choices = list("path/speciesID.fa*" = 1,
                               "path/speciesID/speciesID.fa*" = 2),
                selected = "Path/speciesID.fasta"
            ),
            selectInput(
                "fileExt", "File extension:",
                choices = list("fa" = "fa",
                               "fasta" = "fasta",
                               "fas" = "fas",
                               "txt" = "txt"),
                selected = "fa"
            ),
            selectInput(
                "idFormat",
                "ID format:",
                choices = list(">speciesID:seqID" = 1,
                               ">speciesID@seqID" = 2,
                               ">speciesID|seqID" = 3,
                               ">seqID" = 4),
                selected = 4
            )
        )
    )
}

# * popup for adding new taxa from input file --------------------------
addTaxaModal <- function() {
    modalDialog(
        title = "Add new taxa",
        size = "m",
        HTML(
            "<p><em>Use this form to add taxon that does not exist in NCBI
                taxonomy database (or alternatively you can manually prepare the
                <span style=\"text-decoration: underline;\">
                <span style=\"color: #ff0000; text-decoration: underline;\">
                PhyloProfile/data/newTaxa.txt file with the following
                description for each field).</em></p>
                <p><span style=\"color: #ff0000;\"><em><strong>
                NOTE: ID and name of new taxon must be
                <span style=\"text-decoration: underline;\">different</span>
                from any existing NCBI taxa.</strong></em></span></p>"
        ),
        textInput(
            "newID",
            "ID (must be a number and greater than 999999005,
                e.g. 999999901)",
            999999901,
            width = 500
        ),
        textInput(
            "newName",
            "Name (e.g. Saccharomyces cerevisiae strain ABC)",
            "",
            width = 500
        ),
        textInput(
            "newRank",
            "Rank (e.g. strain, species, order, etc.)",
            "species",
            width = 500
        ),
        textInput(
            "newParent",
            "Parent ID (NCBI taxonomy ID of the next higher rank,
                e.g. 4932 (S.cerevisiae species))",
            4932,
            width = 500
        ),
        actionButton("newAdd", "Add new taxon"),
        hr(),
        fileInput("newTaxaFile",
                  "Or upload file contains IDs for new taxa"),
        HTML(
            "<p><em>Taxonomy file for new taxa has to be a tab-delimited
                text file and has the following header (please follow the rule
                above):</em></p><p>ncbiID &nbsp;fullName &nbsp;rank
                &nbsp;parentID</p>"
        ),
        uiOutput("wrongNewTaxaUI"),
        hr(),
        actionButton("newDone", "Finish adding", class = "btn btn-warning")
    )
}

# * popup for getting taxa from OMA browser ----------------------------
getOmaModal <- function(){
    modalDialog(
        title = "Get OMA data",
        size = "s",
        selectInput(
            "selectedOmaType",
            label = "Select type of OMA orthologs:",
            choices = list("HOG", "OG"),# "PAIR"),
            selected = "HOG"
        ),
        actionButton("getDataOma", "Get data", class = "btn btn-danger"),
        downloadButton("downloadFilesOma", "Save data"),
        br(),
        em("This windows will close automatically when eveything is done!",
           style = "color:red")
    )
}

# * popup for setting ortholog ID format -------------------------------
seqIdFormatModal <- function() {
    modalDialog(
        title = "Sequence ID format",
        size = "s",
        selectInput(
            "seqIdFormat",
            "ID format:",
            choices = list(
                "BIONF format (seed|taxon|ortho)" = 1,
                "seqID" = 2,
                "something<separator>seqID" = 3,
                "seqID<separator>something" = 4
            ),
            selected = 1
        ),
        selectInput(
            "separator",
            "Separator:",
            choices = list(
                "none" = "none",
                "|" = 1,
                "@" = 2,
                "#" = 3,
                ";" = 4
            ),
            selected = 1
        ),
        em("Please note! Only these ID formats are accepted!")
    )
}

# * popup for setting Gene age plot configurations ---------------------
geneAgeConfigModal <- function() {
    modalDialog(
        title = "Plot appearance configuration",
        size = "s",
        sliderInput(
            "geneAgeWidth",
            "Width zoom (*600px)",
            min = 0.1,
            max = 5,
            step = 0.1,
            value = 1,
            width = "100%"
        ),
        sliderInput(
            "geneAgeHeight", "Height zoom (*150px)",
            min = 0.1,
            max = 5,
            step = 0.1,
            value = 1,
            width = "100%"
        ),
        sliderInput(
            "geneAgeText", "Text size zoom",
            min = 0.1,
            max = 5,
            step = 0.1,
            value = 1,
            width = "100%"
        ),
        br(),
        hr(),
        actionButton("resetGeneAgeProtConfig", "Reset", class = "btn btn-danger")    
    )
}

# * popup for setting Group compariosn plot configurations -------------
gcPlotConfigModal <- function() {
    modalDialog(
        title = "Plot appearance configuration",
        size = "l",
        column(
            6,
            createTextSize("xSizeGC", "X-axis label size (px)", 12, '100%')
        ),
        column(
            6,
            createTextSize("ySizeGC", "Y-axis label size (px)", 12, '100%')
        ),
        column(
            6,
            createTextSize("titleSizeGC", "Title size (px)", 15, '100%')
        ),
        column(
            6,
            createTextSize("legendSizeGC", "Legend label size (px)",
                           12, '100%')
        ),
        hr(),
        column(
            12,
            HTML("<p><span style=\"color: #ff0000;\"><strong><em>Options
                     for variable plot only:</em></strong></span></p>")
        ),
        column(
            6, createPlotSize("widthVarGC", "Width (px)", 600)
        ),
        column(
            6, createPlotSize("heightVarGC", "Height (px)", 400)
        ),
        column(
            6,
            selectInput(
                "legendGC", label = "Legend position:",
                choices = list("Right" = "right",
                               "Bottom" = "bottom",
                               "Hide" = "none"),
                selected = "right",
                width = '100%'
            )
        ),
        column(
            6,
            selectInput(
                "mValueGC", label = "Show mean/median point:",
                choices = list("Mean" = "mean",
                               "Median" = "median"),
                selected = "mean",
                width = '100%')
        ),
        hr(),
        column(
            12,
            HTML("<p><span style=\"color: #ff0000;\"><strong><em>Options
                     for feature plot only:</em></strong></span></p>")
        ),
        column(
            6, createPlotSize("widthFeatureGC", "Width (px)", 600)
        ),
        column(
            6, createPlotSize("heightFeatureGC", "Height (px)", 400)
        ),
        column(
            6,
            radioButtons(
                inputId = "xAxisGC",
                label = "Flip coordinates",
                choices = list("Yes", "No"),
                selected = "No",
                inline = TRUE
            )
        ),
        column(
            6,
            sliderInput(
                "angleGC", "Angle of the X-axis label",
                min = 0,
                max = 180,
                step = 1,
                value = 60,
                width = 250
            )
        ),
        hr(),
        column(
            12,
            HTML("<p><span style=\"color: #ff0000;\"><strong><em>Names
                     of in-group and out-group taxa:</em></strong></span></p>")
        ),
        column(
            6,
            textInput("inGroupName",
                      NULL,
                      value = "In-group",
                      width = "100%",
                      placeholder = "Name of in-group taxa")
        ),
        column(
            6,
            textInput("outGroupName",
                      NULL,
                      value = "Out-group",
                      width = "100%",
                      placeholder = "Name of in-group taxa")
        ),
        br(),
        hr(),
        actionButton("resetConfigGC", "Reset", class = "btn btn-danger"),
    )
}

# * popup for input in-group and out-group on Group comparison ---------
gcUploadModal <- function() {
    modalDialog(
        title = "Upload files for group comparison",
        size = "l",
        h4(strong("Sequence(s) that need to be compared")),
        HTML("<p><em>Please use the same sequence IDs as being shown in the
                 phylogenetic profiles!</em></p>"),
        fileInput("gcFile", NULL),
        hr(),
        h4(strong("In-group and out-group taxa")),
        HTML("<p><em>Please upload an <strong>tab-delimited ID list</strong>
                 of in-group and out-group taxa. The ID must follow the format
                 \"<span style=\"text-decoration: underline;\">
                 <strong>ncbi12345</strong></span>\".
                 Example:</em></p>
                     <p><code>
                     ncbi12345 &nbsp; &nbsp;in-group</code></p>
                     <p><code>ncbi24242 &nbsp; &nbsp;in-group</code></p>
                     <p><code>ncbi33333 &nbsp; &nbsp;out-group
                 </code></p>"),
        
        fileInput("taxonGroupGC", NULL),
        conditionalPanel(
            condition = "output.checkTaxonGroupGC == false",
            h5(strong("Invalid taxa were found:")),
            DT::dataTableOutput("invalidTaxonGroupGC")
        )
    )
}


# * popup for selecting a whole taxon clade -------------------
taxaSelectModal <- function(id, title, checkboxId, checkboxLabel) {
    modalDialog(
        title = title,
        size = "s",
        easyClose = TRUE,
        footer = NULL,
        selectTaxonRankUI(id),
        checkboxInput(
            checkboxId,
            strong(checkboxLabel, style = "color:red"),
            value = FALSE
        )
    )
}

# * popup for selecting profile plot colors -------------------
colorSettingsModal <- function() {
    modalDialog(
        title = "Set colors for profile",
        size = "s",
        easyClose = TRUE,
        footer = NULL,
        
        colourpicker::colourInput(
            "lowColorVar1", "Low variable 1 (dot)", value = "#FF8C00"
        ),
        colourpicker::colourInput(
            "midColorVar1", "Mid variable 1 (dot)", value = "#40ABCF"
        ),
        colourpicker::colourInput(
            "highColorVar1", "High variable 1 (dot)", value = "#164294"
        ),
        numericInput(
            "midVar1", "Midpoint variable 1", 
            min = 0, max = 1, step = 0.01, value = 0.5
        ),
        actionButton(
            "defaultColorVar1", "Default", style = "padding:4px; font-size:100%"
        ),
        hr(),
        
        colourpicker::colourInput(
            "lowColorVar2", "Low variable 2 (background)", value = "#CC8D8D"
        ),
        colourpicker::colourInput(
            "midColorVar2", "Mid variable 2 (background)", value = "#FFFFFF"
        ),
        colourpicker::colourInput(
            "highColorVar2", "High variable 2 (background)", value = "#616587"
        ),
        numericInput(
            "midVar2", "Midpoint variable 2", 
            min = 0, max = 1, step = 0.01, value = 1
        ),
        actionButton(
            "defaultColorVar2", "Default", style = "padding:4px; font-size:100%"
        ),
        hr(),
        
        colourpicker::colourInput(
            "paraColor", "Color for inparalogs", value = "#07d000"
        ),
        actionButton(
            "defaultColorPara", "Default", style = "padding:4px; font-size:100%"
        )
    )
}


# * popup for selecting profile plot appearance (except colors) ----------------
plotConfigModal <- function(suffix = "", defaults, resetId, applyId) {
    modalDialog(
        title = "Plot appearance configuration",
        size = "m",
        easyClose = TRUE,
        footer = NULL,
        
        fluidRow(
            column(
                6, 
                createTextSize(
                    idWithSuffix("xSize", suffix), "X-axis label size (px)", 
                    defaults$xSize, 100
                )
            ),
            column(
                6, 
                createTextSize(
                    idWithSuffix("ySize", suffix), "Y-axis label size (px)", 
                    defaults$ySize, 100
                )
            ),
            column(
                6, 
                createTextSize(
                    idWithSuffix("legendSize", suffix), "Legend label size (px)",
                    8, 150
                )
            ),
            column(6, selectInput(
                idWithSuffix(
                    if (suffix=="") "mainLegend" else "selectedLegend", ""
                ),
                "Legend position:",
                choices = list(
                    "Right" = "right", "Left" = "left", "Top" = "top",
                    "Bottom" = "bottom", "Hide" = "none"
                ),
                selected = defaults$legend,
                width = 150
            )),
            
            column(
                6,
                createTextSize(
                    idWithSuffix("groupLabelSize", suffix), 
                    "Group label size (px)", 7, 100
                )
            ),
            column(
                6, 
                numericInput(
                    idWithSuffix("groupLabelDist", suffix), 
                    "Height for group label",
                    min = 1, max = 2, step = 0.1, 
                    value = defaults$groupLabelDist, width = 100
                )
            ),
            column(12,
                   HTML("<strong>Angle for taxonomic group label</strong>:<br>"),
                   sliderInput(
                       idWithSuffix("groupLabelAngle", suffix), "", 
                       min = 0, max = 90, step = 10, 
                       value = defaults$groupLabelAngle, width = 250
                   ),
                   br()
            ),
            column(12,
                   HTML("<strong>Angle for x-axis label</strong>:<br>"),
                   sliderInput(
                       idWithSuffix("xAngle", suffix), "", 
                       min = 0, max = 90, step = 10, 
                       value = defaults$xAngle, width = 250
                   ),
                   br()
            ),
            column(12,
                   HTML("<strong>Zooming factor (α) for dots on profile</strong>:<br>"),
                   sliderInput(
                       idWithSuffix("dotZoom", suffix), "", 
                       min = -1, max = 3, step = 0.1, 
                       value = defaults$dotZoom, width = 250
                   ),
                   HTML("<em>dot size = (1+α)*defaultSize<br>defaultSize=[0:5]</em>"),
                   uiOutput(idWithSuffix("dotSizeInfo", suffix)),
                   br()
            )
        ),
        br(), hr(),
        actionButton(resetId, "Reset", class = "btn btn-danger"),
    )
}
