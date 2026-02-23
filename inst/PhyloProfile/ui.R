#' Import function files
sourceFiles = list.files(path = "R", pattern = "*.R$", full.names = TRUE)
lapply(sourceFiles, source, .GlobalEnv)

#' MAIN UI =====================================================================
shinyUI(
    fluidPage(
        includeCSS("www/custom.css"),
        tags$style(type = "text/css", "body {padding-top: 80px;}"),
        shinyjs::useShinyjs(),
        use_bs_popover(),

        # Application title
        titlePanel("", windowTitle = "PhyloProfile"),

        # TOP WELLPANEL FOR PLOT CONFIGURATION ---------------------------------
        conditionalPanel(
            condition = "input.tabs=='Main profile'",
            wellPanel(
                fluidRow(
                    column(
                        2,
                        radioButtons(
                            inputId = "xAxis",
                            label = "Choose type of x-axis:",
                            choices = list("taxa", "genes"),
                            selected = "taxa",
                            inline = TRUE
                        ),
                        radioButtons(
                            inputId = "geneIdType",
                            label = "Display genes using:",
                            choices = list(
                                "Gene IDs" = "geneID","Gene names" = "geneName"
                            ),
                            selected = "geneID",
                            inline = TRUE
                        ),
                        checkboxInput(
                            "keepOrder",
                            strong(em("Retain gene order")),
                            value = TRUE,
                            width = NULL
                        ) %>%
                            bs_embed_popover(
                                title = "",
                                content = "Do not change gene order while filtering data",
                                placement = "bottom",
                                trigger = "hover"
                            )
                    ),
                    column(
                        1,
                        createPlotSize("width", "Width (px)", 900),
                        checkboxInput(
                            "autoSizing",
                            strong(em("Auto sizing")),
                            value = TRUE,
                            width = NULL
                        )
                    ),
                    column(
                        1, createPlotSize("height", "Height (px)", 900),
                        actionButton("mainPlotConfig", "Appearance")
                    ),
                    column(
                        2, uiOutput("var1Cutoff.ui")
                    ),
                    column(
                        2, uiOutput("var2Cutoff.ui")
                    ),
                    column(
                        2, uiOutput("percentCutoff.ui"),
                        actionButton(
                            "applyFilter",
                            label = "Apply filter",
                            icon = icon("check"),
                            class = "btn btn-warning"
                        )
                    ),
                    column(
                        2,
                        numericInput(
                            "coortholog",
                            "Max co-orthologs",
                            min = 1,
                            max = 999,
                            step = 1,
                            value = 999,
                            width = 150
                        ),
                        actionButton(
                            "resetMain",
                            label = "Reset cutoffs",
                            icon = icon("backward"),
                            class = "btn btn-danger"
                        )
                    )
                )
            )
        ),

        conditionalPanel(
            condition = "input.tabs=='Customized profile'",
            wellPanel(
                fluidRow(
                    column(
                        2,
                        radioButtons(
                            inputId = "xAxisSelected",
                            label = "Choose type of x-axis:",
                            choices = list("taxa", "genes"),
                            selected = "taxa",
                            inline = TRUE
                        ),
                        hr(),
                    ),
                    column(
                        1,
                        createPlotSize("selectedWidth", "Width (px)", 900),
                        checkboxInput(
                            "selectedAutoSizing",
                            strong(em("Auto sizing")),
                            value = TRUE,
                            width = NULL
                        )
                    ),
                    column(
                        1, createPlotSize("selectedHeight", "Height (px)", 900),
                        actionButton("selectedPlotConfig", "Appearance")
                    ),
                    column(
                        2, uiOutput("var1Filter.ui")
                    ),
                    column(
                        2, uiOutput("var2Filter.ui")
                    ),
                    column(
                        2, uiOutput("percentFilter.ui"),
                        actionButton(
                            "applyFilterCustom",
                            label = "Apply filter",
                            icon = icon("check"),
                            class = "btn btn-warning"
                        )

                    ),
                    column(
                        2,
                        uiOutput("coorthologFilter.ui"),
                        actionButton(
                            "resetSelected",
                            "Reset cutoffs",
                            icon = icon("backward"),
                            class = "btn btn-danger"
                        )
                    )
                )
            )
        ),

        # MAIN NARVARPAGE TABS -------------------------------------------------
        navbarPage(
            em(strong("PhyloProfile v2.2.3")),
            id = "tabs",
            collapsible = TRUE,
            inverse = TRUE,
            fluid = TRUE,
            position = "fixed-top",

            # INPUT TAB --------------------------------------------------------
            tabPanel(
                "Input & settings",
                # * 1st column -------------------------------------------------
                column(
                    4,
                    # ** Main input --------------------------------------------
                    strong(h4("Main input:")),
                    conditionalPanel(
                        condition = "input.do",
                        em(
                            strong(
                                "RELOAD THIS TOOL TO UPLOAD A NEW INPUT
                                FILE!!!", style = "color:red"
                            )
                        )
                    ),

                    selectInput(
                        "demoData", label = h5("Use online demo data:"),
                        choices = list(
                            "None" = "none",
                            "AMPK-TOR" = "ampk-tor",
                            "BUSCO Arthropoda" = "arthropoda"
                        ),
                        selected = "none",
                        width = "80%"
                    ),

                    uiOutput("noInternetMsg"),
                    uiOutput("demoDataDescribe"),
                    uiOutput("mainInputFile.ui"),
                    uiOutput("inputCheck.ui"),

                    fluidRow(
                        column(
                            6,
                            conditionalPanel(
                                condition = "output.checkOmaInput",
                                actionButton(
                                    "openOmaWindows", "Get data from OMA"
                                ),
                                br()
                            )
                        )
                    ),

                    # ** Variable 1 --------------------------------------------
                    fluidRow(
                        column(
                            4, uiOutput("var1ID.ui")
                        ),
                        column(
                            4,
                            selectInput(
                                "var1AggregateBy",
                                label = h5("Aggregate by:"),
                                choices = list(
                                    "Max" = "max",
                                    "Min" = "min",
                                    "Mean" = "mean",
                                    "Median" = "median"
                                ),
                                selected = "max",
                                width = 130
                            )
                        ),
                        column(
                            4,
                            selectInput(
                                "var1Relation", label = h5("Relationship:"),
                                choices = list(
                                    "Prot-Prot" = "protein",
                                    "Prot-Spec" = "species"
                                ),
                                selected = "protein",
                                width = 130
                            ) %>%
                                bs_embed_popover(
                                    title = "",
                                    content = paste(
                                        "select if variable is the value between",
                                        "seed protein - ortholog protein",
                                        "or seed protein - search taxon",
                                        sep = " "
                                    ),
                                    placement = "bottom",
                                    trigger = "hover"
                                )
                        )
                    ),

                    # ** Variable 2 --------------------------------------------
                    fluidRow(
                        column(
                            4, uiOutput("var2ID.ui")
                        ),
                        column(
                            4,
                            selectInput(
                                "var2AggregateBy",
                                label = h5("Aggregate by:"),
                                choices = list(
                                    "Max" = "max",
                                    "Min" = "min",
                                    "Mean" = "mean",
                                    "Median" = "median"
                                ),
                                selected = "max",
                                width = 130
                            )
                        ),
                        column(
                            4,
                            selectInput(
                                "var2Relation", label = h5("Relationship:"),
                                choices = list(
                                    "Prot-Prot" = "protein",
                                    "Prot-Spec" = "species"
                                ),
                                selected = "protein",
                                width = 130
                            )
                        )
                    ),
                    hr(),

                    # ** Domain input ------------------------------------------
                    strong(h4("Additional annotation input:")),
                    radioButtons(
                        inputId = "annoLocation", label = "",
                        choices = list("from file", "from folder"),
                        selected = "from file",
                        inline = TRUE
                    ),

                    uiOutput("domainInputFile.ui")
                ),

                # * 2nd column -------------------------------------------------
                column(
                    3,
                    uiOutput("fileExistMsgUI"),
                    uiOutput("inputMsgUI"),

                    # ** List of new taxa --------------------------------------
                    conditionalPanel(
                        condition = "output.unkTaxaStatus == 'unknown' ||
                        output.unkTaxaStatus == 'ncbi' ||
                        output.unkTaxaStatus == 'invalid'",
                        strong(h4("New taxa were found:")),
                        DT::dataTableOutput("unkTaxaFull"),
                        br(),
                        downloadButton("unkTaxa.download", "Download ID list")
                    ),

                    # ** Other input options -----------------------------------
                    conditionalPanel(
                        condition = "output.unkTaxaStatus == 0",
                        strong(h4("Choose genes of interest:")),
                        radioButtons(
                            inputId = "geneListSelected",
                            label = "",
                            choices = list("all", "from file"),
                            selected = "all",
                            inline = TRUE
                        ),

                        conditionalPanel(
                            condition = "input.geneListSelected == 'from file'",
                            fileInput("geneList", "")
                        ),
                        uiOutput("totalGeneNumber.ui"),
                        column(
                            6,
                            numericInput(
                                "stIndex",
                                h5("Show from:"),
                                min = 1,
                                max = 1600,
                                value = 1,
                                width = 130
                            ) %>%
                                bs_embed_popover(
                                    title = "",
                                    content = "Set start index for sequence range",
                                    placement = "bottom",
                                    trigger = "hover"
                                )
                        ),

                        column(
                            6,
                            numericInput(
                                "endIndex",
                                h5("...to:"),
                                min = 1,
                                max = 1600,
                                value = 1000,
                                width = 130
                            ) %>%
                                bs_embed_popover(
                                    title = "",
                                    content = "Set end index for sequence range",
                                    placement = "bottom",
                                    trigger = "hover"
                                )
                        ),
                        hr(),

                        strong(h4("Sequence source:")),
                        column(
                            6,
                            selectInput(
                                "seedSource",
                                label = h5("Seeds:"),
                                choices = list(
                                    "NCBI" = "ncbi",
                                    "UniProt" = "uniprot",
                                    "OrthoDB" = "orthodb",
                                    "OMA" = "oma",
                                    "User-defined" = "user"
                                ),
                                selected = "uniprot",
                                width = 130
                            ),
                            conditionalPanel(
                                condition = "input.seedSource == 'orthodb'",
                                textInput(
                                    "orthodbSeedVer",
                                    h5("OrthoDB version"),
                                    value = "",
                                    placeholder = "e.g. 10-1"
                                ) %>%
                                    bs_embed_popover(
                                        title = "",
                                        content = "Leave blank for latest version",
                                        placement = "bottom",
                                        trigger = "hover"
                                    )
                            )
                        ),
                        column(
                            6,
                            selectInput(
                                "orthoSource",
                                label = h5("Orthologs:"),
                                choices = list(
                                    "NCBI" = "ncbi",
                                    "UniProt" = "uniprot",
                                    "OrthoDB" = "orthodb",
                                    "OMA" = "oma",
                                    "User-defined" = "user"
                                ),
                                selected = "ncbi",
                                width = 130
                            ),
                            conditionalPanel(
                                condition = "input.orthoSource == 'orthodb'",
                                textInput(
                                    "orthodbOrthoVer",
                                    h5("OrthoDB version"),
                                    value = "",
                                    placeholder = "e.g. 10-1"
                                ) %>%
                                    bs_embed_popover(
                                        title = "",
                                        content = "Leave blank for latest version",
                                        placement = "bottom",
                                        trigger = "hover"
                                    )
                            )
                        ),
                        actionButton(
                            "selectSequenceID", "Ortholog ID format",
                            style = "padding:4px; font-size:100%"
                        ),
                        h5(""),
                        hr(),

                        strong(h4("Other optional input:")),

                        actionButton("fastaUpload", "FASTA file(s)"),
                        h5(""),

                        actionButton("uploadGeneCategory", "Gene categories"),
                        h5(""),

                        actionButton("uploadGeneName", "Gene names"),
                        hr(),

                        strong(h4("General configuration:")),
                        column (
                            6,
                            strong("Colors"),
                            br(), br(),
                            actionButton(
                                "setColor", "Change colors",
                                style = "padding:4px; font-size:100%"
                            )
                        ),
                        column(
                            6,
                            strong("Font"),
                            selectInput(
                                "font","",
                                choices = sort(unique(systemfonts::system_fonts()$family)),
                                selected = "Arial"
                            ) %>%
                                bs_embed_popover(
                                    title = "",
                                    content = "Note: This will be applied for all plots!",
                                    placement = "bottom",
                                    trigger = "hover"
                                )
                        ),
                        h5(""),
                        hr()
                    )
                ),

                # * 3rd column -------------------------------------------------
                column(
                    4,
                    # ** Location for taxonomy files ---------------------------
                    strong(h4("Taxonomy DB location:")),
                    radioButtons(
                        inputId = "taxDbLoc", label = "",
                        choices = list("Default", "User-defined"),
                        selected = "Default",
                        inline = TRUE
                    ),
                    conditionalPanel(
                        condition = "input.taxDbLoc == 'User-defined'",
                        shinyFiles::shinyDirButton(
                            "taxDbDir", "Select taxonomy DB" ,
                            title = "Please select a folder",
                            buttonType = "default", class = NULL
                        ),
                        br(), br(),
                        uiOutput("userTaxDBwarning")
                    ),
                    verbatimTextOutput("taxDbPath"),
                    hr(),

                    # ** Msg for parsing new taxa ------------------------------
                    conditionalPanel(
                        condition = "output.unkTaxaStatus == 'unknown' ||
                        output.unkTaxaStatus == 'ncbi' ||
                        output.unkTaxaStatus == 'invalid'",

                        conditionalPanel(
                            condition = "output.unkTaxaStatus == 'invalid'",
                            HTML(
                                "<p><em>Some new taxa have
                                <span style=\"color: #ff0000;\">invalid IDs
                                </span> (either in newTaxa.txt or in the main
                                profile input or both). IDs of non-NCBI taxa
                                have to be greater than 999999005.</em></p>
                                <p><em>Please replace those IDs before
                                continuing!</em></p>"
                            )
                        ),

                        conditionalPanel(
                            condition = "output.unkTaxaStatus == 'unknown'",
                            HTML(
                                '<p><em>NCBI taxonomy information of some taxa
                                can neither</em></p><ul><li><em>be retrieved
                                from NCBI (<span style="color: #0000ff;
                                ">Source="ncbi"</span>) nor </em></li>
                                <li><em>be found in
                                <span style="color: #ff0000;">
                                PhyloProfile/data/newTaxa.txt</span>&nbsp;(<span
                                style="color: #0000ff;">Source="new"</span>)
                                file</em></li></ul><p><strong><em>Please add
                                taxonomy information for those unknown taxa and
                                <span style="color: #ff0000;"> reload the tool
                                </span> to continue!</em></strong></p>'
                            ),
                            h5(""),
                            actionButton(
                                "addTaxa",
                                "Add taxonomy info",
                                icon = icon("plus"),
                                class = "btn btn-warning"
                            )
                        ),

                        conditionalPanel(
                            condition = "output.unkTaxaStatus == 'ncbi'",
                            HTML(
                                '<p><em>NCBI taxonomy information of some taxa
                                can either</em></p><ul><li><em>be retrieved
                                from NCBI (<span style="color: #0000ff;
                                ">Source="ncbi"</span>) or </em></li>
                                <li><em>be found in
                                <span style="color: #ff0000;">
                                PhyloProfile/data/newTaxa.txt</span>&nbsp;(<span
                                style="color: #0000ff;">Source="new"</span>)
                                file</em></li></ul><p><strong><em>Click here to
                                get required taxonomy information for those
                                taxa!</em></strong></p>'
                            ),
                            h5(""),
                            actionButton(
                                inputId = "butParse",
                                label = "Get taxonomy info",
                                class = "btn btn-warning",
                                disabled = FALSE
                            ),

                            hr(),
                            uiOutput("endParsingMsg"),
                            tableOutput("invalidID.output"),
                            hr(),
                            conditionalPanel(
                                condition = "output.unkTaxaStatus
                                            == 'invalid'",
                                downloadButton(
                                    "invalidID.download", "Download invalid IDs"
                                )
                            )
                        )
                    ),

                    # ** List of ranks & available taxa ------------------------
                    conditionalPanel(
                        condition = "output.unkTaxaStatus == 0",
                        strong(h4("Seed (super)taxon:")),
                        br(),

                        strong(h5("Select taxonomy rank:")),
                        uiOutput("rankSelect"),

                        strong(h5("Choose (super)taxon of interest:")),
                        selectizeInput(
                            "inSelect", "", choices = NULL, selected = NULL
                        ),
                        uiOutput("taxonSelectionMsg"),

                        hr(),

                        # ** Sort gene IDs options -----------------------------
                        strong(h4("Order seed IDs")),
                        radioButtons(
                            inputId = "orderGenes",
                            label = "",
                            choices = list(
                                "none", "alphabetically",
                                "profile similarity", "user defined"
                            ),
                            selected = "profile similarity",
                            inline = TRUE
                        ),
                        conditionalPanel(
                            condition = "input.orderGenes
                                        == 'user defined'",
                            uiOutput("inputSortedGenes.ui"),
                            uiOutput("checkSortedGenes.ui")
                        ),

                        # ** Sort taxa options ---------------------------------
                        strong(h4("Order taxa")),
                        radioButtons(
                            inputId = "orderTaxa",
                            label = "",
                            choices = list(
                                "automatically", "by user defined tree",
                                "by a sorted list"
                            ),
                            selected = "automatically",
                            inline = TRUE
                        ),
                        conditionalPanel(
                            condition = "input.orderTaxa
                                        == 'by user defined tree'",
                            uiOutput("inputTree.ui") %>%
                                bsplus::bs_embed_popover(
                                    title = "",
                                    content = "in newick format",
                                    placement = "bottom",
                                    trigger = "hover"
                                ),
                            uiOutput("checkNewick.ui")
                        ),
                        conditionalPanel(
                            condition = "input.orderTaxa
                                        == 'by a sorted list'",
                            uiOutput("inputSortedTaxa.ui"),
                            uiOutput("checkSortedTaxa.ui")
                        ),

                        checkboxInput(
                            "showAllTaxa",
                            strong("Display all input taxa"),
                            value = FALSE,
                            width = NULL
                        ) %>%
                            bsplus::bs_embed_popover(
                                title = "",
                                content = "Including taxa with no orthologs",
                                placement = "bottom",
                                trigger = "hover"
                            ),
                        hr(),

                        actionButton(
                            "do", "PLOT",
                            class = "btn btn-danger btn-lg"
                        ),
                        h5("")
                    )
                )
            ),

            # MAIN PROFILE TAB =================================================
            tabPanel(
                "Main profile",
                sidebarLayout(
                    # * sidebar panel for profile highlight --------------------
                    sidebarPanel(
                        width = 4,
                        column(
                            12,
                            style = "padding:0px;",
                            column(
                                4,
                                style = "padding:0px;",
                                radioButtons(
                                    inputId = "plotMode",
                                    label = "Plot mode:",
                                    choices = list(
                                        "Normal" = "normal", "Fast" = "fast"
                                    ),
                                    selected = "normal",
                                    inline = TRUE
                                )  %>%
                                    bsplus::bs_embed_popover(
                                        title = "",
                                        content = "Fast mode is only recommended for large data",
                                        placement = "top",
                                        trigger = "hover"
                                    )
                            ),
                            column(
                                8,
                                uiOutput("colorVar.ui")
                            )
                        ),
                        column(
                            12,
                            style = "padding:0px;",
                            strong("Select gene to highlight:")
                        ),
                        column(
                            12,
                            fluidRow(
                                column(
                                    8,
                                    style = "padding:0px;",
                                    selectizeInput(
                                        "geneHighlight","", NULL, multiple=TRUE,
                                        options = list(placeholder = 'none')
                                    )
                                ),
                                column(
                                    4,
                                    fileInput(
                                        "geneHighlightFile", "", width = "100%"
                                    )
                                )
                            )
                        ),
                        column(
                            12,
                            style = "padding:0px;",
                            strong(
                                "Select (super)taxon to highlight:"
                            )
                        ),
                        column(
                            8,
                            style = "padding:0px;",
                            selectizeInput(
                                "taxonHighlight","", NULL, multiple=TRUE,
                                options = list(placeholder = 'none')
                            )
                        ),
                        column(
                            4,
                            h3(""),
                            actionButton("taxonHighlightBrowse", "Browse...")
                        ),
                        column(
                            12,
                            checkboxInput(
                                "colorByGroup",
                                strong("Highlight genes by categories"),
                                value = FALSE
                            ),
                            checkboxInput(
                                "colorByOrthoID",
                                strong("Highlight duplicated ortholog IDs"),
                                value = FALSE
                            ) %>%
                                bsplus::bs_embed_popover(
                                    title = "",
                                    content = paste(
                                        "Please check in the Clustering profiles",
                                        "function, if the profiles are clustered",
                                        "using ortho IDs"
                                    ),
                                    placement = "bottom",
                                    trigger = "hover"
                                ),
                            hr()
                        ),
                        uiOutput("superRankSelect.ui"),
                        hr(),
                        checkboxInput(
                            "autoUpdate",
                            strong(em("Auto update plot")),
                            value = FALSE,
                            width = NULL
                        ),
                        actionButton(
                            "updateBtn", "Update appearance",
                            icon = icon("sync"), class = "btn btn-warning"
                        )
                    ),
                    # * main panel for profile plot ----------------------------
                    mainPanel(
                        conditionalPanel(
                            condition = "input.do > 0",
                            createProfilePlotUI("mainProfile")
                        )
                    )
                )
            ),

            # CUSTOMIZED PROFILE TAB ===========================================
            tabPanel(
                "Customized profile",
                sidebarLayout(
                    # * sidebar panel for subseting data -----------------------
                    sidebarPanel(
                        width = 4,
                        column(
                            12,
                            style = "padding:0px;",
                            strong("Select sequence(s) of interest:")
                        ),

                        column(
                            12,
                            fluidRow(
                                column(
                                    8,
                                    style = "padding:0px;",
                                    selectizeInput(
                                        "inSeq","", NULL, multiple=TRUE,
                                        choices = c('all'), selected = "all"
                                    )
                                ),
                                column(
                                    4,
                                    fileInput("customFile", "", width = "100%")
                                )
                            )
                        ),

                        column(
                            12,
                            style = "padding:0px;",
                            strong(
                                "Select (super)taxon/(super)taxa of interest:"
                            )
                        ),
                        column(
                            12,
                            fluidRow(
                                column(
                                    8,
                                    style = "padding:0px;",
                                    uiOutput("cusTaxa.ui")
                                ),
                                column(
                                    4,
                                    h3(""),
                                    actionButton("cusTaxa", "Browse...")
                                )
                            )
                        ),
                        uiOutput("cusSuperRankSelect.ui"),

                        h5(""),
                        actionButton(
                            "plotCustom", "Update apperance",
                            icon = icon("sync"), class = "btn btn-warning"
                        )
                    ),

                    # * main panel for customized profile plot -----------------
                    mainPanel(
                        conditionalPanel(
                            condition = "output.sameProfile == true",
                            h4(
                                "Please select subset of genes and/
                                or taxa for customized profile!"
                            )
                        ),

                        conditionalPanel(
                            condition = "input.do > 0 | input.plotCustom > 0",
                            createProfilePlotUI("customizedProfile")
                        )
                    )
                )
            ),

            # DIMENSIONALITY REDUCTION reduction TAB ===========================
            tabPanel(
                "Dimension Reduction",
                # * Top panel for plot configuration ---------------------------
                wellPanel(
                    fluidRow(
                        column(
                            2,
                            column(
                                6,
                                style = "padding:0px;",
                                selectInput(
                                    "reductionTechnique", label = "Technique",
                                    choices = list("UMAP" = "umap",
                                                   "t-SNE" = "tsne",
                                                   "PCA" = "pca"),
                                    selected = "umap"
                                )
                            ),
                            column(
                                6,
                                conditionalPanel(
                                    condition = 'input.reductionTechnique == "tsne"',
                                    numericInput(
                                        "tsneIter", "# iter.",
                                        min = 100, max = 5000, step = 100,
                                        value = 1000
                                    )
                                )
                            ),
                            column(
                                12,
                                style = "padding:0px;",
                                radioButtons(
                                    "dimRedType", "Apply on",
                                    c("Taxa" = "taxa", "Genes" = "genes"),
                                    inline = TRUE
                                ),
                                selectInput(
                                    "dimRedDataType", label = "using",
                                    choices = list(
                                        "Presence/Absence data" = "binary",
                                        "Numeric scores" = "nonbinary"
                                    ),
                                    selected = "binary"
                                )
                            )
                        ),
                        column(
                            4,
                            column(
                                3,
                                createPlotSize(
                                    "dimRedPlot.width", "Plot width", 900
                                ),
                                createPlotSize(
                                    "dimRedPlot.height", "Plot height", 400
                                )
                            ),
                            column(
                                4,
                                createTextSize(
                                    "dimRedPlot.textsize", "Text size", 12
                                ),
                                selectInput(
                                    "dimRedPlot.legend",
                                    label = "Legend position:",
                                    choices = list("Right" = "right",
                                                   "Left" = "left",
                                                   "Top" = "top",
                                                   "Bottom" = "bottom",
                                                   "Hide" = "none"),
                                    selected = "bottom"
                                )
                            ),
                            column(
                                5,
                                sliderInput(
                                    "dimRedPlot.dotzoom", "Dot size zooming",
                                    min = -3, max = 10, step = 1, value = 0
                                )
                            )
                        ),
                        column(
                            4,
                            column(
                                6,
                                radioButtons(
                                    "dimRedGroupLabelsBy",
                                    "Summarize as [Other] by",
                                    choices = c("taxa", "genes"),
                                    inline = TRUE
                                ),
                                em(paste("If frequency smaller or higher than",
                                         "Freq cutoff, labels will be grouped",
                                         "as [Other]"))
                            ),
                            column(
                                6,
                                sliderInput(
                                    "dimRedLabelNr", "Freq cutoff", min = 0,
                                    max = 99, step = 1, value = c(5,99),
                                    width = 200
                                ),
                                radioButtons(
                                    "dimRedPlotType", "Plot dimensions",
                                    choices = list(
                                        "2D" = "ggplot", "3D" = "plotly"
                                    ),
                                    selected = "ggplot", inline = TRUE
                                )
                            )
                        ),
                        column(
                            2,
                            sliderInput(
                                "dimRedDotAlpha", "Transparent level", min = 0,
                                max = 1, step = 0.05, value = 0.0, width = 200
                            ),
                            numericInput(
                                "randomSeed", "Random seed",
                                min = 0, max = 9999, step = 1, value = 123
                            )
                        )
                    )
                ),
                sidebarLayout(
                    # * Sidebar panel for data filter --------------------------
                    sidebarPanel(
                        tableOutput("dimRedHoverInfo"),
                        hr(),
                        selectInput(
                            "dimRedRank", label = "Taxonomy rank for labels",
                            choices = getTaxonomyRanks(),
                            selected = "phylum"
                        ),
                        hr(),
                        textInput(
                            "dimRedGroupHigherRank",
                            "Group labels into higher rank",
                            value = "",
                            placeholder = paste(
                                "Type taxon names in higher rank",
                                "(e.g.: Fungi;Metazoa)"
                            )
                        ),
                        uiOutput("dimRedGroupHigherRank.warning"),
                        uiOutput("dimRedCustomLabel.ui"),
                        actionButton(
                            "dimRedApplyChangeLables",
                            "Change labels",
                            icon = icon("play"),
                            class = "btn btn-success"
                        ),
                        actionButton(
                            "dimRedResetLables",
                            "Reset labels",
                            icon = icon("rotate-left"),
                            class = "btn btn-secondary"
                        ) %>%
                            bsplus::bs_embed_popover(
                                title = "",
                                content = "Click `Change labels` after reset!",
                                placement = "bottom",
                                trigger = "hover"
                            ),
                        hr(),
                        column(
                            12,
                            style = "padding:0px;",
                            strong("Choose labels to hide")
                        ),
                        column(
                            12,
                            fluidRow(
                                column(
                                    8,
                                    style = "padding:0px;",
                                    selectInput(
                                        "excludeDimRedTaxa", label = "",
                                        choices = "", selected = NULL,
                                        multiple = TRUE
                                    )
                                ),
                                column(
                                    4,
                                    fileInput(
                                        "excludeDimRedTaxaFile", "", width = "100%"
                                    )
                                )
                            )
                        ),
                        column(
                            12,
                            style = "padding:0px;",
                            strong("Choose labels to highlight")
                        ),
                        column(
                            12,
                            fluidRow(
                                column(
                                    8,
                                    style = "padding:0px;",
                                    selectInput(
                                        "highlightDimRedTaxa", label = "",
                                        choices = "", selected = NULL,
                                        multiple = TRUE
                                    )
                                ),
                                column(
                                    4,
                                    fileInput(
                                        "highlightDimRedTaxaFile", "", width = "100%"
                                    )
                                )
                            )
                        ),
                        hr(),
                        selectInput(
                            "colorPalleteDimRed",
                            "Color pallete",
                            choices = c(
                                "Paired", "Set1", "Set2", "Set3",
                                "Accent", "Dark2"
                            ),
                            selected = "Dark2"
                        ),
                        hr(),
                        selectInput(
                            "dimRedFilterVar",
                            "Choose variable for data filtering",
                            choices = c("Var1" = "var1", "Var2" = "var2",
                                        "Both" = "both"),
                            selected = "both"
                        ),
                        sliderInput(
                            "dimRedCutoff", "Filter cutoff", min = 0, max = 1,
                            step = 0.05, value = 0, width = '100%'
                        ),
                        hr(),
                        strong("Add following data to Customized profile"),
                        checkboxInput("addSpecDimRed", em("Selected taxa")) %>%
                            bsplus::bs_embed_popover(
                                title = "",
                                content = paste("Only available when working with the",
                                                "lowest taxonomy rank!"),
                                placement = "bottom",
                                trigger = "hover"
                            ),
                        checkboxInput("addGeneDimRed", em("Selected genes")),
                        uiOutput("addDimRedCustomProfileCheck.ui")
                    ),
                    # * Main panel for DIM red plot and tables -----------------
                    mainPanel(
                        conditionalPanel(
                            condition = 'input.dimRedPlotType == "ggplot"',
                            column(
                                12,
                                em(
                                    "Drag to select, then double click to zoom in. Double click empty space to zoom out",
                                    style = "color:darkblue"
                                )
                            )
                        ),
                        conditionalPanel(
                            condition = 'input.dimRedPlotType == "plotly"',
                            column(
                                12,
                                em(
                                    "Click on data point to select",
                                    style = "color:darkblue"
                                )
                            )
                        ),
                        uiOutput("dimRedPlot.ui"),
                        br(),
                        column(
                            10,
                            column(
                                3,
                                actionButton(
                                    inputId = "plotDimRed",
                                    label = "PLOT",
                                    class = "btn btn-danger"
                                )
                            ),
                            conditionalPanel(
                                condition = 'input.dimRedPlotType == "plotly"',
                                column(
                                    3,
                                    actionButton(
                                        "clear", "CLEAR selected points"
                                    )
                                )
                            ),
                            column(
                                3,
                                conditionalPanel(
                                    condition = 'input.dimRedPlotType == "plotly"',

                                    downloadButton(
                                        "dimRedDownloadPlot3D", "Download 3D plot",
                                        class = "butDL"
                                    )
                                ),
                                conditionalPanel(
                                    condition = 'input.dimRedPlotType == "ggplot"',

                                    downloadButton(
                                        "dimRedDownloadPlot", "Download 2D plot",
                                        class = "butDL"
                                    )
                                )
                            ),
                            column(
                                3,
                                downloadButton(
                                    "dimRedDownloadData", "Download data",
                                    class = "butDL"
                                ) %>%
                                    bsplus::bs_embed_popover(
                                        title = "",
                                        content = paste("Use plotDimRed() to manually create",
                                                        "dimension reduction plot!"),
                                        placement = "bottom",
                                        trigger = "hover"
                                    )
                            )
                        ),
                        br(),
                        uiOutput("dimRedTable.ui")
                    )
                )
            ),

            # FUNCTION TAB =====================================================
            navbarMenu(
                "Function",
                # * Profile clustering ----------------------------------------
                tabPanel(
                    "Profile clustering", value = "profile_clustering",
                    h4(strong("Profile clustering")),
                    uiOutput("descClusteringUI"),

                    wellPanel(
                        fluidRow(
                            column(
                                2, uiOutput("selectProfileType")
                            ),
                            column(
                                3, uiOutput("selectDistMethod")
                            ),

                            column(
                                3,
                                selectInput(
                                    "clusterMethod",
                                    label = h5("Cluster method:"),
                                    choices = list(
                                        "single" = "single",
                                        "complete" = "complete",
                                        "average (UPGMA)" = "average",
                                        "mcquitty (WPGMA)" = "mcquitty",
                                        "median (WPGMC)" = "median",
                                        "centroid (UPGMC)" = "centroid"
                                    ),
                                    selected = "complete"
                                )
                            ),

                            column(
                                1,
                                createPlotSize(
                                    "clusterPlot.width", "Width (px)", 600
                                )
                            ),
                            column(
                                1,
                                createPlotSize(
                                    "clusterPlot.height", "Height (px)", 600
                                )
                            ),
                            column(
                                2,
                                checkboxInput(
                                    "applyCluster",
                                    em(strong(
                                        "Apply clustering to profile plot",
                                        style = "color:darkblue"
                                    )),
                                    value = TRUE
                                ) %>%
                                    bsplus::bs_embed_popover(
                                        title = "",
                                        content = "Uncheck this to sort genes by alphabet",
                                        placement = "bottom",
                                        trigger = "hover"
                                    ),

                                uiOutput("applyClusterCheck.ui"),

                                checkboxInput(
                                    "addClusterCustomProfile",
                                    strong(em(
                                        "Add selected genes to Customized
                                        profile",
                                        style = "color:red"
                                    )),
                                    value = FALSE,
                                    width = NULL
                                ),
                                uiOutput("addClusterCustomProfileCheck.ui")
                            )
                        )
                    ),

                    clusterProfileUI("profileClustering")
                ),

                # * Distribution analysis --------------------------------------
                tabPanel(
                    "Distribution analysis", value = "distribution_analysis",
                    h4(strong("Distribution analysis")),
                    uiOutput("descDistributionUI"),

                    wellPanel(
                        fluidRow(
                            column(
                                2,
                                selectInput(
                                    "dataset.distribution", "Select data",
                                    choices = c("Main data", "Customized data"),
                                    selected = "Main data"
                                ),
                                uiOutput("selected.distribution")
                            ),
                            column(
                                2, uiOutput("var1Dist.ui")
                            ),
                            column(
                                2, uiOutput("var2Dist.ui")
                            ),
                            column(
                                2, uiOutput("percentDist.ui")
                            ),
                            column(
                                2,
                                createTextSize(
                                    "distTextSize", "Label size", 12, 100
                                )
                            ),
                            column(
                                2,
                                createPlotSize(
                                    "distWidth", "Width (px)", 600
                                )
                            )
                        )
                    ),
                    analyzeDistributionUI("distPlot")
                ),


                # * Gene age estimation ----------------------------------------
                tabPanel(
                    "Gene age estimation", value = "gene_age_estimation",
                    h4(strong("Gene age estimation")),
                    uiOutput("descGeneAgeUI"),

                    wellPanel(
                        fluidRow(
                            column(
                                2, uiOutput("var1Age.ui")
                            ),
                            column(
                                2, uiOutput("var2Age.ui")
                            ),
                            column(
                                2, uiOutput("percentAge.ui")
                            ),
                            column(
                                2,
                                strong("Appearance"),
                                actionButton("geneAgeProtConfig", "Plot config")
                            ),
                            column(
                                4,
                                checkboxInput(
                                    "addGeneAgeCustomProfile",
                                    strong(em(
                                        "Add selected genes to Customized
                                        profile",
                                        style = "color:red"
                                    )),
                                    value = FALSE,
                                    width = NULL
                                ),
                                uiOutput("addGeneAgeCustomProfileCheck.ui")
                            )
                        )
                    ),
                    plotGeneAgeUI("geneAge")
                ),

                # * Core gene identification  ----------------------------------
                tabPanel(
                    "Core gene identification", value = "core_gene_ident",
                    h4(strong("Core gene identification")),
                    uiOutput("descCoreGeneUI"),

                    wellPanel(
                        fluidRow(
                            column(
                                3, uiOutput("var1Core.ui")
                            ),
                            column(
                                3, uiOutput("var2Core.ui")
                            ),
                            column(
                                3, uiOutput("percentCore.ui")
                            ),
                            column(
                                3,
                                sliderInput(
                                    "coreCoverage",
                                    "Core taxa coverage",
                                    min = 0,
                                    max = 100,
                                    value = 90,
                                    step = 5
                                )
                            ),
                            column(
                                12,
                                uiOutput("taxaListCore.ui"),
                                actionButton("browseTaxaCore", "Browse")
                            )
                        )
                    ),
                    hr(),

                    column(
                        4,
                        downloadButton(
                            "coreGeneTableDownload", "Download gene list"
                        ),
                        checkboxInput(
                            "addCoreGeneCustomProfile",
                            strong(em("Add core genes to Customized profile",
                                      style = "color:red")),
                            value = FALSE,
                            width = NULL
                        ),
                        uiOutput("addCoreGeneCustomProfileCheck.ui")
                    ),
                    identifyCoreGeneUI("coreGene")
                ),

                # * Group Comparison  ------------------------------------------
                tabPanel(
                    "Group comparison", value = "group_comparison",
                    h4(strong("Group comparison")),
                    uiOutput("descGCUI"),
                    wellPanel(
                        fluidRow(
                            column(
                                3,
                                uiOutput("taxaListGC"),
                                checkboxInput(
                                    "useCommonAncestor",
                                    "Use common ancestor",
                                    value = TRUE,
                                    width = NULL
                                ) %>%
                                    bsplus::bs_embed_popover(
                                        title = "",
                                        content = paste(
                                            "All taxa that have the same common",
                                            "ancestor with the selected taxa above",
                                            "will be considered as the in-group"
                                        ),
                                        placement = "bottom",
                                        trigger = "hover"
                                    ),
                                h5(strong("Select in-group by supertaxon")),
                                actionButton("taxaGC", "Browse")
                            ),
                            column(
                                3,
                                uiOutput("listGenesGC"),
                                h5(strong(
                                    "Upload sequence(s) / In-group taxa"
                                )),
                                actionButton("uploadGC", "Upload")
                            ),
                            column(
                                2,
                                uiOutput("variableGC"),
                                selectInput(
                                    "compareType", "Compare using:",
                                    choices = c(
                                        "Statistical tests", "Median values"
                                    ),
                                    selected = "Statistical tests"
                                )
                            ),
                            column(
                                2,
                                sliderInput(
                                    "significance",
                                    label = "Significance level:",
                                    min = 0,
                                    max = 1,
                                    step = 0.01,
                                    value = 0.05,
                                    width = 300
                                ) %>%
                                    bsplus::bs_embed_popover(
                                        title   = "",
                                        content = "P-value cut-off of the statistic test,
                                                OR cut-off of delta means between 2 groups",
                                        placement = "bottom",
                                        trigger = "hover"
                                    ),
                                checkboxInput(
                                    "addGCGenesCustomProfile",
                                    strong(em(
                                        "Add candidate gene(s) to Customized
                                        profile",
                                        style = "color:red"
                                    )),
                                    value = FALSE,
                                    width = NULL
                                ),
                                uiOutput("addGCCustomProfileCheck")
                            ),
                            column(
                                2,
                                actionButton("gcPlotConfig", "Plot config") %>%
                                    bsplus::bs_embed_popover(
                                        title   = "",
                                        content = "Change the appearance of the plots",
                                        placement = "bottom",
                                        trigger = "hover"
                                    ),
                                hr(),
                                actionButton(
                                    "doCompare", "COMPARE!",
                                    class = "btn btn-danger"
                                ),
                                conditionalPanel(
                                    condition = 'input.rankSelect == "strain"',
                                    HTML('<div style="color: red; font-weight: bold;">
                                      Warning: "strain" is not a supported taxonomy rank
                                      for this function.<br>
                                      Please choose another rank.
                                    </div>')
                                ),
                                h5(),
                                actionButton(
                                    "updateGC",
                                    "Update plot",
                                    icon = icon("sync"),
                                    class = "btn btn-warning"
                                )
                            )
                        )
                    ),
                    groupComparisonUI("groupComparison")
                ),

                # * NCBI taxonomy data -----------------------------------------
                tabPanel(
                    "NCBI taxonomy data", value = "ncbi_tax_data",
                    uiOutput("descNcbiTaxDbUI"),
                    column(
                        3,
                        radioButtons(
                            inputId = "taxDB",
                            label = "Choose a task:",
                            choices = list(
                                "Update NCBI taxonomy DB" = "update",
                                "Reset NCBI taxonomy DB" = "reset",
                                "Export taxonomy DB files" = "export",
                                "Import taxonomy DB files" = "import"
                            )
                        )
                    ),
                    column(
                        9,
                        conditionalPanel(
                            condition = "input.taxDB=='update'",
                            h4(strong("Update NCBI taxonomy")),
                            actionButton(
                                "doUpdateNcbi",
                                "Do update",
                                icon = icon("wrench"),
                                class = "btn btn-warning"
                            ),
                            hr(),
                            verbatimTextOutput("updateNCBITaxStatus")
                        ),
                        conditionalPanel(
                            condition = "input.taxDB=='reset'",
                            h4(strong("Reset taxonomy data")),
                            uiOutput("taxResetWarning.ui"),
                            br(),
                            actionButton(
                                "doResetTax",
                                "Do reset",
                                icon = icon("wrench"),
                                class = "btn btn-warning"
                            ),
                            hr(),
                            verbatimTextOutput("resetTaxonomyDataStatus")
                        ),
                        conditionalPanel(
                            condition =
                                "input.taxDB=='export'",
                            h4(strong("Export current taxonomy files")),
                            uiOutput("taxExportWarning.ui"),
                            br(),
                            shinyFiles::shinyDirButton(
                                "taxDirOut",
                                "Select output directory" ,
                                title = paste(
                                    "Please select output directory"
                                ),
                                buttonType = "default", class = NULL
                            ),
                            br(),
                            uiOutput("taxDirOut.ui"),
                            br(),
                            actionButton(
                                "doExportTax",
                                "Do export",
                                icon = icon("file-export"),
                                class = "btn btn-warning",
                            ),
                            hr(),
                            verbatimTextOutput("exportTaxonomyDataStatus")
                        ),
                        conditionalPanel(
                            condition =
                                "input.taxDB=='import'",
                            h4(strong("Import your own taxonomy files")),
                            uiOutput("taxImportWarning.ui"),
                            br(),
                            shinyFiles::shinyDirButton(
                                "taxDir",
                                "Select input directory" ,
                                title = paste(
                                    "Please select directory that contains
                                    the taxonomy files"
                                ),
                                buttonType = "default", class = NULL
                            ),
                            br(),
                            uiOutput("taxDir.ui"),
                            br(),
                            actionButton(
                                "doImportTax",
                                "Do import",
                                icon = icon("file-import"),
                                class = "btn btn-warning"
                            ),
                            hr(),
                            verbatimTextOutput("importTaxonomyDataStatus")
                        )
                    )
                )
            ),

            # DATA DOWNLOAD TAB ================================================
            navbarMenu(
                "Export data",

                # * Export data ------------------------------------------------
                downloadFilteredMainUI("filteredMainDownload"),
                downloadFilteredCustomizedUI("filteredCustomizedDownload"),

                # * Export processed data --------------------------------------
                tabPanel(
                    "Processed data",
                    h4(strong("Download processed data")),
                    uiOutput("descDownloadProcessedDataUI"),
                    strong("Output dir:"),
                    br(), br(),
                    shinyFiles::shinyDirButton(
                        "procDataOutDir",
                        "Select output directory" ,
                        title = paste(
                            "Please select output directory"
                        ),
                        buttonType = "default", class = NULL
                    ),
                    br(),br(),
                    uiOutput("procDataOutDir.ui"),
                    br(),
                    actionButton(
                        "doDownloadProcData",
                        "Download",
                        icon = icon("file-export"),
                        class = "btn btn-warning"
                    ),
                    hr(),
                    verbatimTextOutput("downloadProcDataStatus")
                ),

                # * Export plot settings ---------------------------------------
                tabPanel(
                    "Export plot settings",
                    h4(strong("Export plot settings")),
                    uiOutput("descExportSettingUI"),
                    radioButtons(
                        inputId = "exportSetting",
                        label = "Export as:",
                        choices = list(
                            "a list" = "list",
                            "an Rscript" = "rscript"
                        )
                    ),
                    hr(),
                    strong("Output dir:"),
                    br(), br(),
                    shinyFiles::shinyDirButton(
                        "settingDir",
                        "Select output directory" ,
                        title = paste(
                            "Please select output directory"
                        ),
                        buttonType = "default", class = NULL
                    ),
                    br(), br(),
                    strong("File name:"),
                    uiOutput("settingFile.ui"),
                    uiOutput("settingDir.ui"),
                    br(),
                    actionButton(
                        "doExportSetting",
                        "Do export",
                        icon = icon("file-export"),
                        class = "btn btn-warning"
                    ),
                    hr(),
                    verbatimTextOutput("exportSettingStatus")
                )
            ),

            # HELP TAB =========================================================
            navbarMenu(
                "Help",
                tabPanel(
                    a(
                        "Wiki",
                        href = "https://github.com/BIONF/PhyloProfile/wiki",
                        target = "_blank"
                    )
                ),
                tabPanel(
                    a(
                        "About",
                        href = "https://BIONF.github.io/PhyloProfile/",
                        target = "_blank"
                    )
                )
            )
        ),

        # POINT INFO BOX =======================================================
        conditionalPanel(
            condition =
                "input.tabs=='Main profile' ||
                input.tabs=='Customized profile'",
            absolutePanel(
                bottom = 5, left = 30,
                fixed = TRUE,
                draggable = TRUE,
                h5("Point's info:"),
                verbatimTextOutput("pointInfo"),
                conditionalPanel(
                    condition = "output.pointInfoStatus == 0",
                    actionButton(
                        "detailedBtn",
                        "Detailed plot",
                        class = "btn btn-success"
                    ),
                    actionButton(
                        "doDomainPlotMain",
                        "Domain plot",
                        class = "btn btn-success"
                    ),
                    actionButton(
                        "highlightMain",
                        "Highlight",
                        class = "btn btn-success"
                    )
                ),
                style = "opacity: 0.80"
            )
        )
    )
)
