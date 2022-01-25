#' Import function files
sourceFiles = list.files(path = "R", pattern = "*.R$", full.names = TRUE)
lapply(sourceFiles, source, .GlobalEnv)

#' MAIN UI ====================================================================
shinyUI(
    fluidPage(

        tags$style(type = "text/css", "body {padding-top: 80px;}"),
        useShinyjs(),

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
                        hr(),
                        checkboxInput(
                            "autoUpdate",
                            strong(em("Auto update plot")),
                            value = FALSE,
                            width = NULL
                        )
                    ),
                    column(
                        1,
                        createPlotSize("width", "Width (px)", 600),
                        checkboxInput(
                            "autoSizing",
                            strong(em("Auto sizing")),
                            value = TRUE,
                            width = NULL
                        )
                    ),
                    column(
                        1, createPlotSize("height", "Height (px)", 600),
                        actionButton("mainPlotConfig", "Appearance")
                    ),
                    column(
                        2, uiOutput("var1Cutoff.ui")
                    ),
                    column(
                        2, uiOutput("var2Cutoff.ui")
                    ),
                    column(
                        2, uiOutput("percentCutoff.ui")
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
                        bsButton(
                            "resetMain",
                            "Reset cutoffs",
                            style = "danger",
                            icon = icon("backward")
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
                        checkboxInput(
                            "autoUpdateSelected",
                            strong(em("Auto update plot")),
                            value = TRUE,
                            width = NULL
                        )
                    ),
                    column(
                        1,
                        createPlotSize("selectedWidth", "Width (px)", 600),
                        checkboxInput(
                            "selectedAutoSizing",
                            strong(em("Auto sizing")),
                            value = TRUE,
                            width = NULL
                        )
                    ),
                    column(
                        1, createPlotSize("selectedHeight", "Height (px)", 600),
                        actionButton("selectedPlotConfig", "Appearance")
                    ),
                    column(
                        2, uiOutput("var1Filter.ui")
                    ),
                    column(
                        2, uiOutput("var2Filter.ui")
                    ),
                    column(
                        2, uiOutput("percentFilter.ui")
                    ),
                    column(
                        2,
                        uiOutput("coorthologFilter.ui"),
                        bsButton(
                            "resetSelected",
                            "Reset cutoffs",
                            style = "danger",
                            icon = icon("backward")
                        )
                    )
                )
            )
        ),

        # MAIN NARVARPAGE TABS -------------------------------------------------
        navbarPage(
            em(strong("PhyloProfile v1.8.5")),
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
                                bsButton("openOmaWindows", "Get data from OMA"),
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
                            ),
                            bsPopover(
                                "var1Relation",
                                "",
                                paste(
                                    "select if variable is the value between ",
                                    "seed protein - ortholog protein",
                                    " or seed protein - search taxon",
                                    sep = "<br>"
                                ),
                                "bottom"
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
                                    "Prot-Prot" = "protein", "Prot-Spec" = "species"
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

                    uiOutput("domainInputFile.ui"),

                    hr(),
                    uiOutput("downloadDemo.ui")
                ),

                # * 2nd column -------------------------------------------------
                column(
                    3,
                    bsAlert("fileExistMsgUI"),
                    bsAlert("inputMsgUI"),

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
                            fileInput("list", "")
                        ),

                        hr(),
                        checkboxInput(
                            "ordering",
                            strong("Order sequence IDs"),
                            value = TRUE
                        ),

                        hr(),
                        HTML("<b>Order taxa</b>"),

                        radioButtons(
                            inputId = "orderTaxa",
                            label = "",
                            choices = list(
                                "automatically", "by user defined tree"
                            ),
                            selected = "automatically",
                            inline = TRUE
                        ),

                        bsPopover(
                            "orderTaxa", "", "in newick format", "bottom"
                        ),

                        conditionalPanel(
                            condition = "input.orderTaxa
                                        == 'by user defined tree'",
                            uiOutput("inputTree.ui")
                        ),

                        uiOutput("checkNewick.ui"),
                        hr(),

                        strong(h4("Other optional input:")),

                        bsButton("fastaUpload", "FASTA file(s)"),
                        h5(""),

                        bsButton("uploadGeneCategory", "Gene categories"),
                        h5(""),
                        hr(),

                        strong(h4("Color configuration:")),
                        actionButton(
                            "setColor", "Change colors",
                            style = "padding:4px; font-size:100%"
                        ),
                        hr()
                    )
                ),

                # * 3rd column -------------------------------------------------
                column(
                    4,
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
                            bsButton(
                                "addTaxa",
                                "Add taxonomy info",
                                disabled = FALSE,
                                style = "warning"
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
                            bsButton(
                                "butParse",
                                "Get taxonomy info",
                                disabled = FALSE,
                                style = "warning"
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
                        br(),

                        strong(h5("Choose (super)taxon of interest:")),
                        uiOutput("select"),
                        br(),

                        bsButton(
                            "do",
                            "PLOT",
                            type = "action",
                            style = "danger",
                            size = "large",
                            disabled = FALSE
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
                        uiOutput("totalGeneNumber.ui"),

                        column(
                            4,
                            numericInput(
                                "stIndex",
                                "Show from:",
                                min = 1,
                                max = 1600,
                                value = 1,
                                width = 100
                            ),
                            style = "padding:0px;"
                        ),

                        column(
                            4,
                            numericInput(
                                "endIndex",
                                "...to:",
                                min = 1,
                                max = 1600,
                                value = 1000,
                                width = 100
                            ),
                            style = "padding:0px;"
                        ),

                        column(
                            4, uiOutput("highlightGeneUI")
                        ),

                        bsPopover(
                            "highlightGeneUI",
                            "",
                            "Select gene to highlight",
                            "bottom"
                        ),

                        bsPopover(
                            "stIndex",
                            "",
                            "Set start index for sequence range",
                            "bottom"
                        ),

                        bsPopover(
                            "endIndex",
                            "",
                            "Set end index for sequence range",
                            "bottom"
                        ),

                        br(),
                        uiOutput("highlightTaxonUI"),

                        checkboxInput(
                            "colorByGroup",
                            strong("Highlight genes by categories"),
                            value = FALSE
                        ),

                        conditionalPanel(
                            condition = "input.autoUpdate == false",
                            bsButton(
                                "updateBtn",
                                "Update plot",
                                style = "warning",
                                icon("refresh")
                            )
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
                                    uiOutput("geneIn")
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
                                    uiOutput("taxaIn")
                                ),
                                column(
                                    4,
                                    h3(""),
                                    bsButton("cusTaxa", "Browse...")
                                )
                            )
                        ),

                        h5(""),
                        conditionalPanel(
                            condition = "input.autoUpdateSelected == false",
                            bsButton(
                                "plotCustom",
                                "Update plot",
                                style = "warning",
                                icon("refresh")
                            )
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
                            condition = "input.do > 0",
                            createProfilePlotUI("customizedProfile")
                        )
                    )
                )
            ),

            # FUNCTION TAB =====================================================
            navbarMenu(
                "Function",
                # * Profiles clustering ----------------------------------------
                tabPanel(
                    "Profiles clustering",
                    h4(strong("Profiles clustering")),
                    bsAlert("descClusteringUI"),

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
                                    value = FALSE
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
                    "Distribution analysis",
                    h4(strong("Distribution analysis")),
                    bsAlert("descDistributionUI"),

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
                    "Gene age estimation",
                    h4(strong("Gene age estimation")),
                    bsAlert("descGeneAgeUI"),

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
                                bsButton("geneAgeProtConfig", "Plot config")
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
                    "Core gene identification",
                    h4(strong("Core gene identification")),
                    bsAlert("descCoreGeneUI"),

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
                                    value = 100,
                                    step = 5
                                )
                            ),
                            column(
                                12,
                                uiOutput("taxaListCore.ui"),
                                bsButton("browseTaxaCore", "Browse")
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
                    "Group comparison",
                    h4(strong("Group comparison")),
                    bsAlert("descGCUI"),
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
                                ),
                                bsPopover(
                                    "useCommonAncestor",
                                    "",
                                    paste0(
                                        "All taxa that have the same common ",
                                        "ancestor with the selected taxa above",
                                        " will be considered as the in-group"
                                    ),
                                    "bottom"
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
                                bsButton("uploadGC", "Upload")
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
                                popify(
                                    sliderInput(
                                        "significance",
                                        paste("Significance level:"),
                                        min = 0,
                                        max = 1,
                                        step = 0.01,
                                        value = c(0.05),
                                        width = 300
                                    ),
                                    "",
                                    "P-value cut-off of the statistic test, OR
                                    cut-off of delta means between 2 groups"
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
                                popify(
                                    actionButton(
                                        "gcPlotConfig", "Plot config"
                                    ),
                                    "",
                                    "Change the appearance of the plots"
                                ),
                                hr(),
                                bsButton(
                                    "doCompare", "COMPARE!", style = "danger"
                                ),
                                h5(),
                                bsButton(
                                    "updateGC",
                                    "Update plot",
                                    style = "warning",
                                    icon("refresh")
                                )
                            )
                        )
                    ),
                    groupComparisonUI("groupComparison")
                ),
                
                # * NCBI taxonomy data -----------------------------------------
                tabPanel(
                    "NCBI taxonomy data",
                    bsAlert("descNcbiTaxDbUI"),
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
                            bsButton(
                                "doUpdateNcbi",
                                "Do update",
                                style = "warning",
                                icon("wrench")
                            ),
                            hr(),
                            verbatimTextOutput("updateNCBITaxStatus")
                        ),
                        conditionalPanel(
                            condition = "input.taxDB=='reset'",
                            h4(strong("Reset taxonomy data")),
                            bsButton(
                                "doResetTax",
                                "Do reset",
                                style = "warning",
                                icon("wrench")
                            ),
                            hr(),
                            verbatimTextOutput("resetTaxonomyDataStatus")
                        ),
                        conditionalPanel(
                            condition =
                                "input.taxDB=='export'",
                            h4(strong("Export current taxonomy files")),
                            shinyDirButton(
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
                            bsButton(
                                "doExportTax",
                                "Do export",
                                style = "warning",
                                icon("file-export")
                            ),
                            hr(),
                            verbatimTextOutput("exportTaxonomyDataStatus")
                        ),
                        conditionalPanel(
                            condition =
                                "input.taxDB=='import'",
                            h4(strong("Import your own taxonomy files")),
                            shinyDirButton(
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
                            bsButton(
                                "doImportTax",
                                "Do import",
                                style = "warning",
                                icon("file-import")
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
                
                # * Export plot settings ---------------------------------------
                tabPanel(
                    "Export plot settings",
                    h4(strong("Export plot settings")),
                    bsAlert("descExportSettingUI"),
                    radioButtons(
                        inputId = "exportSetting",
                        label = "as:",
                        choices = list(
                            "a list" = "list",
                            "an Rscript" = "rscript"
                        )
                    ),
                    hr(),
                    strong("Output dir:"),
                    br(), br(),
                    shinyDirButton(
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
                    bsButton(
                        "doExportSetting",
                        "Do export",
                        style = "warning",
                        icon("file-export")
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

        # LIST OF POP-UP WINDOWS ===============================================

        # * popup for getting taxa from OMA browser ----------------------------
        bsModal(
            "getOmaDataWindows",
            "Get OMA data",
            "openOmaWindows",
            size = "small",
            selectInput(
                "selectedOmaType",
                label = "Select type of OMA orthologs:",
                choices = list("HOG", "OG"),# "PAIR"),
                selected = "HOG"
            ),
            bsButton("getDataOma", "Get data", style = "danger"),
            downloadButton("downloadFilesOma", "Save data"),
            br(),
            em("This windows will close automatically when eveything is done!",
               style = "color:red")
        ),

        # * popup for adding new taxa from input file --------------------------
        bsModal(
            "addTaxaWindows",
            "Add new taxa",
            "addTaxa",
            size = "medium",
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
                "Rank (e.g. \"norank\" (for strain), species, order, etc.)",
                "norank",
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
            bsAlert("wrongNewTaxa"),
            hr(),
            bsButton(
                "newDone", "Finish adding", style = "warning", disabled = TRUE
            )
        ),

        # * popup for confirming parsing taxa from input file ------------------
        bsModal(
            "parseConfirm",
            "Get taxonomy info",
            "butParse",
            size = "small",
            HTML(
                '<p>Fetching Missing Taxonomy Information and
                Post-processing.</p><p><em>This windows will close
                automatically when eveything is done. Please wait...</em></p>
                <p><strong><span style="color: #ff0000;">PLEASE RELOAD THIS
                TOOL WHEN FINISHED!!!</span></strong></p>'
            )
        ),

        # * popup for plotting detailed plot -----------------------------------
        bsModal(
            "modalBs",
            "Detailed plot",
            "detailedBtn",
            size = "large",
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
            bsButton(
                "doDomainPlot", "Show domain architecture", disabled = TRUE
            ),
            uiOutput("checkDomainFiles"),
            br(),
            h4("Sequence:"),
            verbatimTextOutput("fasta"),
            br(),
            h4("Links:"),
            uiOutput("dbLink")
        ),

        # * popup for plotting domain architecture plot ------------------------
        bsModal(
            "plotArchi",
            "Domain architecture",
            "doDomainPlot",
            size = "large",
            fluidRow(
                column(
                    2, createPlotSize("archiHeight", "Plot height(px)", 400)
                ),
                column(
                    2, createPlotSize("archiWidth", "Plot width(px)", 800)
                ),
                column(
                    2,
                    createTextSize("titleArchiSize", "Title size(px)", 14, 150)
                ),
                column(
                    2,
                    createTextSize("labelArchiSize","Domain ID size(px)",12,150)
                )
            ),
            uiOutput("test.ui"),
            createArchitecturePlotUI("archiPlot")
        ),

        # * popup for setting plot colors (profiles) ---------------------------
        bsModal(
            "color",
            "Set colors for profile",
            "setColor",
            size = "small",
            colourpicker::colourInput(
                "lowColorVar1",
                "Low variable 1 (dot)",
                value = "#FF8C00"
            ),            
            colourpicker::colourInput(
                "midColorVar1",
                "Mid variable 1 (dot)",
                value = "#40ABCF"
            ),
            colourpicker::colourInput(
                "highColorVar1",
                "High variable 1 (dot)",
                value = "#164294"
            ),
            numericInput(
                "midVar1",
                "Mitpoint varriable 1",
                min = 0,
                max = 1,
                step = 0.01,
                value = 0.5
            ),
            actionButton(
                "defaultColorVar1",
                "Default",
                style = "padding:4px; font-size:100%"
            ),
            hr(),
            colourpicker::colourInput(
                "lowColorVar2",
                "Low variable 2 (background)",
                value = "#CC8D8D"
            ),
            colourpicker::colourInput(
                "midColorVar2",
                "Mid variable 2 (background)",
                value = "#FFFFFF"
            ),
            colourpicker::colourInput(
                "highColorVar2",
                "High variable 2 (background)",
                value = "#616587"
            ),
            numericInput(
                "midVar2",
                "Mitpoint varriable 2",
                min = 0,
                max = 1,
                step = 0.01,
                value = 1
            ),
            actionButton(
                "defaultColorVar2",
                "Default",
                style = "padding:4px; font-size:100%"
            ),
            hr(),
            colourpicker::colourInput(
                "paraColor",
                "Color for inparalogs",
                value = "#07d000"
            ),
            actionButton(
                "defaultColorPara",
                "Default",
                style = "padding:4px; font-size:100%"
            )
        ),

        # * popup for FASTA upload ---------------------------------------------
        bsModal(
            "fastaUploadBs",
            "FASTA upload",
            "fastaUpload",
            size = "small",
            selectInput(
                "inputType", "Choose location for:",
                c("Concatenated fasta file", "Fasta folder")
            ),
            hr(),
            uiOutput("defaultColorPara.ui"),
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
        ),

        # * popup for upload gene category -------------------------------------
        bsModal(
            "uploadGeneCategoryBs",
            "Upload gene categories",
            "uploadGeneCategory",
            size = "small",
            fileInput("geneCategory", "")
        ),

        # * popup for setting Main plot configurations -------------------------
        bsModal(
            "mainPlotConfigBs",
            "Plot appearance configuration",
            "mainPlotConfig",
            size = "small",
            column(
                6, createTextSize("xSize", "X-axis label size (px)", 14, 100)
            ),
            column(
                6, createTextSize("ySize", "Y-axis label size (px)", 14, 100)
            ),
            column(
                6,
                createTextSize("legendSize", "Legend label size (px)", 8, 150)
            ),
            column(
                6,
                selectInput(
                    "mainLegend", label = "Legend position:",
                    choices = list("Right" = "right",
                                   "Left" = "left",
                                   "Top" = "top",
                                   "Bottom" = "bottom",
                                   "Hide" = "none"),
                    selected = "right",
                    width = 150
                )
            ),
            column(
                12,
                HTML("<strong>Angle for x-axis label</strong>:<br>"),
                sliderInput(
                    "xAngle",
                    "",
                    min = 0,
                    max = 90,
                    step = 10,
                    value = 60,
                    width = 250
                ),
                br()
            ),
            column(
                12,
                HTML("<strong>Zooming factor (α) for dots on
                    profile</strong>:<br>"),
                sliderInput(
                    "dotZoom", "",
                    min = -1,
                    max = 3,
                    step = 0.1,
                    value = 0,
                    width = 250
                ),
                HTML("<em>dot size = (1+α)*defaultSize<br>defaultSize
                    =[0:5]</em>"),
                uiOutput("dotSizeInfo"),
                br()
            ),
            br(),
            hr(),
            bsButton("resetMainConfig", "Reset", style = "danger"),
            bsButton("applyMainConfig", "Done", style = "warning")
        ),

        # * popup for setting Customized plot configurations -------------------
        bsModal(
            "selectedPlotConfigBs",
            "Plot appearance configuration",
            "selectedPlotConfig",
            size = "small",
            column(
                6,
                createTextSize("xSizeSelect", "X-axis label size (px)", 14, 100)
            ),
            column(
                6,
                createTextSize("ySizeSelect", "Y-axis label size (px)", 14, 100)
            ),

            column(
                6,
                createTextSize("legendSizeSelect", "Legend label size (px)",
                               8, 150)
            ),
            column(
                6,
                selectInput(
                    "selectedLegend", label = "Legend position:",
                    choices = list("Right" = "right",
                                   "Left" = "left",
                                   "Top" = "top",
                                   "Bottom" = "bottom",
                                   "Hide" = "none"),
                    selected = "right",
                    width = 150
                )
            ),
            column(
                12,
                HTML("<strong>Angle for x-axis label</strong>:<br>"),
                sliderInput(
                    "xAngleSelect", "",
                    min = 0,
                    max = 90,
                    step = 10,
                    value = 60,
                    width = 250
                ),
                br()
            ),
            column(
                12,
                HTML("<strong>Zooming factor (α) for dots on
                     profile</strong>:<br>"),
                sliderInput(
                    "dotZoomSelect", "",
                    min = -1,
                    max = 3,
                    step = 0.1,
                    value = 0,
                    width = 250
                ),
                HTML("<em>dot size = (1+α)*defaultSize<br>
                     defaultSize=[0:5]</em>"),
                uiOutput("dotSizeInfoSelect"),
                br()
            ),
            br(),
            hr(),
            bsButton("resetSelectedConfig", "Reset", style = "danger"),
            bsButton("applySelectedConfig", "Done", style = "warning")
        ),

        # * popup for setting Gene age plot configurations ---------------------
        bsModal(
            "geneAgeProtConfigBs",
            "Plot appearance configuration",
            "geneAgeProtConfig",
            size = "small",
            sliderInput(
                "geneAgeWidth",
                "Width zoom (*600px)",
                min = 0,
                max = 5,
                step = 0.1,
                value = 1,
                width = "100%"
            ),
            sliderInput(
                "geneAgeHeight", "Height zoom (*150px)",
                min = 0,
                max = 5,
                step = 0.1,
                value = 1,
                width = "100%"
            ),
            sliderInput(
                "geneAgeText", "Text size zoom",
                min = 0,
                max = 5,
                step = 0.1,
                value = 1,
                width = "100%"
            ),
            br(),
            hr(),
            bsButton(
                "resetGeneAgeProtConfig",
                "Reset",
                style = "danger"
            )
        ),

        # * popup for setting Group compariosn plot configurations -------------
        bsModal(
            "gcPlotConfigBs",
            "Plot appearance configuration",
            "gcPlotConfig",
            size = "large",
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
            bsButton("resetConfigGC", "Reset", style = "danger"),
            bsButton("applyConfigGC", "Done", style = "warning")
        ),

        # * popup for select taxa on Customized Profile ------------------------
        bsModal(
            "cusTaxaBs",
            "Select taxon/taxa of interest",
            "cusTaxa",
            size = "small",
            selectTaxonRankUI("selectTaxonRank"),
            checkboxInput(
                "applyCusTaxa",
                strong("Apply to customized profile",
                       style = "color:red"),
                value = FALSE
            )
        ),

        # * popup for select taxa on Core gene finding -------------------------
        bsModal(
            "browseTaxaCoreBs",
            "Select taxon/taxa of interest",
            "browseTaxaCore",
            size = "small",
            selectTaxonRankUI("selectTaxonRankCore"),
            checkboxInput(
                "applyCoreTaxa",
                strong("Apply", style = "color:red"),
                value = FALSE
            )
        ),

        # * popup for select taxa on Group comparison --------------------------
        bsModal(
            "taxaGCBs",
            "Select taxon/taxa of interest",
            "taxaGC",
            size = "small",
            selectTaxonRankUI("selectTaxonRankGC"),
            checkboxInput(
                "applyTaxonGC",
                strong("Apply",
                       style = "color:red"),
                value = FALSE
            )
        ),

        # * popup for input in-group and out-group on Group comparison ---------
        bsModal(
            "uploadGCBs",
            "Upload files for group comparison",
            "uploadGC",
            size = "large",
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
                    bsButton(
                        "detailedBtn",
                        "Detailed plot",
                        style = "success",
                        disabled = FALSE
                    )
                ),
                style = "opacity: 0.80"
            )
        )
    )
)
