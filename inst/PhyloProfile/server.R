#' Import function files
sourceFiles = list.files( path = "R", pattern = "*.R$", full.names = TRUE)
lapply(sourceFiles, source, .GlobalEnv)
library(PhyloProfile)

checkpoint0 <- Sys.time()
rtCheck <- FALSE

#' set size limit for input (9999mb)
options(
    scipen = 999,
    shiny.maxRequestSize = 9999 * 1024 ^ 2, # size limit for input 9999mb
    htmlwidgets.TOJSON_ARGS = list(na = 'string')
)

#' MAIN SERVER =================================================================
shinyServer(function(input, output, session) {
    homePath = c(wd='~/')
    # Automatically stop a Shiny app when closing the browser tab
    session$allowReconnect(TRUE)

    nameFullFile <- paste0(
        getwd(), "/data/preProcessedTaxonomy.txt"
    )
    currentNCBIinfo <- NULL
    if (file.exists(nameFullFile))
        currentNCBIinfo <- as.data.frame(data.table::fread(nameFullFile))

    fastModeCutoff <- 600
    # =========================== INITIAL CHECKING  ============================
    # * check for internet connection ------------------------------------------
    observe({
        if (hasInternet() == FALSE) shinyjs::toggleState("demoData")
    })

    output$noInternetMsg <- renderUI({
        if (hasInternet() == FALSE && input$demoData != "preCalcDt") {
            strong(
                em("Internet connection is required for using demo data!"),
                style = "color:red"
            )
        } else return()
    })

    # * check for the existence of taxonomy files ------------------------------
    observe({
        fileExist <- file.exists("data/preProcessedTaxonomy.txt")
        if (fileExist == FALSE) {
            msg <- paste0(
                "Please wait while preprocessed data are being downloaded!!!"
            )
            shinyBS::createAlert(
                session, "fileExistMsgUI", "fileExistMsg", title = "",
                content = msg,
                append = FALSE
            )
        } else shinyBS::closeAlert(session, "fileExistMsg")
    })

    observe({
        if (!file.exists(isolate("data/rankList.txt"))) {
            withProgress(message = '1/6 rankList.txt...', value = 0.5, {
                data(rankList)
                write.table(
                    rankList, file = "data/rankList.txt",
                    col.names = FALSE,
                    row.names = FALSE,
                    quote = FALSE,
                    sep = "\t"
                )
            })
        }
    })

    observe({
        if (!file.exists(isolate("data/idList.txt"))) {
            withProgress(message = '2/6 idList.txt...', value = 0.5, {
                data(idList)
                write.table(
                    idList, file = "data/idList.txt",
                    col.names = FALSE,
                    row.names = FALSE,
                    quote = FALSE,
                    sep = "\t"
                )
            })
        }
    })

    observe({
        if (!file.exists(isolate("data/taxonNamesReduced.txt"))) {
            withProgress(message = '3/6 taxonNamesReduced.txt...', value = 0.5,{
                data(taxonNamesReduced)
                write.table(
                    taxonNamesReduced, file = "data/taxonNamesReduced.txt",
                    col.names = TRUE,
                    row.names = FALSE,
                    quote = FALSE,
                    sep = "\t"
                )
            })
        }
    })

    observe({
        if (!file.exists(isolate("data/taxonomyMatrix.txt"))) {
            withProgress(message = '4/6 taxonomyMatrix.txt...', value = 0.5, {
                data(taxonomyMatrix)
                write.table(
                    taxonomyMatrix, file = "data/taxonomyMatrix.txt",
                    col.names = TRUE,
                    row.names = FALSE,
                    quote = FALSE,
                    sep = "\t"
                )
            })
        }
    })

    observe({
        if (!file.exists(isolate("data/preProcessedTaxonomy.txt"))) {
            withProgress(
                message = '5/6 preProcessedTaxonomy.txt...', value = 0.5, {
                    if (hasInternet() == TRUE) {
                        preProcessedTaxonomy <- processNcbiTaxonomy()
                        write.table(
                            preProcessedTaxonomy,
                            file = "data/preProcessedTaxonomy.txt",
                            col.names = TRUE,
                            row.names = FALSE,
                            quote = FALSE,
                            sep = "\t"
                        )
                        currentNCBIinfo <- preProcessedTaxonomy
                    } else {
                        system(
                            "cp data/newTaxa.txt data/preProcessedTaxonomy.txt"
                        )
                    }
                }
            )
            shinyBS::closeAlert(session, "fileExistMsg")
        }
    })

    observe({
        if (!file.exists(isolate("data/preCalcTree.nw"))) {
            withProgress(
                message = '6/6 preCalcTree.nw...', value = 0.5, {
                    ppTaxonomyMatrix <- getTaxonomyMatrix()
                    preCalcTree <- createUnrootedTree(ppTaxonomyMatrix)
                    ape::write.tree(preCalcTree, file = "data/preCalcTree.nw")
                }
            )
            shinyBS::closeAlert(session, "fileExistMsg")
        }
    })

    # ========================== LOAD CONFIG FILE ==============================
    configFile <- .GlobalEnv$configFile
    i_mainInput <- i_domainInput <- i_fastaInput <- NULL
    i_treeInput <- i_sortedTaxaInput <- NULL
    i_rank <- i_refspec <- NULL
    i_cluster <- FALSE
    i_profileType <- i_distMethod <- i_clusterMethod <- NULL
    i_xAxis <- NULL
    i_colorByGroup <- i_orderGenes <- i_geneCategory <- i_geneName <- NULL
    if (!is.null(configFile)) {
        if (file.exists(configFile)) {
            configs <- yaml::read_yaml(configFile)
            i_mainInput <- configs$mainInput
            i_domainInput <- configs$domainInput
            i_treeInput <- configs$treeInput
            i_sortedTaxaInput <- configs$sortedTaxaInput
            i_fastaInput <- configs$fastaInput
            i_rank <- configs$rank
            i_refspec <- configs$refspec
            i_cluster <- configs$clusterProfile
            i_profileType <- configs$profileTypeClustering
            i_distMethod <- configs$distMethodClustering
            i_clusterMethod <- configs$clusterMethod
            i_xAxis <- configs$xAxis
            i_colorByGroup <- configs$colorByGroup
            i_orderGenes <- configs$orderGenes
            i_geneCategory <- configs$geneCategory
            i_geneName <- configs$geneName
        } else configFile <- NULL

        if (!file.exists(i_mainInput))
            stop(paste("Main input file", i_mainInput ,"not found!"))
    }
    if (!is.logical(i_cluster)) i_cluster <- FALSE

    # * render main input ------------------------------------------------------
    observe({
        if (!is.null(i_mainInput)) {
            updateSelectInput(
                session, "demoData", label = "",
                choices = list("Pre-calculated" = "preCalcDt"),
                selected = "preCalcDt"
            )
        }
    })
    # * change order taxa by user input tree or sorted taxon list---------------
    observe({
        if (!is.null(i_treeInput)) {
            updateRadioButtons(
                session,
                inputId = "orderTaxa",
                label = "",
                choices = list(
                    "automatically", "by user defined tree/sorted list",
                    "by a sorted list"
                ),
                selected = "by user defined tree/sorted list",
                inline = TRUE
            )
        }
    })
    observe({
        if (!is.null(i_sortedTaxaInput)) {
            updateRadioButtons(
                session,
                inputId = "orderTaxa",
                label = "",
                choices = list(
                    "automatically", "by user defined tree/sorted list",
                    "by a sorted list"
                ),
                selected = "by a sorted list",
                inline = TRUE
            )
        }
    })
    # * change x-axis ----------------------------------------------------------
    observe({
        if (!is.null(i_xAxis)) {
            if (i_xAxis == "taxa" || i_xAxis == "genes") {
                updateRadioButtons(
                    inputId = "xAxis",
                    label = "Choose type of x-axis:",
                    choices = list("taxa", "genes"),
                    selected = i_xAxis,
                    inline = TRUE
                )
            }
        }
    })
    # * render colorByGroup option ---------------------------------------------
    observe({
        if (!is.null(i_colorByGroup)) {
            updateCheckboxInput(
                session, "colorByGroup", value = as.logical(i_colorByGroup)
            )
        }
    })
    # * render orderGenes option -----------------------------------------------
    observe({
        if (!is.null(i_orderGenes)) {
            updateCheckboxInput(
                session, "orderGenes", value = i_orderGenes
            )
        }
    })
    # * prepare gene categories ------------------------------------------------
    if (!is.null(i_geneCategory) && !file.exists(i_geneCategory)){
        stop(paste("Gene categories file ", i_geneCategory ,"not found!"))
    }
    # * prepare gene names -----------------------------------------------------
    if (!is.null(i_geneName) && !file.exists(i_geneName)){
        stop(paste("Gene names file ", i_geneName ,"not found!"))
    }
    # * apply clustering profiles ----------------------------------------------
    if (is.null(i_profileType) ||
        (i_profileType != "binary" && i_profileType != "var1" &&
         i_profileType != "var2")
    ) i_profileType <- "binary"

    if (i_profileType != "binary") {
        if (is.null(i_distMethod) ||
            (i_distMethod != "mutualInformation" &&
             i_distMethod != "distanceCorrelation")
        ) i_distMethod <- "mutualInformation"
    } else {
        if (is.null(i_distMethod) ||
            (i_distMethod != "euclidean" && i_distMethod != "maximum" &&
             i_distMethod != "manhattan" && i_distMethod != "canberra" &&
             i_distMethod != "binary" && i_distMethod != "pearson" &&
             i_distMethod != "mutualInformation" &&
             i_distMethod != "distanceCorrelation")
        ) i_distMethod <- "euclidean"
    }

    observe({
        if (i_cluster) {
            updateCheckboxInput(
                session,
                "applyCluster",
                "Apply clustering to profile plot",
                value = TRUE
            )

            if (i_clusterMethod) {
                updateSelectInput(
                    session,
                    "clusterMethod",
                    label = "Cluster method:",
                    choices = list(
                        "single" = "single",
                        "complete" = "complete",
                        "average (UPGMA)" = "average",
                        "mcquitty (WPGMA)" = "mcquitty",
                        "median (WPGMC)" = "median",
                        "centroid (UPGMC)" = "centroid"
                    ),
                    selected = i_clusterMethod
                )
            } else {
                updateSelectInput(
                    session,
                    "clusterMethod",
                    label = "Cluster method:",
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
            }
        }
    })

    # ======================== INPUT & SETTINGS TAB ============================
    # * Render input message ---------------------------------------------------
    observe({
        filein <- input$mainInput
        if (is.null(filein) & input$demoData == "none") {
            msg <- paste0(
                "PhyloProfile is ready to use! Please upload
        <em>an input file</em>, <em>a folder of pre-processed data</em>, or
        select <em>a demo data</em><br /> to begin!
        To learn more about the <em>input data</em>, please visit
        <span style=\"text-decoration: underline;\">
        <a href=\"https://github.com/BIONF/PhyloProfile/wiki/Input-Data\">
        <span style=\"color: #ff0000;\">our wiki</span></span></a>."
            )
            shinyBS::createAlert(
                session, "inputMsgUI", "inputMsg", title = "",
                content = msg,
                append = FALSE
            )
        } else shinyBS::closeAlert(session, "inputMsg")
    })

    # * check the validity of input file and render inputCheck.ui --------------
    output$inputCheck.ui <- renderUI({
        filein <- input$mainInput
        if (is.null(filein)) return()
        inputType <- inputType()

        if (inputType[1] == "noGeneID") {
            shinyBS::updateButton(session, "do", disabled = TRUE)
            HTML(
                "<font color=\"red\"><em><strong>ERROR: Unsupported input
                format.<a
                href=\"https://github.com/BIONF/PhyloProfile/wiki/Input-Data\"
                target=\"_blank\">Click here for more
                info</a></em></strong></font>"
            )
        } else if (inputType[1] == "emptyCell") {
            shinyBS::updateButton(session, "do", disabled = TRUE)
            em(strong("ERROR: Rows have unequal length",style= "color:red"))
        } else if (inputType[1] == "moreCol") {
            shinyBS::updateButton(session, "do", disabled = TRUE)
            em(strong(
                "ERROR: More columns than column names", style = "color:red"
            ))
        } else if (inputType[1] == "invalidFormat") {
            shinyBS::updateButton(session, "do", disabled = TRUE)
            em(strong(
                "ERROR: Invalid format", style = "color:red"
            ))
        } else {
            validType = c("xml", "fasta", "wide", "long", "oma")
            if (!(inputType[1] %in% validType)) {
                shinyBS::updateButton(session, "do", disabled = TRUE)
                invalidOma <- paste(inputType, collapse = "; ")
                msg <- paste0("ERROR: Invalid IDs found! ", invalidOma)
                em(strong(msg, style = "color:red"))
            } else {
                shinyBS::updateButton(session, "do", disabled = FALSE)
                return()
            }
        }
    })

    # * render download link for Demo online files -----------------------------
    output$downloadDemo.ui <- renderUI({
        if (input$demoData == "ampk-tor" || input$demoData == "arthropoda") {
            em(a(
                "Click here to download demo files",
                href = "https://github.com/BIONF/phyloprofile-data",
                target = "_blank"
            ))
        }
    })

    output$mainInputFile.ui <- renderUI({
        if (input$demoData == "arthropoda") {
            url <- paste0(
                "https://raw.githubusercontent.com/BIONF/",
                "phyloprofile-data/master/arthropoda.zip"
            )
            strong(a("Download demo data", href = url, target = "_blank"))
        } else if (input$demoData == "ampk-tor") {
            url <- paste0(
                "https://raw.githubusercontent.com/BIONF/phyloprofile-data/",
                "master/ampk-tor.zip"
            )
            strong(a("Download demo data", href = url, target = "_blank"))
        } else if (input$demoData == "preCalcDt") {
            em("Input from config file")
        } else {
            list(
                radioButtons(
                    "mainInputType", "Upload input:",
                    choices = c("file","folder"), selected = "file",
                    inline = TRUE
                ),
                conditionalPanel(
                    condition = "input.mainInputType == 'file'",
                    fileInput("mainInput", "")
                ),
                conditionalPanel(
                    condition = "input.mainInputType == 'folder'",
                    shinyFiles::shinyDirButton(
                        "mainInputDir", "Browse folder" ,
                        title = "Please select a folder",
                        buttonType = "default", class = NULL
                    ),
                    br(), br(),
                    verbatimTextOutput("mainInputDirInfo")
                )
            )
        }
    })

    output$domainInputFile.ui <- renderUI({
        if (input$demoData == "arthropoda") {
            strong("Download demo data (link above)")
        } else if (input$demoData == "ampk-tor") {
            strong("Download demo data (link above)")
        } else if (input$demoData == "preCalcDt") {
            if(!is.null(i_domainInput)) {
                em("Input from config file")
            } else {
                if (input$annoLocation == "from file") {
                    fileInput("fileDomainInput", "")
                } else textInput(
                    "domainPath", "", "",
                    placeholder = "Give full path to domain directory"
                )
            }
        } else {
            if (input$annoLocation == "from file") {
                fileInput("fileDomainInput", "")
            } else textInput(
                "domainPath", "", "",
                placeholder = "Give full path to domain directory"
            )
        }
    })

    # * render description for Demo data ---------------------------------------
    output$demoDataDescribe <- renderUI({
        if (input$demoData == "none") return()
        else if (input$demoData == "ampk-tor") {
            url <- paste0(
                "https://github.com/BIONF/phyloprofile-data/blob/master/",
                "ampk-tor.md"
            )
            em(a("Data description", href = url, target = "_blank"))
        } else if (input$demoData == "arthropoda"){
            url <- paste0(
                "https://github.com/BIONF/phyloprofile-data/blob/master/",
                "arthropoda.md"
            )
            em(a("Data description", href = url, target = "_blank"))
        }
    })

    # * check main input folder ------------------------------------------------
    getMainInputDir <- reactive({
        shinyFiles::shinyDirChoose(
            input, "mainInputDir", roots = homePath, session = session
        )
        outputPath <- shinyFiles::parseDirPath(homePath, input$mainInputDir)
        return(replaceHomeCharacter(as.character(outputPath)))
    })

    checkMainInputDir <- reactive({
        mainDir <- getMainInputDir()
        rdsFiles <- list.files(
            path = mainDir, pattern = "\\.rds$", ignore.case = TRUE
        )
        reqFiles <- c(
            "fullData.rds", "longDf.rds", "preData.rds", "sortedtaxaList.rds"
        )
        availFile <- intersect(reqFiles, rdsFiles)
        return(paste0(mainDir, "/", availFile))
    })

    output$mainInputDirInfo <- renderText({
        req(getMainInputDir())
        if (length(checkMainInputDir()) < 4) {
            paste(
                "Required RDS file(s) not found!",
                "Please upload a single input file!"
            )
        } else {
            shinyBS::closeAlert(session, "inputMsg")
            getMainInputDir()
        }
    })

    # * check OMA input --------------------------------------------------------
    output$checkOmaInput <- reactive({
        filein <- input$mainInput
        if (is.null(filein)) return()
        inputType() == "oma"
    })
    outputOptions(output, "checkOmaInput", suspendWhenHidden = FALSE)

    # * download OMA data after parsing ----------------------------------------
    output$downloadFilesOma <- downloadHandler(
        filenname <- function() {
            "omaDataToPhyloprofileInput.zip"
        },
        content <- function(file) {
            write.table(
                getMainInput(), "phyloprofile.txt",
                sep = "\t",
                row.names = FALSE,
                col.names = TRUE,
                quote = FALSE
            )

            write.table(
                getAllFastaOma(finalOmaDf()), "fasta.txt",
                sep = "\t",
                row.names = FALSE,
                col.names = FALSE,
                quote = FALSE
            )

            write.table(
                getAllDomainsOma(finalOmaDf()), "domain.txt",
                sep = "\t",
                row.names = FALSE,
                col.names = FALSE,
                quote = FALSE
            )

            zip(
                zipfile = file,
                files = c("phyloprofile.txt", "domain.txt", "fasta.txt")
            )
        },
        contentType = "application/zip"
    )

    # * close OMA parsing popup windows -------------------------------------
    observeEvent(input$getDataOma, {
        toggleModal(session, "getOmaDataWindows", toggle = "close")
        shinyBS::updateButton(session, "getDataOma", disabled = TRUE)
        shinyjs::toggleState("mainInput")
        shinyjs::toggleState("fileDomainInput")
        shinyjs::toggleState("fastaUpload")
    })

    # * render textinput for Variable 1 & 2 ------------------------------------
    output$var1ID.ui <- renderUI({
        longDataframe <- getMainInput()
        if (is.null(longDataframe)) {
            textInput(
                "var1ID",
                h5("1st variable:"),
                value = "Variable 1",
                width = "100%",
                placeholder = "Name of first variable"
            )
        } else {
            textInput(
                "var1ID", h5("1st variable:"),
                value = colnames(longDataframe)[4],
                width = "100%",
                placeholder = "Name of first variable"
            )
        }
    })

    output$var2ID.ui <- renderUI({
        longDataframe <- getMainInput()
        if (is.null(longDataframe)) {
            textInput(
                "var2ID",
                h5("2st variable:"),
                value = "Variable 2",
                width = "100%",
                placeholder = "Name of second variable"
            )
        } else {
            textInput(
                "var2ID", h5("2st variable:"),
                value = colnames(longDataframe)[5],
                width = "100%",
                placeholder = "Name of second variable"
            )
        }
    })

    # * update colorVar for heatmapPlottingFast --------------------------------
    output$colorVar.ui <- renderUI({
        choices <- setNames( c("var1","var2"), c(input$var1ID, input$var2ID))
        if (input$plotMode == "fast")
            radioButtons( #selectInput(
                "colorVar", "Color dots by:", choices = choices,
                selected = "var1", inline = TRUE
            )
    })

    # * get total number of genes ----------------------------------------------
    output$totalGeneNumber.ui <- renderUI({
        geneList <- getMainInput()
        out <- as.list(levels(factor(geneList$geneID)))

        listIn <- input$geneList
        if (!is.null(listIn)) {
            geneListDf <- read.table(file = listIn$datapath, header = FALSE)
            out <- as.list(unique(geneListDf$V1))
        }
        if (length(out) > 0) {
            strong(paste0("Total number of genes:  ", length(out)))
        }
    })

    # * check the existance of the input concatenate fasta file ----------------
    output$concatFasta.existCheck <- renderUI({
        req(input$concatFasta)
        f <- input$concatFasta$datapath
        if (!file.exists(f)) {
            helpText("File not exists!!")
        } else {
            if (length(readLines(f, n = 1)) == 0) {
                helpText("is not a fasta file!!")
            } else {
                firstLine <- readLines(f, n = 1)
                a <- substr(firstLine, 1, 1)
                if (a == ">") {
                    HTML(
                        '<p><span style="color: #0000ff;">
                        <strong>Please click CLOSE to comfirm!
                        </strong></span></p>'
                    )
                } else helpText("is not a fasta file!!")
                }
            }
    })

    # * render location of taxonomy DB -----------------------------------------
    getUserTaxDBpath <- reactive({
        shinyFiles::shinyDirChoose(
            input, "taxDbDir", roots = homePath, session = session
        )
        outputPath <- shinyFiles::parseDirPath(homePath, input$taxDbDir)
        return(replaceHomeCharacter(as.character(outputPath)))
    })

    checkUserTaxDB <- reactive({
        req(getUserTaxDBpath())
        missingFiles <- c()
        if (!(file.exists(paste0(getUserTaxDBpath(),"/idList.txt"))))
            missingFiles <- c(missingFiles, "idList.txt")
        if (!(file.exists(paste0(getUserTaxDBpath(),"/rankList.txt"))))
            missingFiles <- c(missingFiles, "rankList.txt")
        if (!(file.exists(paste0(getUserTaxDBpath(),"/taxonNamesReduced.txt"))))
            missingFiles <- c(missingFiles, "taxonNamesReduced.txt")
        if (!(file.exists(paste0(getUserTaxDBpath(),"/taxonomyMatrix.txt"))))
            missingFiles <- c(missingFiles, "taxonomyMatrix.txt")
        if (length(missingFiles) == 0) {
            if (!(file.exists(paste0(getUserTaxDBpath(),"/preCalcTree.nw")))) {
                withProgress(
                    message = 'Create preCalcTree.nw...', value = 0.5, {
                        ppTaxonomyMatrix <-getTaxonomyMatrix(getUserTaxDBpath())
                        preCalcTree <- createUnrootedTree(ppTaxonomyMatrix)
                        ape::write.tree(
                            preCalcTree,
                            file = paste0(getUserTaxDBpath(),"/preCalcTree.nw")
                        )
                    }
                )
            }
        }
        return(missingFiles)
    })

    output$userTaxDBwarning <- renderUI({
        req(getUserTaxDBpath())
        if (length(checkUserTaxDB() > 0)) {
            msg <- paste0(
                "These files are missing: ",
                paste(checkUserTaxDB(), collapse = "; "),
                ". Default DB will be used!"
            )
            em(msg)
        }
    })

    getTaxDBpath <- reactive({
        defaultTaxDB <- system.file(
            "PhyloProfile", "data", package = "PhyloProfile", mustWork = TRUE
        )
        if (!file.exists(paste0(defaultTaxDB, "/taxonomyMatrix.txt"))) {
            if (file.exists(paste0(getwd(), "/data/taxonomyMatrix.txt"))) {
                return(paste0(getwd(), "/data"))
            }
        }
        if (length(getUserTaxDBpath()) > 0) {
            if (length(checkUserTaxDB() > 0))
                return(defaultTaxDB)
            else
                return(getUserTaxDBpath())
        } else {
            return(defaultTaxDB)
        }
    })

    output$taxDbPath <- renderText({
        return(getTaxDBpath())
    })

    # * render upload input for sorted gene list -------------------------------
    output$inputSortedGenes.ui <- renderUI({
        req(input$mainInputType)
        if (
            !(((input$mainInputType == "file" & is.null(input$mainInput)) |
                (input$mainInputType == "folder" & is.null(getMainInputDir()))
            ) & input$demoData == "none")
        ) {
            fileInput("inputSortedGenes", "")
        }
    })

    checkInputSortedGenes <- reactive({
        filein <- input$inputSortedGenes
        if (!(is.null(filein))) {
            sortedGeneInputDf <- read.table(
                file = filein$datapath,
                header = FALSE,
                check.names = FALSE,
                comment.char = "",
                fill = FALSE
            )
            sortedGenes <- unique(sortedGeneInputDf$V1)
            allGenes <- as.character(unique(getMainInput()[,c("geneID")]))
            if (!(all(sortedGeneInputDf$V1 %in% allGenes))) {
                if (length(setdiff(allGenes, sortedGenes)) > 0) {
                    return(
                        list(missing = setdiff(allGenes, sortedGenes))
                    )
                } else if (length(setdiff(sortedGenes, allGenes)) > 0){
                    return(
                        list(more = setdiff(sortedGenes, allGenes))
                    )
                } else {
                    # actually this condition will never happen
                    return(
                        list(
                            diff = union(
                                setdiff(allGenes, sortedGenes),
                                setdiff(sortedGenes, allGenes)
                            )
                        )
                    )
                }
            } else return(list(sortedGenes = sortedGenes))
        }
    })

    output$checkSortedGenes.ui <- renderUI({
        checkStatus <- checkInputSortedGenes()
        if (length(checkStatus) == 0) return(strong(em("Please upload file!")))
        if (names(checkStatus[1]) == "missing") {
            msg <- paste0(
                "<p><em><span style=\"color: #ff0000;\">ERROR: ",
                "Some genes missings from the list! ",
                "(", checkStatus[1],")",
                "</span></em></p>"
            )
        } else if (names(checkStatus[1]) == "more"){
            msg <- paste0(
                "<p><em><span style=\"color: #ff0000;\">WARNING: ",
                "Some genes from the list not present in the input! ",
                "(", checkStatus[1],")",
                "</span></em></p>"
            )
        } else if (names(checkStatus[1]) == "diff") {
            msg <- paste0(
                "<p><em><span style=\"color: #ff0000;\">ERROR: ",
                "Gene IDs from the list are not identical with the ",
                "IDs in the phylogenetic profile input! ",
                "(", checkStatus[1],")",
                "</span></em></p>"
            )
        } else msg <- NULL
        HTML(msg)
    })

    # * check the validity of input tree file and render checkNewick.ui --------
    output$inputTree.ui <- renderUI({
        filein <- input$mainInput
        if (!(is.null(filein) & input$demoData == "none")) {
            fileInput("inputTree", "")
        }
    })

    checkNewickID <- reactive({
        tree <- NULL
        if (input$demoData == "preCalcDt") {
            if (!is.null(i_treeInput)) {
                tree <- read.table(
                    file = i_treeInput,
                    header = FALSE,
                    check.names = FALSE,
                    comment.char = "",
                    fill = FALSE
                )
            }
        } else {
            req(input$mainInput)
            req(input$inputTree)
            filein <- input$inputTree
            tree <- read.table(
                file = filein$datapath,
                header = FALSE,
                check.names = FALSE,
                comment.char = "",
                fill = FALSE
            )
        }
        if (!is.null(tree)) {
            checkNewick <- checkNewick(tree, inputTaxonID())
            if (checkNewick == 0) shinyBS::updateButton(session, "do", disabled = FALSE)
            return(checkNewick)
        } else return()
    })

    output$checkNewick.ui <- renderUI({
        checkNewick <- checkNewickID()
        req(checkNewick)
        if (checkNewick == 1) {
            shinyBS::updateButton(session, "do", disabled = TRUE)
            HTML("<p><em><span style=\"color: #ff0000;\"><strong>
            ERROR: Parenthesis(-es) missing!</strong></span></em></p>")
        } else if (checkNewick == 2) {
            shinyBS::updateButton(session, "do", disabled = TRUE)
            HTML("<p><em><span style=\"color: #ff0000;\"><strong>
            ERROR: Comma(s) missing!</strong></span></em></p>")
        } else if (checkNewick == 3) {
            shinyBS::updateButton(session, "do", disabled = TRUE)
            HTML("<p><em><span style=\"color: #ff0000;\"><strong>
            ERROR: Tree contains singleton!</strong></span></em></p>")
        } else if (checkNewick == 0) {
            return()
        } else {
            shinyBS::updateButton(session, "do", disabled = TRUE)
            strong(
                em(paste0(checkNewick, " not exist in main input file!")),
                style = "color:red"
            )
        }
    })

    # * render upload input for sorted taxon list ------------------------------
    output$inputSortedTaxa.ui <- renderUI({
        filein <- input$mainInput
        if (!(is.null(filein) & input$demoData == "none")) {
            fileInput("inputSortedTaxa", "")
        }
    })

    output$checkSortedTaxa.ui <- renderUI({
        msg <- paste0(
            "<p><em>Check <a href=\"https://github.com/BIONF/PhyloProfile/",
            "wiki/Input-Data#list-of-sorted-taxa\">this link</a> for more ",
            "details!</em></p>"
        )
        sortedTaxonInputDf <- NULL
        if (input$demoData == "preCalcDt") {
            if (!is.null(i_sortedTaxaInput)) {
                sortedTaxonInputDf <- read.table(
                    file = i_sortedTaxaInput,
                    header = FALSE,
                    check.names = FALSE,
                    comment.char = "",
                    fill = FALSE
                )
            }
        } else {
            req(input$mainInput)
            filein <- input$inputSortedTaxa
            if (!is.null(filein)) {
                sortedTaxonInputDf <- read.table(
                    file = filein$datapath,
                    header = FALSE,
                    check.names = FALSE,
                    comment.char = "",
                    fill = FALSE
                )
            }
        }
        sortedTaxonList <- sortedTaxonInputDf$V1
        if (length(sortedTaxonList) > 0) {
            if ("FALSE" %in% grepl("ncbi", sortedTaxonList)) {
                shinyBS::updateButton(session, "do", disabled = TRUE)
                msg <- paste0("<p><em><span style=\"color: #ff0000;\"><strong>
            ERROR: Invalid input!</strong></span></em></p>")
            } else {
                msg <- ""
            }
        }
        HTML(msg)
    })

    observeEvent(input$inputSortedTaxa,  ({
        if (!is.null(input$inputSortedTaxa)) {
            shinyjs::toggleState("rankSelect")
            shinyjs::toggleState("inSelect")
        }
    }))

    # * reset profile plot colors ----------------------------------------------
    observeEvent(input$defaultColorVar2, {
        shinyjs::reset("lowColorVar2")
        shinyjs::reset("midColorVar2")
        shinyjs::reset("highColorVar2")
        shinyjs::reset("midVar2")
    })

    observeEvent(input$defaultColorVar1, {
        shinyjs::reset("lowColorVar1")
        shinyjs::reset("midColorVar1")
        shinyjs::reset("highColorVar1")
        shinyjs::reset("midVar1")
    })

    observeEvent(input$defaultColorPara, {
        shinyjs::reset("paraColor")
    })

    # * render list of taxonomy ranks ------------------------------------------
    output$rankSelect <- renderUI({
        if (input$demoData == "arthropoda") {
            selectInput(
                "rankSelect", label = "",
                choices = getTaxonomyRanks(),
                selected = "class"
            )
        } else if (input$demoData == "ampk-tor") {
            selectInput(
                "rankSelect", label = "",
                choices = getTaxonomyRanks(),
                selected = "species"
            )
        } else if (input$demoData == "preCalcDt") {
            selectedRank <- "species"
            if (!is.null(i_rank)) {
                selectedRank <- i_rank
            }
            selectInput(
                "rankSelect", label = "",
                choices = getTaxonomyRanks(),
                selected = selectedRank
            )
        } else {
            req(input$mainInputType)
            if (input$mainInputType == "folder") {
                if (length(checkMainInputDir()) == 4) {
                    tmpDf <- readRDS(
                        paste0(getMainInputDir(),"/sortedtaxaList.rds")
                    )
                    selectedRank <- unique(tmpDf$rank)
                    if (length(selectedRank) > 1) selectedRank <- "strain"
                    selectInput(
                        "rankSelect", label = "",
                        choices = selectedRank, #getTaxonomyRanks(),
                        selected = selectedRank
                    )
                }
            } else {
                longDataframe <- getMainInput()
                req(longDataframe)
                lowestRank <- getLowestRank(longDataframe, getTaxDBpath())
                if (lowestRank %in% getTaxonomyRanks()) {
                    selectInput(
                        "rankSelect", label = "",
                        choices = getTaxonomyRanks(),
                        selected = lowestRank
                    )
                } else {
                    selectInput(
                        "rankSelect", label = "",
                        choices = getTaxonomyRanks(),
                        selected = "species"
                    )
                }
            }
        }
    })

    # * predict reference species ----------------------------------------------
    refSpec <- reactive({
        req(getMainInput())
        if (rtCheck) {
            checkpoint121 <- Sys.time()
            print(paste("checkpoint121 - before predict refspec",checkpoint121))
        }
        tmp <- as.data.table(
            getMainInput() %>% dplyr::select(geneID, ncbiID) %>% distinct()
        )
        
        refspecID <- tmp[, .N, by = ncbiID][order(-N)][1, ncbiID]
        refsecpName <- getInputTaxaName(
            input$rankSelect, refspecID, getTaxDBpath()
        )
        if (rtCheck) {
            checkpoint122 <- Sys.time()
            print(
                paste("checkpoint122 - predict refspec done", checkpoint122, 
                      " --- ",  checkpoint122 - checkpoint121)
            )
        }
        return(refsecpName)
    })

    # * render list of (super)taxa ---------------------------------------------
    observe({
        choice <- inputTaxonName()
        choice$fullName <- as.factor(choice$fullName)
        if (input$demoData == "arthropoda") {
            hellemDf <- data.frame(
                "name" = c(
                    "Drosophila melanogaster",
                    "Drosophila melanogaster",
                    "Drosophila",
                    "Drosophilidae",
                    "Diptera",
                    "Insecta",
                    "Arthropoda",
                    "Metazoa",
                    "Eukaryota"
                ),
                "rank" = c(
                    "strain",
                    "species",
                    "genus",
                    "family",
                    "order",
                    "class",
                    "phylum",
                    "kingdom",
                    "superkingdom"
                )
            )
            rankName <- input$rankSelect
            updateSelectizeInput(
                session, "inSelect", "", server = TRUE,
                choices = as.list(levels(choice$fullName)),
                selected = hellemDf$name[hellemDf$rank == rankName]
            )
        } else if (input$demoData == "ampk-tor") {
            humanDf <- data.frame(
                "name" = c(
                    "Homo sapiens",
                    "Homo sapiens",
                    "Homo",
                    "Hominidae",
                    "Primates",
                    "Mammalia",
                    "Chordata",
                    "Metazoa",
                    "Eukaryota"
                ),
                "rank" = c(
                    "strain",
                    "species",
                    "genus",
                    "family",
                    "order",
                    "class",
                    "phylum",
                    "kingdom",
                    "superkingdom"
                )
            )
            rankName <- input$rankSelect
            updateSelectizeInput(
                session, "inSelect", "", server = TRUE,
                choices = as.list(levels(choice$fullName)),
                selected = humanDf$name[humanDf$rank == rankName]
            )
        } else if (input$demoData == "preCalcDt") {
            refspec <- levels(choice$fullName)[1]
            if (!is.null(i_refspec)) {
                if (i_refspec %in% levels(choice$fullName))
                    refspec <- i_refspec
            }
            updateSelectizeInput(
                session, "inSelect", "", server = TRUE,
                choices = as.list(levels(choice$fullName)),
                selected = refspec
            )
        } else {
            if (input$mainInputType == "folder") {
                req(getMainInputDir())
                sortedTaxa <- readRDS(
                    paste0(getMainInputDir(),"/sortedtaxaList.rds")
                )
                refSpec <- sortedTaxa %>%
                    filter(supertaxon == levels(sortedTaxa$supertaxon)[1]) %>%
                    select(supertaxonID)
                updateSelectizeInput(
                    session, "inSelect", "", server = TRUE,
                    choices = as.list(levels(choice$fullName)),
                    selected = unique(
                        choice$fullName[choice$ncbiID == unique(refSpec$supertaxonID)]
                    )
                )
            } else {
                predRefspec <- refSpec()
                selectedRefspec <- levels(choice$fullName)[1]
                if (nrow(predRefspec) > 0) {
                    if (predRefspec$fullName %in% levels(choice$fullName))
                        selectedRefspec <- predRefspec$fullName
                }
                if (length(choice$fullName) > 0) {
                    # remove cache first
                    updateSelectizeInput(
                        session, "inSelect", choices = NULL, server = TRUE
                    )
                    # then update the list
                    updateSelectizeInput(
                        session, "inSelect", "", server = TRUE,
                        choices = as.list(levels(choice$fullName)),
                        selected = selectedRefspec
                    )
                }
            }
        }
    })

    # * enable "PLOT" button ---------------------------------------------------
    observeEvent(input$rankSelect,  ({
        if (input$rankSelect == "") {
            shinyBS::updateButton(session, "do", disabled = TRUE)
        } else {
            unkTaxa <- unkTaxa()
            if (length(unkTaxa) == 0) {
                shinyBS::updateButton(session, "do", disabled = FALSE)
            }
        }
    }))
    # * move to main tab when "PLOT" button has been clicked -------------------
    observe({
        # use tabsetPanel "id" argument to change tabs
        if (input$do > 0) {
            updateTabsetPanel(session, "tabs", selected = "Main profile")
        }
    })

    # * disable main input, genelist input and demo data checkbox --------------
    observe({
        if (input$do > 0) {
            shinyjs::disable("mainInputType")
            shinyjs::disable("mainInput")
            shinyjs::disable("mainInputDir")
            shinyjs::disable("geneListSelected")
            shinyjs::disable("demoData")
            shinyjs::disable("rankSelect")
            shinyjs::disable("taxDbLoc")
        }
    })

    observe({
        req(input$mainInputType)
        if (input$mainInputType == "folder") {
            shinyjs::disable("inSelect")
            shinyjs::disable("orderGenes")
            shinyjs::disable("orderTaxa")
            shinyjs::disable("geneListSelected")
        } else {
            shinyjs::enable("inSelect")
            shinyjs::enable("rankSelect")
            shinyjs::enable("orderGenes")
            shinyjs::enable("orderTaxa")
            shinyjs::enable("geneListSelected")
        }
    })

    # * enable/disable update plot button --------------------------------------
    observe({
        if (input$autoUpdate == TRUE) {
            shinyBS::updateButton(session, "applyFilter", disabled = TRUE)
            shinyBS::updateButton(session, "updateBtn", disabled = TRUE)
        } else {
            shinyBS::updateButton(session, "applyFilter", disabled = FALSE)
            shinyBS::updateButton(session, "updateBtn", disabled = FALSE)
        }
    })

    # =========================== RENDER FILTER SLIDEBARS ======================

    # * render filter slidebars for Main plot ----------------------------------
    output$var1Cutoff.ui <- renderUI({
        createSliderCutoff(
            "var1", paste(input$var1ID, "cutoff:"), 0.0, 1.0, input$var1ID
        )
    })

    output$var2Cutoff.ui <- renderUI({
        createSliderCutoff(
            "var2", paste(input$var2ID, "cutoff:"), 0.0, 1.0, input$var2ID
        )
    })

    output$percentCutoff.ui <- renderUI({
        createSliderCutoff(
            "percent", "% of present taxa:", 0.0, 1.0, "percent"
        )
    })

    # * render filter slidebars for Customized plot ----------------------------
    output$var1Filter.ui <- renderUI({
        req(input$var1)
        createSliderCutoff(
            "var1cus",
            paste(input$var1ID, "cutoff:"),
            input$var1[1], input$var1[2], input$var1ID
        )
    })

    output$var2Filter.ui <- renderUI({
        req(input$var2)
        createSliderCutoff(
            "var2cus",
            paste(input$var2ID, "cutoff:"),
            input$var2[1], input$var2[2], input$var2ID
        )
    })

    output$percentFilter.ui <- renderUI({
        req(input$percent)
        createSliderCutoff(
            "percent2",
            "% of present taxa:",
            input$percent[1], input$percent[2], "percent"
        )
    })

    output$coorthologFilter.ui <- renderUI({
        numericInput(
            "coortholog2",
            "Max co-orthologs",
            min = 1,
            max = 999,
            step = 1,
            value = input$coortholog,
            width = 150
        )
    })

    # * render filter slidebars for Distribution plot --------------------------
    output$var1Dist.ui <- renderUI({
        createSliderCutoff(
            "var1Dist",
            paste(input$var1ID, "cutoff:"),
            input$var1[1], input$var1[2], input$var1ID
        )
    })

    output$var2Dist.ui <- renderUI({
        createSliderCutoff(
            "var2Dist",
            paste(input$var2ID, "cutoff:"),
            input$var2[1], input$var2[2], input$var2ID
        )
    })

    output$percentDist.ui <- renderUI({
        createSliderCutoff(
            "percentDist",
            "% of present taxa:",
            input$percent[1], input$percent[2], "percent"
        )
    })

    # * render filter slidebars for Gene age estimation plot -------------------
    output$var1Age.ui <- renderUI({
        createSliderCutoff(
            "var1Age",
            paste(input$var1ID, "cutoff:"),
            input$var1[1], input$var1[2], input$var1ID
        )
    })

    output$var2Age.ui <- renderUI({
        createSliderCutoff(
            "var2Age",
            paste(input$var2ID, "cutoff:"),
            input$var2[1], input$var2[2], input$var2ID
        )
    })

    output$percentAge.ui <- renderUI({
        createSliderCutoff(
            "percentAge",
            "% of present taxa:",
            input$percent[1], input$percent[2], "percent"
        )
    })

    # * render filter slidebars for Core gene finding function -----------------
    output$var1Core.ui <- renderUI({
        createSliderCutoff(
            "var1Core", paste(input$var1ID, "cutoff:"), 0.0, 1.0,
            input$var1ID
        )
    })

    output$var2Core.ui <- renderUI({
        createSliderCutoff(
            "var2Core", paste(input$var2ID, "cutoff:"), 0.0, 1.0,
            input$var2ID
        )
    })

    output$percentCore.ui <- renderUI({
        createSliderCutoff(
            "percentCore",
            "% of present taxa:",
            0, 1, "percent"
        )
    })

    # * update value for filter slidebars of Main Plot -------------------------
    # ** based on customized profile
    observe({
        newVar1 <- input$var1cus
        updateSliderCutoff(
            session,
            "var1", paste(input$var1ID, "cutoff:"), newVar1, input$var1ID
        )
    })

    observe({
        newVar2 <- input$var2cus
        updateSliderCutoff(
            session,
            "var2", paste(input$var2ID, "cutoff:"), newVar2, input$var2ID
        )
    })

    observe({
        newPercent <- input$percent2
        updateSliderCutoff(
            session,
            "percent", "% of present taxa:", newPercent, "percent"
        )
    })

    observe({
        newCoortholog <- input$coortholog2
        updateNumericInput(
            session,
            "coortholog",
            value = newCoortholog
        )
    })

    # ** based on "Distribution analysis"
    observe({
        newVar1 <- input$var1Dist
        updateSliderCutoff(
            session,
            "var1", paste(input$var1ID, "cutoff:"), newVar1, input$var1ID
        )
    })

    observe({
        newVar2 <- input$var2Dist
        updateSliderCutoff(
            session,
            "var2", paste(input$var2ID, "cutoff:"), newVar2, input$var2ID
        )
    })

    observe({
        newPercent <- input$percentDist
        updateSliderCutoff(
            session,
            "percent", "% of present taxa:", newPercent, "percent"
        )
    })

    # ** based on "Gene age estimation"
    observe({
        newVar1 <- input$var1Age
        updateSliderCutoff(
            session,
            "var1", paste(input$var1ID, "cutoff:"), newVar1, input$var1ID
        )
    })

    observe({
        newVar2 <- input$var2Age
        updateSliderCutoff(
            session,
            "var2", paste(input$var2ID, "cutoff:"), newVar2, input$var2ID
        )
    })

    observe({
        newPercent <- input$percentAge
        updateSliderCutoff(
            session,
            "percent", "% of present taxa:", newPercent, "percent"
        )
    })

    # * reset cutoffs of Main plot ---------------------------------------------
    observeEvent(input$resetMain, {
        shinyjs::reset("var1")
        shinyjs::reset("var2")
        shinyjs::reset("percent")
        shinyjs::reset("coortholog")
    })

    # * reset cutoffs of Customized plot ---------------------------------------
    observeEvent(input$resetSelected, {
        shinyjs::reset("var1")
        shinyjs::reset("var2")
        shinyjs::reset("percent")
        shinyjs::reset("coortholog")
    })

    # ========================= PARSING UNKNOWN TAXA ===========================
    # * get list of "unknown" taxa in main input -------------------------------
    unkTaxa <- reactive({
        withProgress(message = 'Checking for unknown taxa...', value = 0.5, {
            longDataframe <- getMainInput()
            req(longDataframe)

            inputTaxa <- levels(longDataframe$ncbiID)
            inputTaxa <- unlist(strsplit(inputTaxa, split = "\t"))

            if (inputTaxa[1] == "geneID") {
                # remove "geneID" element from vector inputTaxa
                inputTaxa <- inputTaxa[-1]
            }

            taxDB <- getTaxDBpath()
            if (!file.exists(isolate(paste0(taxDB, "/rankList.txt")))) {
                return(inputTaxa)
            } else {
                info <- file.info(paste0(taxDB, "/rankList.txt"))
                if (info$size == 0) {
                    return(inputTaxa)
                } else {
                    rankListFile <- paste0(taxDB, "/rankList.txt")
                    allTaxa <- as.factor(
                        unlist(fread(file = rankListFile, select = 1))
                    )

                    # list of unknown taxa
                    unkTaxa <- inputTaxa[!(inputTaxa %in% allTaxa)]
                    if (identical(unkTaxa, character(0))) return()

                    # get non-ncbi taxa
                    unkTaxa <- data.frame(TaxonID = unkTaxa)
                    unkTaxa$id <- as.numeric(substring(unkTaxa$TaxonID, 5))
                    unkTaxa$Source <- "ncbi"

                    nameFullFile <- paste0(
                        getwd(), "/data/preProcessedTaxonomy.txt"
                    )
                    ncbiTaxa <- as.factor(
                        unlist(fread(file = nameFullFile, select = 1))
                    )

                    ncbiID <- levels(ncbiTaxa)
                    maxNCBI <- max(sort(as.numeric(ncbiID[ncbiID != "ncbiID"])))

                    unkTaxaId <- c()
                    if (nrow(unkTaxa[!(unkTaxa[,"id"] %in% ncbiTaxa),]) > 0) {
                        unkTaxaId <- unkTaxa[!(unkTaxa$id %in% ncbiTaxa),]$id
                        unkTaxa[unkTaxa$id %in% unkTaxaId,]$Source <- "unknown"
                    }

                    newTaxaFile <- paste0(taxDB, "/newTaxa.txt")
                    newTaxaDf <- fread(file = newTaxaFile, select = 1)
                    newTaxa <- as.numeric(newTaxaDf$ncbiID)

                    if (nrow(unkTaxa[unkTaxa$id %in% newTaxa,]) > 0) {
                        unkTaxa[unkTaxa$id %in% newTaxa,]$Source <- "new"
                    }

                    # check invalid tax IDs in newTaxa.txt file
                    if (any(as.numeric(newTaxa) < as.numeric(maxNCBI))) {
                        invalidNewtaxa <- length(
                            newTaxa[as.numeric(newTaxa) < as.numeric(maxNCBI)]
                        )
                        unkTaxa[nrow(unkTaxa) + 1,] <-
                            c(
                                paste(invalidNewtaxa, "IDs in newTaxa.txt"),
                                as.character(invalidNewtaxa),"invalid"
                            )
                    }

                    # check for invalid taxon IDs
                    if (any(as.numeric(unkTaxaId) < as.numeric(maxNCBI))) {
                        unkTaxa[
                            unkTaxa$id %in% unkTaxaId
                            & as.numeric(unkTaxa$id) < as.numeric(maxNCBI),
                        ]$Source <- "invalid"
                    }

                    # return list of unkTaxa
                    return(unkTaxa)
                }
            }
        })
    })

    # * check the status of unkTaxa --------------------------------------------
    output$unkTaxaStatus <- reactive({
        unkTaxa <- unkTaxa()
        if (length(unkTaxa) > 0) {
            if ("invalid" %in% unkTaxa$Source) return("invalid")
            if ("unknown" %in% unkTaxa$Source) return("unknown")
            else return("ncbi")
        } else return(0)
    })
    outputOptions(output, "unkTaxaStatus", suspendWhenHidden = FALSE)

    # * render list of unkTaxa -------------------------------------------------
    output$unkTaxaFull <-
        DT::renderDataTable(options = list(searching = FALSE, pageLength = 10),{
            if (length(unkTaxa()) > 0) {
                tb <- unkTaxa()
                tb[, c("TaxonID", "Source")]
            }
        })

    # * download list of unkTaxa -----------------------------------------------
    output$unkTaxa.download <- downloadHandler(
        filename = function() {
            c("unknownTaxa.txt")
        },
        content = function(file) {
            dataOut <- unkTaxa()
            dataOut <- dataOut[, c("TaxonID", "Source")]
            write.table(
                dataOut, file, sep = "\t", row.names = FALSE, quote = FALSE
            )
        }
    )

    # * update the form for adding new taxa ------------------------------------
    newTaxa <- reactiveValues()
    newTaxa$Df <- data.frame(
        "ncbiID" = numeric(),
        "fullName" = character(),
        "rank" = character(),
        "parentID" = numeric(),
        stringsAsFactors = FALSE
    )
    newIndex <- reactiveValues()
    newIndex$value <- 1

    observeEvent(input$newAdd, {
        newTaxa$Df[newIndex$value, ] <- c(
            input$newID, input$newName, input$newRank, input$newParent
        )
        newIndex$value <- newIndex$value + 1
        updateTextInput(session, "newID", value = as.numeric(input$newID) + 1)
        updateTextInput(session, "newName", value = "")
        updateTextInput(session, "newRank", value = "norank")
        updateTextInput(session, "newParent", value = "")
        shinyjs::enable("newDone")
    })

    # * get info for new taxa from uploaded file -------------------------------
    newTaxaFromFile <- reactive({
        filein <- input$newTaxaFile
        req(filein)

        tmpDf <- read.table(
            file = filein$datapath,
            sep = "\t",
            header = TRUE,
            check.names = FALSE,
            comment.char = ""
        )
        if (ncol(tmpDf) != 4) {
            shinyBS::createAlert(
                session, "wrongNewTaxa", "wrongNewTaxaMsg",
                content = "Wrong format. Please check your file!",
                append = FALSE
            )
            shinyjs::disable("newDone")
            return()
        } else {
            shinyBS::createAlert(
                session, "wrongNewTaxa", "wrongNewTaxaMsg",
                content = "Click Finish adding to continue!",
                append = FALSE
            )
            shinyjs::enable("newDone")
            colnames(tmpDf) <- c("ncbiID", "fullName", "rank", "parentID")
            newTaxa$Df <- tmpDf
            return(newTaxa$Df)
        }
    })

    observeEvent(input$newTaxaFile, {
        newTaxaFromFile()
    })

    # * close adding taxa windows ----------------------------------------------
    observeEvent(input$newDone, {
        toggleModal(session, "addTaxaWindows", toggle = "close")
        write.table(
            newTaxa$Df, "data/newTaxa.txt",
            sep = "\t",
            eol = "\n",
            row.names = FALSE,
            quote = FALSE
        )
    })

    # * check if data is loaded and "parse" button is clicked and confirmed ----
    v1 <- reactiveValues(parse = FALSE)
    observeEvent(input$butParse, {
        toggleModal(session, "parseConfirm", toggle = "close")
        v1$parse <- input$butParse
        shinyBS::updateButton(session, "butParse", disabled = TRUE)
        shinyjs::toggleState("newTaxaAsk")
        shinyjs::toggleState("mainInput")
    })

    # * create rankList, idList, taxonNamesReduced, taxonomyMatrix -------------
    # * and preCalcTree files
    invalidID <- reactive({
        req(input$mainInputType)
        if (input$mainInputType == "folder") {
            return()
        }

        filein <- input$mainInput
        req(filein)
        inputType <- inputType()

        if (inputType == "xml" |
            inputType == "long" |
            inputType == "wide" |
            inputType == "fasta" |
            inputType == "oma") {

            if (v1$parse == FALSE) return()
            else {
                taxDB <- getTaxDBpath()
                inputDf <- read.table(
                    file = filein$datapath,
                    sep = "\t",
                    header = TRUE,
                    check.names = FALSE,
                    comment.char = ""
                )

                # get list of taxa need to be parsed (taxa mising taxonomy info)
                if (v1$parse == TRUE) {
                    unkTaxaDf <- unkTaxa()
                    unkTaxa <- as.character(substring(unkTaxaDf$TaxonID, 5))
                }

                invalidID <- data.frame(
                    "id" = as.character(),
                    "type" = as.character(),
                    stringsAsFactors = FALSE
                )

                ## join all ncbi taxa and new taxa together
                ncbiTaxonInfo <- fread("data/preProcessedTaxonomy.txt")
                newTaxaFromFile <- fread(
                    paste0(taxDB, "/newTaxa.txt"),
                    colClasses = c("ncbiID" = "character")
                )
                allTaxonInfo <- as.data.frame(
                    rbindlist(list(newTaxaFromFile, ncbiTaxonInfo))
                )

                ## check missing ids
                if (any(!(unkTaxa %in% allTaxonInfo$ncbiID))) {
                    invalidMissing <-
                        unkTaxa[!(unkTaxa %in% allTaxonInfo$ncbiID)]
                    invalidIDTmp <- data.frame(
                        "id" = invalidMissing,
                        "type" = rep("missing", length(invalidMissing))
                    )
                    invalidID <- rbindlist(list(invalidID, invalidIDTmp))
                }

                ## check IDs & names from newTaxa that are present in
                ## taxonNamesFull
                if (nrow(newTaxaFromFile[newTaxaFromFile$ncbiID
                                         %in% ncbiTaxonInfo$ncbiID,]) > 0) {
                    invalidID <- newTaxaFromFile[
                        newTaxaFromFile$ncbiID %in% ncbiTaxonInfo$ncbiID,
                        ]$ncbiID
                    invalidIDTmp <- data.frame(
                        "id" = invalidID,
                        "type" = rep("id", length(invalidID))
                    )
                    invalidID <- rbindlist(list(invalidID, invalidIDTmp))

                    newTaxaFromFile <- newTaxaFromFile[
                        !(newTaxaFromFile$ncbiID %in% ncbiTaxonInfo$ncbiID),]
                }

                if (nrow(newTaxaFromFile[newTaxaFromFile$fullName
                                         %in% ncbiTaxonInfo$fullName,]) > 0) {
                    invalidName <- newTaxaFromFile[
                        newTaxaFromFile$fullName %in% ncbiTaxonInfo$fullName,
                        ]$ncbiID
                    invalidIDTmp <- data.frame(
                        "id" = invalidName,
                        "type" = rep("name", length(invalidName))
                    )
                    invalidID <- rbindlist(list(invalidID, invalidIDTmp))
                }

                if (nrow(invalidID) > 0) return(invalidID)

                ## parse taxonomy info
                withProgress(
                    message = "Parsing new taxa...", value = 0, {
                        taxonomyInfo <- getIDsRank(unkTaxa, allTaxonInfo)
                        rankList <- as.data.frame(taxonomyInfo[2])
                        idList <- as.data.frame(taxonomyInfo[1])
                        reducedInfoList <- as.data.frame(taxonomyInfo[3])
                    }
                )

                withProgress(
                    message = "Generating taxonomy file...",
                    value = 0, {
                        # open existing files
                        # (idList, rankList and taxonNamesReduced.txt)
                        ncol <- max(
                            count.fields(paste0(taxDB,"/rankList.txt"),sep="\t")
                        )
                        oldIDList <- fread(
                            paste0(taxDB, "/idList.txt"),
                            sep = "\t",
                            header = FALSE,
                            check.names = FALSE,
                            fill = TRUE,
                            stringsAsFactors = TRUE,
                            na.strings = c("", "NA"),
                            col.names = paste0("X", seq_len(ncol))
                        )
                        oldRankList <- fread(
                            paste0(taxDB, "/rankList.txt"),
                            sep = "\t",
                            header = FALSE,
                            check.names = FALSE,
                            fill = TRUE,
                            stringsAsFactors = TRUE,
                            na.strings = c("", "NA"),
                            col.names = paste0("X", seq_len(ncol))
                        )
                        oldNameList <- fread(
                            paste0(taxDB, "/taxonNamesReduced.txt"),
                            sep = "\t",
                            header = TRUE,
                            check.names = FALSE,
                            fill = TRUE,
                            stringsAsFactors = TRUE
                        )

                        # and append new info into those files
                        newIDList <- rbindlist(
                            list(oldIDList, idList), fill = TRUE
                        )
                        newRankList <- rbindlist(
                            list(oldRankList, rankList), fill = TRUE
                        )
                        newNameList <- rbindlist(
                            list(reducedInfoList, oldNameList), fill = TRUE
                        )
                        newNameList <- newNameList[
                            !duplicated(newNameList$ncbiID),
                        ]

                        # write output files
                        # (idList, rankList and taxonNamesReduced)
                        write.table(
                            newIDList[!duplicated(newIDList), ],
                            file  = paste0(taxDB, "/idList.txt"),
                            col.names = FALSE,
                            row.names = FALSE,
                            quote = FALSE,
                            sep = "\t"
                        )
                        write.table(
                            newRankList[!duplicated(newRankList), ],
                            file = paste0(taxDB, "/rankList.txt"),
                            col.names = FALSE,
                            row.names = FALSE,
                            quote = FALSE,
                            sep = "\t"
                        )
                        write.table(
                            newNameList[!duplicated(newNameList), ],
                            file = paste0(taxDB, "/taxonNamesReduced.txt"),
                            col.names = TRUE,
                            row.names = FALSE,
                            quote = FALSE,
                            sep = "\t"
                        )

                        # create taxonomy matrix (taxonomyMatrix.txt)
                        taxMatrix <- taxonomyTableCreator(
                            paste0(taxDB, "/idList.txt"),
                            paste0(taxDB, "/rankList.txt")
                        )
                        write.table(
                            taxMatrix,
                            file = paste0(taxDB, "/taxonomyMatrix.txt"),
                            sep = "\t",
                            eol = "\n",
                            row.names = FALSE,
                            quote = FALSE
                        )

                        # prepare matrix for calculating distances
                        preCalcTree <- createUnrootedTree(taxMatrix)
                        ape::write.tree(
                            preCalcTree, file = paste0(taxDB, "/preCalcTree.nw")
                        )
                    }
                )
            }
        }
        return()
    })

    # * output invalid NCBI ID -------------------------------------------------
    output$invalidID.output <- renderTable({
        req(invalidID())
        outDf <- invalidID()
        colnames(outDf) <- c("Invalid ID(s)", "Type")
        return(outDf)
    })

    # * download list of invalidID ---------------------------------------------
    output$invalidID.download <- downloadHandler(
        filename = function() {
            c("invalidIDs.txt")
        },
        content = function(file) {
            dataOut <- invalidID()
            colnames(dataOut) <- c("Invalid ID(s)", "Type")
            write.table(
                dataOut, file, sep = "\t", row.names = FALSE, quote = FALSE
            )
        }
    )

    # * render final msg after taxon parsing -----------------------------------
    output$endParsingMsg <- renderUI({
        if (is.null(invalidID())) {
            strong(
                h4("PLEASE RELOAD THIS TOOL WHEN FINISHED!!!"),
                style = "color:red"
            )
        } else {
            HTML(
                '<p><strong><span style="color: #e12525;"> SOME INVALID TAXON
                 IDs HAVE BEEN FOUND!!</span><br /> </strong></p>
                <p><em>Type="<span style="color: #0000ff;">id</span>"/
                <span style="color: #0000ff;">name</span>:
                IDs/names already exist in NCBI!</em></p>
                <p><em>Type="<span style="color: #0000ff;">missing</span>": IDs
                cannot be found in both NCBI and newTaxa.txt file.</em></p>
                <p>For IDs with type of <em><span style="color: #0000ff;">"id"
                </span></em> and <em><span style="color: #0000ff;">"name"</span>
                </em>, please remove them from newTaxa.txt file or
                renamed their IDs and names.</p>
                <p>For IDs with type of <em><span style="color: #0000ff;">
                "missing"</span></em>, please check the validity of them&nbsp;in
                <a href="https://www.ncbi.nlm.nih.gov/taxonomy" target="_blank"
                rel="noopener"> NCBI taxonomy database</a>!</p>'
            )
        }
    })

    # ====================== PROCESSING INPUT DATA =============================
    # * check if data is loaded and "plot" button is clicked -------------------
    v <- reactiveValues(doPlot = FALSE)
    observeEvent(input$do, {
        # 0 will be coerced to FALSE
        # 1+ will be coerced to TRUE
        v$doPlot <- input$do
        filein <- input$mainInput
        if (
            input$mainInputType == "file" & is.null(filein) &
            input$demoData == "none"
        ) {
            v$doPlot <- FALSE
            shinyBS::updateButton(session, "do", disabled = TRUE)
        }
    })
    w <- reactiveValues(doCusPlot = FALSE)
    observeEvent(input$plotCustom, {
        w$doCusPlot <- input$plotCustom
        filein <- input$mainInput
        if (
            input$mainInputType == "file" & is.null(filein) &
            input$demoData == "none"
        ) {
            w$doCusPlot <- FALSE
            shinyBS::updateButton(session, "plotCustom", disabled = TRUE)
        }
    })

    # * check if genes ordered by distances has been selected ------------------
    output$applyClusterCheck.ui <- renderUI({
        if (input$orderGenes != "profile similarity") {
            HTML('<p><em>(Choose "Order seed IDs by profile similarity" in
                 <strong>Input & settings tab</strong>&nbsp;to enable this
                 function)</em></p>')
        }
    })

    # * enable clustering ------------------------------------------------------
    observe({
        if (input$orderGenes != "profile similarity") {
            updateCheckboxInput(
                session, "applyCluster", value = FALSE
            )
            shinyjs::disable("applyCluster")
        } else shinyjs::enable("applyCluster")
    })

    # * get OMA data for input list --------------------------------------------
    getOmaBrowser <- function(idList, orthoType) {
        withProgress(
            message = "Retrieving OMA data",
            value = 0, {
                omaDf <- pbapply::pblapply(
                    idList,
                    function (x) getDataForOneOma(x, orthoType)
                )
            }
        )
        return(data.frame(rbindlist(omaDf, use.names = TRUE)))
    }

    finalOmaDf <- reactive({
        filein <- input$mainInput
        req(filein)
        if (inputType() == "oma") {
            if (input$getDataOma[1] == 0) return()
            omaIDs <- fread(
                file = filein$datapath,
                header = FALSE,
                stringsAsFactors = FALSE,
                select = 1
            )
            return(getOmaBrowser(omaIDs$V1, input$selectedOmaType))
        } else return()
    })

    # * get gene names (if provided) -------------------------------------------
    getGeneNames <- reactive({
        geneNameFile <- input$geneName
        if (!is.null(geneNameFile)) {
            inputNameDt <- read.table(
                file = geneNameFile$datapath,
                sep = "\t",
                header = FALSE,
                check.names = FALSE,
                comment.char = "",
                fill = TRUE
            )
            colnames(inputNameDt) <- c("geneID","geneName")
        } else if (!is.null(i_geneName)){
            inputNameDt <- read.table(
                file = i_geneName,
                sep = "\t",
                header = FALSE,
                check.names = FALSE,
                comment.char = "",
                fill = TRUE
            )
            colnames(inputNameDt) <- c("geneID","geneName")
        } else inputNameDt <- NULL
        return(inputNameDt)
    })
    
    # * check main input format (wide, long, xml, fasta, or invalidFormat) -----
    inputType <- reactive ({
        req(input$mainInputType)
        filein <- input$mainInput
        if (is.null(filein)) return()
        return(checkInputValidity(filein$datapath))
    })

    # * convert main input file in any format into long format dataframe -------
    getMainInput <- reactive({
        req(input$mainInputType)
        if (rtCheck) {
            checkpoint11 <- Sys.time()
            print(paste("checkpoint11 - before reading main input",checkpoint11))
        }
        
        withProgress(message = 'Reading main input...', value = 0.5, {
            if (input$mainInputType == "folder") {
                req(getMainInputDir())
                if (length(checkMainInputDir()) == 4) {
                    longDataframe <- readRDS(
                        paste0(getMainInputDir(),"/longDf.rds")
                    )
                } else return()
            } else {
                if (input$demoData == "arthropoda") {
                    longDataframe <- myData[["EH2547"]]
                } else if (input$demoData == "ampk-tor") {
                    longDataframe <- myData[["EH2544"]]
                } else if (input$demoData == "preCalcDt") {
                    longDataframe <- createLongMatrix(i_mainInput)
                } else {
                    filein <- input$mainInput
                    if (is.null(filein)) return()
                    if (inputType() == "oma") {
                        if (input$getDataOma[1] == 0) return()
                        longDataframe <- createProfileFromOma(finalOmaDf())
                        longDataframe <- as.data.frame(unclass(longDataframe))
                    } else {
                        longDataframe <- createLongMatrix(filein$datapath)
                    }
                }
                
                if (input$demoData == "arthropoda" | input$demoData == "ampk-tor") {
                    # convert geneID, ncbiID and orthoID into factor and
                    # var1, var2 into numeric
                    for (i in seq_len(3)) {
                        longDataframe[, i] <- as.factor(longDataframe[, i])
                    }
                    if (ncol(longDataframe) > 3 & ncol(longDataframe) < 6) {
                        for (j in seq(4, ncol(longDataframe))){
                            longDataframe[,j] <- suppressWarnings(
                                as.numeric(as.character(longDataframe[,j]))
                            )
                        }
                    }
                }
                
                # remove fdogMA orthologs if required
                if (input$showAllTaxa == FALSE)
                    longDataframe <- longDataframe[longDataframe$orthoID != "fdogMA",]
            }
            # update number of genes to plot based on input
            updateNumericInput(
                session, "endIndex",
                value = nlevels(as.factor(longDataframe$geneID))
            )
            # add geneName column if not yet exist
            if (!("geneName" %in% colnames(longDataframe))) {
                if (!(is.null(getGeneNames()))) {
                    # add gene names if specified by uploaded file
                    longDataframe <- longDataframe %>%
                        left_join(getGeneNames(), by = "geneID")
                    longDataframe$geneName[is.na(longDataframe$geneName)] <-
                        longDataframe$geneID[is.na(longDataframe$geneName)]
                } else {
                    longDataframe$geneName <- longDataframe$geneID
                }
            }
            # return
            if (rtCheck) {
                checkpoint12 <- Sys.time()
                print(paste(
                    "checkpoint12 - reading main input done at", checkpoint12, 
                    " --- ",  checkpoint12 - checkpoint11
                ))
            }
            
            return(longDataframe)
        })
    })

    # * parse domain info into data frame --------------------------------------
    getDomainInformation <- reactive({
        withProgress(message = 'Reading domain input...', value = 0.5, {
            if (v$doPlot == FALSE) return()
            if (input$demoData == "none") {
                filein <- input$mainInput
                inputType <- inputType()
            } else inputType <- "demo"

            if (inputType == "oma") {
                domainDf <- getAllDomainsOma(finalOmaDf())
            } else {
                mainInput <- getMainInput()

                if (inputType == "demo") {
                    if (input$demoData == "preCalcDt") {
                        if (!is.null(i_domainInput)) {
                            if (!file.exists(i_domainInput)) {
                                paste(
                                    "Domain file", i_domainInput ,"not found!"
                                )
                                return()
                            }
                            domainDf <- parseDomainInput(
                                NULL,
                                i_domainInput,
                                "file"
                            )
                        } else {
                            if (input$annoLocation == "from file") {
                                inputDomain <- input$fileDomainInput
                                domainDf <- parseDomainInput(
                                    NULL,
                                    inputDomain$datapath,
                                    "file"
                                )
                            } else {
                                # GET INFO BASED ON CURRENT TAB
                                if (input$tabs == "Main profile") {
                                    # info contains groupID, orthoID,
                                    # supertaxon, mVar1, %spec, var2
                                    info <- mainpointInfo()
                                } else if (input$tabs == "Customized profile") {
                                    info <- selectedpointInfo()
                                }
                                domainDf <- parseDomainInput(
                                    info[1],
                                    input$domainPath,
                                    "folder"
                                )
                            }
                        }
                    } else {
                        if (input$demoData == "arthropoda") {
                            domainDf <- myData[["EH2549"]]
                        } else {
                            domainDf <- myData[["EH2546"]]
                        }

                        domainDf$seedID <- as.character(domainDf$seedID)
                        domainDf$orthoID <- as.character(domainDf$orthoID)
                        domainDf$seedID <- gsub("\\|",":",domainDf$seedID)
                        domainDf$orthoID <- gsub("\\|",":",domainDf$orthoID)
                    }
                } else {
                    if (input$annoLocation == "from file") {
                        inputDomain <- input$fileDomainInput
                        domainDf <- parseDomainInput(
                            NULL,
                            inputDomain$datapath,
                            "file"
                        )
                    } else {
                        # GET INFO BASED ON CURRENT TAB
                        if (input$tabs == "Main profile") {
                            # info = groupID,orthoID,supertaxon,mVar1,%spec,var2
                            info <- mainpointInfo()
                        } else if (input$tabs == "Customized profile") {
                            info <- selectedpointInfo()
                        }
                        domainDf <- parseDomainInput(
                            info[1],
                            input$domainPath,
                            "folder"
                        )
                    }
                }
            }
            return(domainDf)
        })
    })

    # * get ID list of input taxa from main input ------------------------------
    inputTaxonID <- reactive({
        if (rtCheck) st <- Sys.time()
        if (input$demoData == "arthropoda" |
            input$demoData == "ampk-tor" |
            length(unkTaxa()) == 0) {
            withProgress(message = 'Getting input taxon IDs...', value = 0.5, {
                longDataframe <- getMainInput()
                inputTaxa <- getInputTaxaID(longDataframe)
                if (rtCheck) {
                    checkpoint21 <- Sys.time()
                    print(paste(
                        "checkpoint21 - get input tax IDs done", checkpoint21, 
                        " --- ",  checkpoint21 - st
                    ))
                }
                return(inputTaxa)
            })
        } else return()
    })

    # * get NAME list of all (super)taxa ---------------------------------------
    inputTaxonName <- reactive({
        req(input$rankSelect)
        if (rtCheck) st <- Sys.time()
        if (is.null(getMainInput()) & input$demoData == "none") return()
        if (length(unkTaxa()) > 0) return()
        if (input$rankSelect == "") return()
        withProgress(message = 'Getting input taxon names...', value = 0.5, {
            inputTaxaName <- getInputTaxaName(
                input$rankSelect, inputTaxonID(), getTaxDBpath()
            )
            if (rtCheck) {
                checkpoint31 <- Sys.time()
                print(paste(
                    "checkpoint31 - get input tax names done", checkpoint31, 
                    " --- ",  checkpoint31 - st
                ))    
            }
            return(inputTaxaName)
        })
    })

    # * sort taxonomy data of input taxa ---------------------------------------
    sortedtaxaList <- reactive({
        req(isTruthy(v$doPlot)|isTruthy(w$doCusPlot))
        req(input$rankSelect)
        req(input$inSelect)
        
        if (rtCheck) {
            checkpoint41 <- Sys.time()
            print(paste("checkpoint41 - before sorting taxa", checkpoint41))
        }
        
        withProgress(message = 'Sorting input taxa...', value = 0.5, {
            if (input$mainInputType == "folder") {
                req(getMainInputDir())
                if (length(checkMainInputDir()) == 4) {
                    sortedOut <- readRDS(
                        paste0(getMainInputDir(),"/sortedtaxaList.rds")
                    )
                } else return()
            } else {
                # get input taxonomy tree
                inputTaxaTree <- NULL
                if (input$demoData == "preCalcDt") {
                    if (!is.null(i_treeInput)) {
                        inputTaxaTree <- ape::read.tree(file = i_treeInput)
                    }
                } else {
                    treeIn <- input$inputTree
                    if (!is.null(treeIn)) {
                        inputTaxaTree <- ape::read.tree(file = treeIn$datapath)
                    } else {
                        preCalcTree <- paste0(getTaxDBpath(), "/preCalcTree.nw")
                        if (file.exists(preCalcTree)) 
                            inputTaxaTree <- ape::read.tree(file = preCalcTree)
                    }
                }
                # get list of sorted taxa
                sortedTaxaFile <- input$inputSortedTaxa
                if (!is.null(sortedTaxaFile)) {
                    sortedTaxonInputDf <- read.table(
                        file = sortedTaxaFile$datapath,
                        header = FALSE,
                        check.names = FALSE,
                        comment.char = "",
                        fill = FALSE
                    )
                    sortedTaxonList <- sortedTaxonInputDf$V1
                } else sortedTaxonList <- NULL
                
                # sort taxonomy matrix based on selected refTaxon
                sortedOut <- sortInputTaxa(
                    taxonIDs = inputTaxonID(),
                    rankName = input$rankSelect,
                    refTaxon = input$inSelect,
                    taxaTree = inputTaxaTree,
                    sortedTaxonList = sortedTaxonList,
                    taxDB = getTaxDBpath()
                )
            }
            # return
            if (rtCheck) {
                checkpoint42 <- Sys.time()
                print(paste(
                    "checkpoint42 - sorting taxa done", checkpoint42, " --- ",  
                    checkpoint42 - checkpoint41
                ))
            }
            
            return(sortedOut)
        })
    })

    # * get list of all input (super)taxa and their ncbi IDs -------------------
    allInputTaxa <- reactive({
        req(isTruthy(v$doPlot)|isTruthy(w$doCusPlot))
        {
            input$applyFilterCustom
            input$applyFilter
        }
        if (rtCheck) st <- Sys.time()
        allTaxa <- sortedtaxaList()
        allTaxa <- allTaxa[,c("supertaxonID", "supertaxon")]
        allTaxa$supertaxon <- factor(
            substr(
                as.character(allTaxa$supertaxon), 8 ,
                nchar(as.character(allTaxa$supertaxon))),
            levels = substr(
                levels(as.factor(allTaxa$supertaxon)), 8,
                nchar(levels(as.factor(allTaxa$supertaxon)))))
        if (rtCheck) {
            checkpoint51 <- Sys.time()
            print(paste(
                "checkpoint51 - get all input (super)taxa done", checkpoint51, 
                " --- ",  checkpoint51 - st
            ))
        }
        if (input$showAllTaxa) {
            return(allTaxa[!duplicated(allTaxa),])
        } else return()
    })

    # * count taxa for each supertaxon -----------------------------------------
    getCountTaxa <- reactive({
        req(sortedtaxaList())
        if (rtCheck) st <- Sys.time()
        taxaCount <- sortedtaxaList() %>% dplyr::group_by(supertaxon) %>%
            dplyr::summarise(n = n(), .groups = "drop")
        if (rtCheck) {
            checkpoint61 <- Sys.time()
            print(paste(
                "checkpoint61 - count taxa done", checkpoint61, " --- ",
                checkpoint61 - st
            ))
        }
        return(taxaCount)
    })

    # * get subset data for plotting -------------------------------------------
    preData <- reactive({
        req(isTruthy(v$doPlot)|isTruthy(w$doCusPlot))
        req(input$mainInputType)
        if (rtCheck) {
            checkpoint71 <- Sys.time()
            print(paste("checkpoint71 - before subseting input", checkpoint71))
        }
        ### if input a folder
        if (input$mainInputType == "folder") {
            req(getMainInputDir())
            if (length(checkMainInputDir()) == 4) {
                return(readRDS(paste0(getMainInputDir(),"/preData.rds")))
            }
            return(NULL)
        }
        ### if input a file
        longDataframe <- getMainInput()
        # isolate start and end gene index
        input$applyFilter
        if (input$autoUpdate == TRUE) {
            startIndex <- input$stIndex
            endIndex <- input$endIndex
        } else {
            startIndex <- isolate(input$stIndex)
            endIndex <- isolate(input$endIndex)
        }
        if (is.na(endIndex)) endIndex <- 1000

        withProgress(message = 'Subseting data...', value = 0.5, {
            longDataframe <- sortGeneIDs(
                longDataframe, input$orderGenes, checkInputSortedGenes()
            )
            # filter preData based on DIM reduction selection
            if (!v$doPlot && (isTruthy(input$addSpecDimRed) || isTruthy(input$addGeneDimRed))) {
                selectedTaxa <- longDataframe$ncbiID
                if (input$addSpecDimRed) {
                    dimRedTaxa <- dimRedSelectedTaxa()
                    selectedTaxa <- unique(head(dimRedTaxa$`NCBI taxon ID`))
                }
                dimRedGenes <- dimRedSelectedGenes()
                selectedGenes <- unique(dimRedGenes$geneID)

                if (length(selectedGenes) > 0) {
                    data <- longDataframe %>%
                        filter(geneID %in% selectedGenes & ncbiID %in% selectedTaxa) %>%
                        droplevels()
                } else {
                    return(NULL)
                }
            } else {
                # Handle custom gene list or subset ID
                geneList <- input$geneList
                if (!is.null(geneList)) {
                    geneListDf <- read.table(file = geneList$datapath, header = FALSE)
                    listGeneOri <- unique(geneListDf$V1)

                    listGene <- listGeneOri[startIndex:min(endIndex, length(listGeneOri))]
                    data <- longDataframe[longDataframe$geneID %in% listGene[!is.na(listGene)], ]
                } else {
                    subsetID <- levels(longDataframe$geneID)[startIndex:endIndex]
                    data <- longDataframe[longDataframe$geneID %in% subsetID, ]
                }
            }
            # ensure a minimum of 5 columns
            while (ncol(data) < 5) {
                data[paste0("newVar", ncol(data) + 1)] <- 1
            }
            # return preData
            if (nrow(data) == 0) return()
            if (ncol(data) < 6) {
                colnames(data) <- c("geneID","ncbiID","orthoID","var1","var2")
            } else {
                colnames(data) <- c(
                    "geneID", "ncbiID", "orthoID", "var1", "var2", "geneName"
                )
            }
            data$geneID <- droplevels(data$geneID)
            if (rtCheck) {
                checkpoint72 <- Sys.time()
                print(paste(
                    "checkpoint72 - subseting input done", checkpoint72, 
                    " --- ",  checkpoint72 - checkpoint71
                ))
            }
            return(data)
        })
    })

    # * creating main dataframe for subset taxa (in species/strain level) ------
    getFullData <- reactive({
        req(isTruthy(v$doPlot)|isTruthy(w$doCusPlot))
        req(input$mainInputType)
        if (rtCheck) {
            checkpoint81 <- Sys.time()
            print(paste("checkpoint81 - before get full data", checkpoint81))
        }
        if (input$mainInputType == "folder") {
            req(getMainInputDir())
            if (length(checkMainInputDir()) == 4) {
                fullMdData <- readRDS(
                    paste0(getMainInputDir(),"/fullData.rds")
                )
                return(fullMdData)
            } else return()
        } else {
            req(preData())
            req(getCountTaxa())
            req(sortedtaxaList())
            {
                input$applyFilter
            }
            withProgress(message = 'Parsing profile data...', value = 0.5, {
                if (input$autoUpdate == TRUE) {
                    coorthologCutoffMax <- input$coortholog
                } else {
                    coorthologCutoffMax <- isolate(input$coortholog)
                }
                if (rtCheck) {
                    checkpoint81a <- Sys.time()
                    print(paste(
                        "checkpoint81a - start get full data", checkpoint81a, 
                        " --- ", checkpoint81a- checkpoint81
                    ))
                }
                fullMdData <- parseInfoProfile(
                    inputDf = preData(),
                    sortedInputTaxa = sortedtaxaList(),
                    taxaCount = getCountTaxa(),
                    coorthoCOMax = coorthologCutoffMax
                )
                if (rtCheck) {
                    checkpoint82 <- Sys.time()
                    print(paste(
                        "checkpoint82 - get full data done", checkpoint82, 
                        " --- ",  checkpoint82- checkpoint81a
                    ))
                }
                return(fullMdData)
            })
        }
    })

    # * filter full data -------------------------------------------------------
    filteredDataHeat <- reactive({
        req(isTruthy(v$doPlot)|isTruthy(w$doCusPlot))
        {
            input$applyFilterCustom
            input$applyFilter
        }
        if (rtCheck) {
            checkpoint91 <- Sys.time()
            print(paste("checkpoint91 - before filter full data", checkpoint91))
        }
        withProgress(message = 'Creating data for plotting...', value = 0.5, {
            # check input file
            req(input$mainInputType)
            if (input$mainInputType == "folder" | input$demoData == "arthropoda" |
                input$demoData == "ampk-tor" | input$demoData == "preCalcDt") {
                filein <- 1
            } else {
                filein <- input$mainInput
            }
            req(filein)

            # get all cutoffs
            if (input$autoUpdate == TRUE) {
                percentCutoff <- input$percent
                coorthologCutoffMax <- input$coortholog
                var1Cutoff <- input$var1
                var2Cutoff <- input$var2
                colorByGroup <- input$colorByGroup
            } else {
                percentCutoff <- isolate(input$percent)
                coorthologCutoffMax <- isolate(input$coortholog)
                var1Cutoff <- isolate(input$var1)
                var2Cutoff <- isolate(input$var2)
                colorByGroup <- isolate(input$colorByGroup)
            }

            # get selected supertaxon name
            split <- strsplit(as.character(input$inSelect), "_")
            inSelect <- as.character(split[[1]][1])

            # get gene categories
            inputCatDt <- NULL
            if (length(colorByGroup) > 0 && colorByGroup == TRUE) {
                # get gene category
                geneCategoryFile <- input$geneCategory
                if (!is.null(geneCategoryFile)) {
                    inputCatDt <- read.table(
                        file = geneCategoryFile$datapath,
                        sep = "\t",
                        header = FALSE,
                        check.names = FALSE,
                        comment.char = "",
                        fill = TRUE
                    )
                    colnames(inputCatDt) <- c("geneID","group")
                } else if (!is.null(i_geneCategory)){
                    inputCatDt <- read.table(
                        file = i_geneCategory,
                        sep = "\t",
                        header = FALSE,
                        check.names = FALSE,
                        comment.char = "",
                        fill = TRUE
                    )
                    colnames(inputCatDt) <- c("geneID","group")
                } else inputCatDt <- NULL
            } else colorByGroup = FALSE

            # create data for heatmap plotting
            if (rtCheck) {
                checkpoint91a <- Sys.time()
                print(paste(
                    "checkpoint91a - start filter full data", checkpoint91a, 
                    " --- ", checkpoint91a- checkpoint91
                ))
            }
            filteredDf <- filterProfileData(
                DF = getFullData(),
                taxaCount = getCountTaxa(),
                refTaxon = inSelect,
                percentCutoff,
                coorthologCutoffMax,
                var1Cutoff,
                var2Cutoff,
                input$var1Relation,
                input$var2Relation,
                groupByCat = colorByGroup,
                catDt = inputCatDt,
                var1AggregateBy = input$var1AggregateBy,
                var2AggregateBy = input$var2AggregateBy
            )
            if (rtCheck) {
                checkpoint92 <- Sys.time()
                print(paste(
                    "checkpoint92 - filter full data done", checkpoint92, 
                    " --- ",  checkpoint92- checkpoint91a
                ))
            }
            return(filteredDf)
        })
    })

    # * heatmap data input -----------------------------------------------------
    dataHeat <- reactive({
        req(isTruthy(v$doPlot)|isTruthy(w$doCusPlot))
        req(filteredDataHeat())
        if (rtCheck) {
            checkpoint101 <- Sys.time()
            print(paste("checkpoint101 - before reduce profile", checkpoint101))
        }
        dataHeat <- reduceProfile(filteredDataHeat())
        if (rtCheck) {
            checkpoint102 <- Sys.time()
            print(paste(
                "checkpoint102 - reduce profile done", checkpoint102, " --- ",
                checkpoint102- checkpoint101
            ))
        }
        return(dataHeat)
    })

    # * clustered heatmap data -------------------------------------------------
    clusteredDataHeat <- reactive({
        req(isTruthy(v$doPlot)|isTruthy(w$doCusPlot))
        if (rtCheck) {
            checkpoint111 <- Sys.time()
            print(paste(
                "checkpoint111 - before cluster heatmap data", checkpoint111
            ))
        }
        dataHeat <- dataHeat()
        if (nlevels(dataHeat$supertaxon) == 1) return(dataHeat)
        if (!is.null(i_clusterMethod)) clusterMethod <- i_clusterMethod
        else clusterMethod <- input$clusterMethod

        if (nlevels(as.factor(dataHeat$geneID)) > 1) {
            withProgress(message = 'Clustering profile data...', value = 0.5, {
                dat <- getProfiles()
                # do clustering based on distance matrix
                row.order <- fastcluster::hclust(
                    as.dist(getDistanceMatrixProfiles()), method = clusterMethod
                )$order

                # re-order distance matrix accoring to clustering
                datNew <- dat[row.order, ]

                # return clustered gene ID list
                clusteredGeneIDs <- as.factor(row.names(datNew))

                # sort original data according to clusteredGeneIDs
                dataHeat$geneID <- factor(
                    dataHeat$geneID, levels = clusteredGeneIDs
                )
                orderedName <- unique(
                    dataHeat[
                        order(match(dataHeat$geneID, clusteredGeneIDs)),
                        "geneName"
                    ]
                )
                dataHeat$geneName <- factor(
                    dataHeat$geneName, levels = orderedName
                )
                dataHeat <- dataHeat[!is.na(dataHeat$geneID),]
                if (rtCheck) {
                    checkpoint112 <- Sys.time()
                    print(paste(
                        "checkpoint112 - cluster heatmap data done", 
                        checkpoint112, " --- ",  checkpoint112 - checkpoint111
                    ))
                }
                return(dataHeat)
            })
        } else return(dataHeat)
    })

    # * switch to fast mode for large data -------------------------------------
    observe({
        dt <- getMainInput()
        if (
            nlevels(as.factor(dt$geneID)) > fastModeCutoff |
            nlevels(as.factor(dt$ncbiID)) > fastModeCutoff
        )
            updateRadioButtons(
                session, "plotMode",
                choices = list("Normal" = "normal", "Fast" = "fast"),
                selected = "fast", inline = TRUE
            )
    })

    # =========================== MAIN PROFILE TAB =============================

    # * update choices for geneIdType ------------------------------------------
    observe({
        if (!(is.null(getGeneNames())))
            updateRadioButtons(session, "geneIdType", selected = "geneName")
    })

    # * render popup for selecting rank and return list of subset taxa ---------
    mainTaxaName <- callModule(
        selectTaxonRank,
        "selectTaxonRankMain",
        rankSelect = reactive(input$rankSelect),
        inputTaxonID = inputTaxonID,
        taxDB = getTaxDBpath
    )

    # * get list of taxa for highlighting --------------------------------------
    # output$taxonHighlight.ui <- renderUI({
    observe({
        req(getMainInput())
        req(v$doPlot)
        if (rtCheck) {
            checkpoint131 <- Sys.time()
            print(paste(
                "checkpoint131 - before render taxa hghlight", checkpoint131
            ))
        }
        filein <- input$mainInput
        if (
            input$demoData == "arthropoda" | input$demoData == "ampk-tor" |
            input$demoData == "preCalcDt"
        ) {
            filein <- 1
        }
        choice <- inputTaxonName()
        out <- levels(factor(choice$fullName))
        if (input$applyMainTaxa == TRUE) {
            outSub <- mainTaxaName()
            updateSelectizeInput(
                session, "taxonHighlight", server = TRUE,
                choices = out, selected = outSub
            )
        } else {
            updateSelectizeInput(
                session, "taxonHighlight", server = TRUE,
                choices = out
            )
        }
        if (rtCheck) {
            checkpoint132 <- Sys.time()
            print(paste(
                "checkpoint132 - render tax highlight done", checkpoint132, 
                " --- ",  checkpoint132 - checkpoint131
            ))
        }
    })

    # * get list of genes for highlighting -------------------------------------
    getAllGenes <- function() {
        req(getMainInput())
        req(v$doPlot)
        df <- getMainInput()
        if (rtCheck) {
            checkpoint141 <- Sys.time()
            print(paste(
                "checkpoint141 - before render gene highlight", checkpoint141
            ))
        }
        idNameDf <- df %>% dplyr::select(geneID, geneName) %>% distinct()
        idNameList <- setNames(
            as.character(idNameDf$geneID), as.character(idNameDf$geneName)
        )
        if (rtCheck) {
            checkpoint142 <- Sys.time()
            print(paste(
                "checkpoint142 - render gene highlight done", checkpoint142, 
                " --- ",  checkpoint142 - checkpoint141
            ))
        }
        return(idNameList)
    }

    observe({
        out <- getAllGenes()
        if (!(is.null(input$geneHighlightFile))) {
            fileHighlight <- input$geneHighlightFile
            highlightList <- read.table(
                file = fileHighlight$datapath, header = FALSE
            )
            updateSelectizeInput(
                session, "geneHighlight", server = TRUE,
                choices = out, selected = intersect(out, highlightList$V1)
            )
        } else {
            updateSelectizeInput(
                session, "geneHighlight", server = TRUE,
                choices = out
            )
        }
    })

    # * disable/enable highlighing orthologs having the same ID ----------------
    observe({
        longDataframe <- getMainInput()
        req(longDataframe)
        req(input$rankSelect)
        lowestRank <- getLowestRank(longDataframe, getTaxDBpath())
        if (!(lowestRank == input$rankSelect))
            shinyjs::disable("colorByOrthoID")
        else
            shinyjs::enable("colorByOrthoID")
    })

    # * render list of superRanks for adding vertical lines --------------------
    output$superRankSelect.ui <- renderUI({
        allRanks <- getTaxonomyRanks()
        selectInput(
            "superRankSelect", label = "Display taxonomic labels for:",
            choices = c(allRanks[!(allRanks %in% c("strain"))]),
            selected = ""
        )
    })

    # * update plot size based on input ----------------------------------------
    observe({
        req(getMainInput())
        longDataframe <- getMainInput()
        if (rtCheck) {
            checkpoint151 <- Sys.time()
            print(paste("checkpoint151 - before update plot size", checkpoint151))
        }
        req(longDataframe)
        req(input$endIndex)
        req(input$plotMode)
        if (input$autoSizing & input$plotMode == "normal") {
            inputSuperTaxon <- inputTaxonName()
            nrTaxa <- nlevels(as.factor(inputSuperTaxon$fullName))
            nrGene <- input$endIndex

            adaptedSize <- adaptPlotSize(
                nrTaxa, nrGene, input$xAxis, input$dotZoom
            )
            req(length(adaptedSize) > 0)
            h <- adaptedSize[1]
            hv <- adaptedSize[2]
            wv <- adaptedSize[3]
            if (h <= 20) {
                updateSelectInput(
                    session, "mainLegend",
                    label = "Legend position:",
                    choices = list("Right" = "right",
                                   "Left" = "left",
                                   "Top" = "top",
                                   "Bottom" = "bottom",
                                   "Hide" = "none"),
                    selected = "top"
                )
                updateNumericInput(session, "width", value = min(20000, wv  + 50))
            } else if (h <= 30) {
                updateNumericInput(session, "width", value = min(20000, wv + 50))
            } else {
                updateNumericInput(session, "width", value = min(20000, wv))
            }
            updateNumericInput(session, "height", value = min(20000, hv))
        }
        if (rtCheck) {
            checkpoint152 <- Sys.time()
            print(paste(
                "checkpoint152 - update plot size done", checkpoint152, 
                " --- ",  checkpoint152 - checkpoint151
            ))
        }
    })

    observe({
        if(!(is.null(input$plotMode))) {
            if (input$plotMode == "fast") {
                updateNumericInput(session, "height", value = 900)
                updateNumericInput(session, "width", value = 900)
                updateNumericInput(session, "selectedHeight", value = 900)
                updateNumericInput(session, "selectedWidth", value = 900)
            }
        }
    })

    # * reset configuration windows of Main plot -------------------------------
    observeEvent(input$resetMainConfig, {
        shinyjs::reset("xSize")
        shinyjs::reset("ySize")
        shinyjs::reset("legendSize")
        shinyjs::reset("xAngle")
        shinyjs::reset("dotZoom")
    })

    # * close configuration windows of Main plot -------------------------------
    observeEvent(input$applyMainConfig, {
        toggleModal(session, "mainPlotConfigBs", toggle = "close")
    })

    # * parameters for the main profile plot -----------------------------------
    getParameterInputMain <- reactive({
        input$updateBtn
        if (rtCheck) {
            checkpoint161 <- Sys.time()
            print(paste(
                "checkpoint161 - before get main plot para", checkpoint161
            ))
        }
        colorByGroup <- input$colorByGroup
        # get category colors
        catColors <- NULL
        if (length(colorByGroup) > 0 && colorByGroup == TRUE) {
            geneCategoryFile <- input$geneCategory
            if (!is.null(geneCategoryFile)) {
                catColors <- getCatColors(geneCategoryFile, type = "file")
            } else if (!is.null(i_geneCategory)){
                catColors <- getCatColors(i_geneCategory, type = "config")
            }
        } else colorByGroup == FALSE
        if (input$autoUpdate == TRUE) {
            inputPara <- list(
                "xAxis" = input$xAxis,
                "geneIdType" = input$geneIdType,
                "var1ID" = input$var1ID,
                "var2ID"  = input$var2ID,
                "midVar1" = input$midVar1,
                "midVar2" = input$midVar2,
                "lowColorVar1" =  input$lowColorVar1,
                "midColorVar1" =  input$midColorVar1,
                "highColorVar1" = input$highColorVar1,
                "lowColorVar2" = input$lowColorVar2,
                "midColorVar2" =  input$midColorVar2,
                "highColorVar2" = input$highColorVar2,
                "paraColor" = input$paraColor,
                "xSize" = input$xSize,
                "ySize" = input$ySize,
                "legendSize" = input$legendSize,
                "mainLegend" = input$mainLegend,
                "dotZoom" = input$dotZoom,
                "xAngle" = input$xAngle,
                "font" = input$font,
                "guideline" = 0,
                "width" = input$width,
                "height" = input$height,
                "colorByGroup" = colorByGroup,
                "catColors" = catColors,
                "colorByOrthoID" = input$colorByOrthoID,
                "groupLabelSize" = input$groupLabelSize,
                "groupLabelDist" = input$groupLabelDist,
                "groupLabelAngle" = input$groupLabelAngle,
                "colorVar" = input$colorVar
            )
        } else {
            inputPara <- isolate(
                list(
                    "xAxis" = input$xAxis,
                    "geneIdType" = input$geneIdType,
                    "var1ID" = input$var1ID,
                    "var2ID"  = input$var2ID,
                    "midVar1" = input$midVar1,
                    "midVar2" = input$midVar2,
                    "lowColorVar1" =  input$lowColorVar1,
                    "midColorVar1" =  input$midColorVar1,
                    "highColorVar1" = input$highColorVar1,
                    "lowColorVar2" = input$lowColorVar2,
                    "midColorVar2" =  input$midColorVar2,
                    "highColorVar2" = input$highColorVar2,
                    "paraColor" = input$paraColor,
                    "font" = input$font,
                    "xSize" = input$xSize,
                    "ySize" = input$ySize,
                    "legendSize" = input$legendSize,
                    "mainLegend" = input$mainLegend,
                    "dotZoom" = input$dotZoom,
                    "xAngle" = input$xAngle,
                    "guideline" = 0,
                    "width" = input$width,
                    "height" = input$height,
                    "colorByGroup" = colorByGroup,
                    "catColors" = catColors,
                    "colorByOrthoID" = input$colorByOrthoID,
                    "groupLabelSize" = input$groupLabelSize,
                    "groupLabelDist" = input$groupLabelDist,
                    "groupLabelAngle" = input$groupLabelAngle,
                    "colorVar" = input$colorVar
                )
            )
        }
        if (rtCheck) {
            checkpoint162 <- Sys.time()
            print(paste(
                "checkpoint162 - get plot main para done", checkpoint162, 
                " --- ",  checkpoint162 - checkpoint161
            ))
        }
        return(inputPara)
    })

    # * render dot size to dotSizeInfo ---------------------------------------
    output$dotSizeInfo <- renderUI({
        req(v$doPlot)

        dataHeat <- dataHeat()
        dataHeat$presSpec[dataHeat$presSpec == 0] <- NA
        presentVl <- dataHeat$presSpec[!is.na(dataHeat$presSpec)]

        minDot <- (floor(min(presentVl) * 10) / 10 * 5) * (1 + input$dotZoom)
        maxDot <- (floor(max(presentVl) * 10) / 10 * 5) * (1 + input$dotZoom)

        em(paste0("current point's size: ", minDot, " - ", maxDot))
    })

    # * plot main profile ------------------------------------------------------
    mainpointInfo <- callModule(
        createProfilePlot, "mainProfile",
        data = dataHeat,
        clusteredDataHeat = clusteredDataHeat,
        applyCluster = reactive(input$applyCluster),
        parameters = getParameterInputMain,
        inSeq = reactive(input$inSeq),
        inTaxa = reactive(input$inTaxa),
        rankSelect = reactive(input$rankSelect),
        inSelect = reactive(input$inSelect),
        taxonHighlight = reactive(input$taxonHighlight),
        geneHighlight = reactive(input$geneHighlight),
        typeProfile = reactive("mainProfile"),
        taxDB = getTaxDBpath,
        superRank = reactive(input$superRankSelect),
        allTaxa = allInputTaxa,
        mode = reactive(input$plotMode)
    )

    # ======================== CUSTOMIZED PROFILE TAB ==========================
    observe({
        if ("all" %in% input$inSeq & "all" %in% input$inTaxa) {
            shinyjs::disable("plotCustom")
            # shinyjs::disable("applyFilterCustom")
        } else {
            shinyjs::enable("plotCustom")
            # shinyjs::enable("applyFilterCustom")
        }
    })

    # * get list of all sequence IDs for customized profile --------------------
    observe({
        fileCustom <- input$customFile
        if (v$doPlot == FALSE) {
            if (input$addGeneDimRed == TRUE) {
                dimRedGenes <- dimRedSelectedGenes()
                req(dimRedGenes)
                if (nrow(dimRedGenes) > 0)
                    return(updateSelectizeInput(
                        session, "inSeq", server = TRUE, "",
                        unique(dimRedGenes$geneID),
                        selected = unique(dimRedGenes$geneID)
                    ))
            }
        } else {
            idNameList <- getAllGenes()
            outAll <- c(setNames("all", "all"), idNameList)
            if (input$addGeneAgeCustomProfile == TRUE) {
                outAll <- as.list(selectedgeneAge())
                outAll <- outAll[[1]]
            } else if (input$addClusterCustomProfile == TRUE) {
                outAll <- as.list(brushedClusterGene())
            } else if (input$addCoreGeneCustomProfile == TRUE) {
                outAll <- as.list(coreGeneDf())
            } else if (input$addGCGenesCustomProfile == TRUE) {
                outAll <- as.list(candidateGenes())
            } else if (input$addGeneDimRed == TRUE) {
                req(dimRedSelectedGenes())
                dimRedGenes <- dimRedSelectedGenes()
                if (nrow(dimRedGenes) > 0) outAll <- unique(dimRedGenes$geneID)
            } else {
                if (!is.null(fileCustom)) {
                    customList <- read.table(
                        file = fileCustom$datapath, header = FALSE
                    )
                    customList$V1 <- as.factor(customList$V1)
                    outAll <- as.list(levels(customList$V1))
                } else {
                    return(updateSelectizeInput(
                        session, "inSeq", server = TRUE, "", outAll,
                        selected = "all"
                    ))
                }
            }
            if (length(outAll) == 0) outAll <- c("all")
            if (outAll[1] == "all") {
                updateSelectizeInput(
                    session, "inSeq", server = TRUE,"", outAll, selected = "all"
                )
            } else {
                updateSelectizeInput(
                    session, "inSeq", server = TRUE, "", outAll, selected=outAll
                )
            }
        }
    })

    # * render popup for selecting rank and return list of subset taxa ---------
    cusTaxaName <- callModule(
        selectTaxonRank,
        "selectTaxonRank",
        rankSelect = reactive(input$rankSelect),
        inputTaxonID = inputTaxonID,
        taxDB = getTaxDBpath
    )

    # * get list of all taxa for customized profile ----------------------------
    output$cusTaxa.ui <- renderUI({
        if (
            input$demoData == "arthropoda" | input$demoData == "ampk-tor" |
            input$demoData == "preCalcDt"  | input$mainInputType == "folder"
        ) {
            filein <- 1
        } else filein <- input$mainInput

        if (is.null(filein) & input$addSpecDimRed != TRUE)
            return(selectInput("inTaxa", "", "all"))
        if (v$doPlot == FALSE) {
            if (input$addSpecDimRed == TRUE) {
                dimRedTaxa <- dimRedSelectedTaxa()
                if (is.null(dimRedTaxa))
                    return(selectInput("inTaxa", "", "all"))
                if (nrow(dimRedTaxa) > 0) out <- unique(dimRedTaxa$`Taxon name`)
                selectizeInput("inTaxa","",out, selected = out, multiple = TRUE)
            } else return(selectInput("inTaxa", "", "all"))
        } else {
            choice <- inputTaxonName()
            out <- c("all", as.list(levels(factor(choice$fullName))))
            if (input$applyCusTaxa == TRUE) {
                out <- cusTaxaName()
                selectizeInput("inTaxa","",out, selected = out, multiple = TRUE)
            } else if (input$addSpecDimRed == TRUE) {
                dimRedTaxa <- dimRedSelectedTaxa()
                if (is.null(dimRedTaxa))
                    return(selectInput("inTaxa", "", "all"))
                if (nrow(dimRedTaxa) > 0) out <- unique(dimRedTaxa$`Taxon name`)
                selectizeInput("inTaxa","",out, selected = out, multiple = TRUE)
            } else {
                selectizeInput("inTaxa","",out,selected = out[1], multiple=TRUE)
            }
        }
    })

    # * render list of superRanks for adding vertical lines --------------------
    output$cusSuperRankSelect.ui <- renderUI({
        allRanks <- getTaxonomyRanks()
        selectInput(
            "cusSuperRankSelect", label = "Display taxonomic labels for:",
            choices = c(allRanks[!(allRanks %in% c("strain"))]),
            selected = ""
        )
    })

    # * check if all genes and all species are selected ------------------------
    output$sameProfile <- reactive({
        if (v$doPlot == FALSE) return(FALSE)
        if (length(input$inSeq[1]) == 0) return(FALSE)
        else {
            if (length(input$inSeq) == 0 || length(input$inTaxa) == 0)
                return(TRUE)
            if ("all" %in% input$inSeq & "all" %in% input$inTaxa) return(TRUE)
        }
    })
    outputOptions(output, "sameProfile", suspendWhenHidden = FALSE)

    # * update customized plot size based on input -----------------------------
    observe({
        longDataframe <- getMainInput()
        req(longDataframe)
        req(input$inTaxa)
        req(input$inSeq)
        if (input$selectedAutoSizing) {
            nrTaxa <- length(input$inTaxa)
            nrGene <- length(input$inSeq)
            if (nrTaxa < 10000 && nrGene < 10000) {
                if ("all" %in% input$inTaxa) {
                    inputSuperTaxon <- inputTaxonName()
                    nrTaxa <- nlevels(as.factor(inputSuperTaxon$fullName))
                }
                if ("all" %in% input$inSeq) {
                    nrGene <- input$endIndex
                }
                # adapte to axis type
                if (input$xAxisSelected == "taxa") {
                    h <- nrGene
                    w <- nrTaxa
                } else {
                    w <- nrGene
                    h <- nrTaxa
                }
                # adapt to dot zoom factor
                if (input$dotZoomSelect < -0.5){
                    hv <- (200 + 12 * h) * (1 + input$dotZoomSelect) + 500
                    wv <- (200 + 12 * w) * (1 + input$dotZoomSelect) + 500
                }  else if ((input$dotZoomSelect < 0)) {
                    hv <- (200 + 12 * h) * (1 + input$dotZoomSelect) + 200
                    wv <- (200 + 12 * w) * (1 + input$dotZoomSelect) + 200
                } else {
                    hv <- (200 + 12 * h) * (1 + input$dotZoomSelect)
                    wv <- (200 + 12 * w) * (1 + input$dotZoomSelect)
                }
                # minimum size
                if (hv < 300) hv <- 300
                if (wv < 300) wv <- 300
                # update plot size based on number of genes/taxa
                hv <- hv + 300
                wv <- wv + 300
                if (h <= 20) {
                    updateSelectInput(
                        session, "selectedLegend",
                        label = "Legend position:",
                        choices = list("Right" = "right",
                                       "Left" = "left",
                                       "Top" = "top",
                                       "Bottom" = "bottom",
                                       "Hide" = "none"),
                        selected = "top"
                    )
                    updateNumericInput(
                        session,
                        "selectedWidth", value = min(20000, wv  + 50)
                    )
                } else if (h <= 30) {
                    updateNumericInput(
                        session,
                        "selectedWidth", value = min(20000, wv + 50)
                    )
                } else {
                    updateNumericInput(
                        session,
                        "selectedWidth", value = min(20000, wv)
                    )
                }
                updateNumericInput(
                    session,
                    "selectedHeight", value = min(20000, hv)
                )
            }
        }
    })

    # * reset configuration windows of Customized plot -------------------------
    observeEvent(input$resetSelectedConfig, {
        shinyjs::reset("xSizeSelect")
        shinyjs::reset("ySizeSelect")
        shinyjs::reset("legendSizeSelect")
        shinyjs::reset("xAngleSelect")
        shinyjs::reset("dotZoomSelect")
    })

    # ** close configuration windows of Customized plot ------------------------
    observeEvent(input$applySelectedConfig, {
        toggleModal(session, "selectedPlotConfigBs", toggle = "close")
    })

    # * parameters for the customized profile plot -----------------------------
    getParameterInputCustomized <- reactive({
        input$plotCustom

        colorByGroup <- input$colorByGroup
        # get category colors
        catColors <- NULL
        if (length(colorByGroup) > 0 && colorByGroup == TRUE) {
            geneCategoryFile <- input$geneCategory
            if (!is.null(geneCategoryFile)) {
                catColors <- getCatColors(geneCategoryFile, type = "file")
            } else if (!is.null(i_geneCategory)){
                catColors <- getCatColors(i_geneCategory, type = "config")
            }
        } else colorByGroup = FALSE

        inputPara <- isolate(
            list(
                "xAxis" = input$xAxisSelected,
                "geneIdType" = input$geneIdType,
                "var1ID" = input$var1ID,
                "var2ID"  = input$var2ID,
                "midVar1" = input$midVar1,
                "midVar2" = input$midVar2,
                "lowColorVar1" =  input$lowColorVar1,
                "midColorVar1" =  input$midColorVar1,
                "highColorVar1" = input$highColorVar1,
                "lowColorVar2" = input$lowColorVar2,
                "midColorVar2" =  input$midColorVar2,
                "highColorVar2" = input$highColorVar2,
                "paraColor" = input$paraColor,
                "xSize" = input$xSizeSelect,
                "ySize" = input$ySizeSelect,
                "legendSize" = input$legendSizeSelect,
                "mainLegend" = input$selectedLegend,
                "dotZoom" = input$dotZoomSelect,
                "xAngle" = input$xAngleSelect,
                "font" = input$font,
                "guideline" = 0,
                "width" = input$selectedWidth,
                "height" = input$selectedHeight,
                "colorByGroup" = colorByGroup,
                "catColors" = catColors,
                "colorByOrthoID" = input$colorByOrthoID,
                "groupLabelSize" = input$groupLabelSizeSelect,
                "groupLabelDist" = input$groupLabelDistSelect,
                "groupLabelAngle" = input$groupLabelAngleSelect
            )
        )
        return(inputPara)
    })

    # * plot customized profile ------------------------------------------------
    selectedpointInfo <- callModule(
        createProfilePlot, "customizedProfile",
        data = dataHeat,
        clusteredDataHeat = clusteredDataHeat,
        applyCluster = reactive(input$applyCluster),
        parameters = getParameterInputCustomized,
        inSeq = reactive(input$inSeq),
        inTaxa = reactive(input$inTaxa),
        rankSelect = reactive(input$rankSelect),
        inSelect = reactive(input$inSelect),
        taxonHighlight = reactive("none"),
        geneHighlight = reactive("none"),
        typeProfile = reactive("customizedProfile"),
        taxDB = getTaxDBpath,
        superRank = reactive(input$cusSuperRankSelect),
        allTaxa = allInputTaxa,
        mode = reactive(input$plotMode)
    )

    # ==================== DIMENSIONALITY REDUCTION PLOT =======================
    shinyjs::disable("addSpecDimRed")
    observe({
        if (input$addSpecDimRed) shinyjs::disable("addGeneDimRed")
        else shinyjs::enable("addGeneDimRed")
    })

    u <- reactiveValues(doCusPlot = FALSE)
    observeEvent(input$plotDimRed, {
        u$doDimRedPlot <- input$plotDimRed
        filein <- input$mainInput
        if (
            input$mainInputType == "file" & is.null(filein) &
            input$demoData == "none"
        ) {
            u$doDimRedPlot <- FALSE
            shinyBS::updateButton(session, "plotDimRed", disabled = TRUE)
        }
        shinyBS::updateButton(session, "plotDimRed", "Update plot")
    })

    # * toggle dimRedRank based on dimRedType genes or taax --------------
    observe({
        if (input$dimRedType == "genes") {
            shinyjs::disable("dimRedRank")
            shinyjs::disable("dimRedGroupHigherRank")
            shinyjs::disable("dimRedApplyChangeLables")
            shinyjs::disable("dimRedResetLables")
            updateRadioButtons(
                session, "dimRedGroupLabelsBy", choices = c("taxa"), inline = TRUE
            )
        } else {
            shinyjs::enable("dimRedRank")
            shinyjs::enable("dimRedGroupHigherRank")
            shinyjs::enable("dimRedApplyChangeLables")
            shinyjs::enable("dimRedResetLables")
            updateRadioButtons(
                session, "dimRedGroupLabelsBy", choices = c("taxa", "genes"),
                inline = TRUE
            )
        }
    })

    # * update variables for data filtering ------------------------------------
    observe({
        req(getMainInput())
        inputDf <- getMainInput()
        if (ncol(inputDf) == 5) {
            choiceList <- c("var1", "var2", "both")
            names(choiceList) <- c(
                colnames(inputDf)[4], colnames(inputDf)[5], "Both"
            )
            updateSelectInput(
                session, "dimRedFilterVar", choices = choiceList, selected="both"
            )
        } else if (ncol(inputDf) == 4) {
            choiceList <- "var1"
            names(choiceList) <- colnames(inputDf)[4]
            updateSelectInput(session, "dimRedFilterVar", choices = choiceList)
        } else if (ncol(inputDf) == 3) {
            shinyjs::disable("dimRedFilterVar")
            shinyjs::disable("dimRedCutoff")
            shinyjs::disable("dimRedDataType")
        }
    })

    # * data for DIM reduction clustering --------------------------------------
    dimRedData <- reactive({
        req(getMainInput())
        req(u$doDimRedPlot)
        if(is.null(getMainInput())) stop("Input data is NULL!")
        
        input$plotDimRed
        dimRedParams <- isolate(list(
            rank = input$dimRedRank,
            type = input$dimRedType,
            filterVar = input$dimRedFilterVar,
            cutoff = input$dimRedCutoff,
            groupLabelsBy = input$dimRedGroupLabelsBy
        ))
        
        withProgress(
            message = "Preparing data for clustering...", value = 0.5, {
                dimRedData <- prepareDimRedData(
                    getMainInput(), dimRedParams$rank, dimRedParams$type, 
                    getTaxDBpath(), dimRedParams$filterVar, dimRedParams$cutoff, 
                    dimRedParams$groupLabelsBy
                )
                return(dimRedData)
            }
        )
    })

    # * read custom labels -----------------------------------------------------
    values <- reactiveValues(
        uploadLabelState = NULL
    )

    output$dimRedCustomLabel.ui <- renderUI({
        if (input$dimRedType == "taxa")
            fileInput("dimRedCustomLabel", "Add customized labels")
    })

    observeEvent(input$input$dimRedCustomLabel, {
        values$uploadLabelState <- 'uploaded'
    })

    observeEvent(input$dimRedResetLables, {
        values$uploadLabelState <- NULL
        updateTextInput(session, "dimRedGroupHigherRank", value = "")
    })

    getCustomLabels <- reactive({
        req(getMainInput())
        if (is.null(values$uploadLabelState))
            return(data.frame(ncbiID = c(), label = c()))
        filein <- input$dimRedCustomLabel
        if (!is.null(filein)) {
            customLabels <- read.table(
                file = filein$datapath, header = FALSE, check.names = FALSE,
                comment.char = "", fill = FALSE
            )
            colnames(customLabels) <- c("ncbiID", "Label")
            mainInput <- getMainInput()
            return(customLabels[customLabels$ncbiID %in% mainInput$ncbiID,])
        } else return(data.frame(ncbiID = c(), label = c()))
    })

    # * apply user-defined labels ----------------------------------------------
    renameLabelsDimRed <- reactive({
        req(dimRedData())
        input$dimRedApplyChangeLables
        dimRedData <- dimRedData()
        isolate({
            if (input$dimRedType == "taxa") {
                # group labels into higher rank
                higherRankTaxa <- unlist(strsplit(input$dimRedGroupHigherRank, ";"))
                higherRankTaxa <- trimws(higherRankTaxa)
                if (length(higherRankTaxa) > 0) {
                    taxMatrix <- getTaxonomyMatrix(getTaxDBpath())
                    nameList <- getNameList(getTaxDBpath())
                    selDf <- data.frame(
                        selRank = nameList$rank[
                            nameList$fullName %in% higherRankTaxa],
                        selID = nameList$ncbiID[
                            nameList$fullName %in% higherRankTaxa],
                        Label = nameList$fullName[
                            nameList$fullName %in% higherRankTaxa]
                    )
                    selDf <- selDf[complete.cases(selDf),]
                    if (nrow(selDf) > 0) {
                        selTaxList <- lapply(
                            seq(nrow(selDf)), function (x) {
                                selRank <- selDf$selRank[x]
                                selID <- selDf$selID[x]
                                if (!(selRank %in% mainTaxonomyRank()))
                                    selRank <- paste0("norank_", selID)
                                selRank <- quo(!! sym(selRank))
                                selTaxDf <- taxMatrix %>%
                                    filter((!!selRank) %in% selID) %>%
                                    select(abbrName, !!selRank)
                                colnames(selTaxDf) <- c("ncbiID", "supertaxonID")
                                selTaxDf$sel_label <- selDf$Label[selDf$selID == selID]
                                return(selTaxDf)
                            }
                        )
                        joinedSelTaxDf <- do.call(rbind, selTaxList)
                        joinedSelTaxDf <- joinedSelTaxDf %>% group_by(ncbiID) %>%
                            filter(supertaxonID == min(supertaxonID))
                        dimRedData <- left_join(
                            dimRedData, joinedSelTaxDf, by = "ncbiID"
                        ) %>% mutate(
                            Label = ifelse(!is.na(sel_label), sel_label, Label)
                        ) %>% dplyr::select(-c(supertaxonID, sel_label))
                    }
                }
                # apply custom labels (if provided)
                customLabels <- getCustomLabels()
                if(nrow(customLabels) > 0) {
                    dimRedData$Label[
                        dimRedData$ncbiID %in% customLabels$ncbiID
                    ] <- customLabels$Label
                }
            }
        })
        return(dimRedData)
    })

    output$dimRedGroupHigherRank.warning <- renderUI({
        req(input$dimRedGroupHigherRank)
        if (length(input$dimRedGroupHigherRank) > 0) {
            list(
                em(paste("Click `Change labels` to apply. If you don't see your",
                         "specified labels, please check for typos!")),
                br(),br()
            )
        }
    })

    # * do DIM reduction -------------------------------------------------------
    dimRedCoord <- reactive({
        req(renameLabelsDimRed())
        input$plotDimRed
        params <- isolate(list(
            type = input$dimRedType,
            dataType = input$dimRedDataType,
            reductionTechnique = input$reductionTechnique,
            tsneIter = input$tsneIter,
            randomSeed = input$randomSeed,
            dim = ifelse(input$dimRedPlotType == "ggplot", 2, 3)
        ))
        set.seed(params$randomSeed)
        withProgress(
            message = "Performing dimension reduction...", value = 0.5, {
                dimRedData.coord <- dimReduction(renameLabelsDimRed(), 
                    params$type, params$dataType, params$randomSeed, 
                    params$reductionTechnique, params$dim, params$tsneIter
                )
                return(dimRedData.coord)
            }
        )
    })

    # * generate list of dimRed plot labels ------------------------------------
    output$dimRedTaxa.ui <- renderUI({
        req(renameLabelsDimRed())
        df <- groupLabelDimRedData(renameLabelsDimRed(), input$dimRedLabelNr)
        list(
            selectInput(
                "excludeDimRedTaxa", "Choose labels to hide", multiple = TRUE,
                c(levels(as.factor(df$Label))), selected = NULL
            ),
            selectInput(
                "highlightDimRedTaxa", "Choose labels to highlight",
                c(levels(as.factor(df$Label))), selected = NULL, multiple = TRUE
            )
        )
    })

    # * update dimRedLabelNr based on the freq of genes/taxa -------------------
    observe({
        req(dimRedData())
        df <- dimRedData()
        freqList <- sort(unique(dimRedData()$n))
        selectFreq <- tail(freqList, 5)[1]
        updateSliderInput(
            session, "dimRedLabelNr", "Freq cutoff", min = freqList[1],
            max = tail(freqList, 1), step = 1,
            value = c(selectFreq, tail(freqList, 1))
        )
    })

    # * create DIM reduction plot ----------------------------------------------
    observe({
        if (input$dimRedPlotType == "plotly") {
            updateSliderInput(
                session, "dimRedPlot.dotzoom", "Dot size zooming",
                min = 0, max = 100, step = 5, value = 0
            )
            updateSliderInput(
                session, "dimRedDotAlpha", "Transparent level", min = 0,
                max = 1, step = 0.05, value = 0
            )
        } else {
            updateSliderInput(
                session, "dimRedPlot.dotzoom", "Dot size zooming",
                min = -3, max = 10, step = 1, value = 0
            )
            updateSliderInput(
                session, "dimRedDotAlpha", "Transparent level", min = 0,
                max = 1, step = 0.05, value = 0.5
            )
        }
    })
    ranges <- reactiveValues(x = NULL, y = NULL)

    dimRedPlotData <- reactive({
        req(getMainInput())
        if(is.null(getMainInput())) stop("Input data is NULL!")
        req(dimRedCoord())
        req(renameLabelsDimRed())
        plotDf <- createDimRedPlotData(
            dimRedCoord(), renameLabelsDimRed(), 
            freqCutoff = input$dimRedLabelNr, 
            excludeTaxa = input$excludeDimRedTaxa,
            currentNCBIinfo = currentNCBIinfo
        )
        return(plotDf)
    })

    output$dimRedPlot <- renderPlot({
        req(getMainInput())
        if(is.null(getMainInput())) stop("Input data is NULL!")
        withProgress(
            message = "Plotting...", value = 0.5, {
                g <- plotDimRed(
                    dimRedPlotData(), legendPos = input$dimRedPlot.legend,
                    colorPalette = input$colorPalleteDimRed,
                    transparent = input$dimRedDotAlpha,
                    textSize = input$dimRedPlot.textsize, font = input$font,
                    highlightTaxa = input$highlightDimRedTaxa,
                    dotZoom = input$dimRedPlot.dotzoom
                )
                g + coord_cartesian(
                    xlim = ranges$x, ylim = ranges$y, expand = TRUE
                )
            }
        )
    })

    output$dimRedPlotly <- renderPlotly({
        req(getMainInput())
        options(htmlwidgets.TOJSON_ARGS = NULL)
        if(is.null(getMainInput())) stop("Input data is NULL!")
        withProgress(
            message = "Plotting...", value = 0.5, {
                g <- plotDimRed3D(
                    dimRedPlotData(), legendPos = input$dimRedPlot.legend,
                    colorPalette = input$colorPalleteDimRed,
                    transparent = input$dimRedDotAlpha,
                    highlightTaxa = input$highlightDimRedTaxa,
                    dotZoom = input$dimRedPlot.dotzoom
                )
                return(g)
            }
        )
    })

    output$dimRedPlot.ui <- renderUI({
        if (input$dimRedPlotType == "ggplot") {
            shinycssloaders::withSpinner(
                plotOutput(
                    "dimRedPlot",
                    height = input$dimRedPlot.height,
                    width = input$dimRedPlot.width,
                    click = "dimRedClick",
                    dblclick = "dimReddblClick",
                    brush = brushOpts(
                        id = "dimRedBrush",
                        delay = input$brushDelay,
                        delayType = input$brushPolicy,
                        direction = input$brushDir,
                        resetOnNew = TRUE
                    ),
                    hover = hoverOpts(
                        id = "dimRedHover",
                        nullOutside = FALSE
                    )
                )
            )
        } else {
            shinycssloaders::withSpinner(
                plotlyOutput(
                    "dimRedPlotly",
                    height = input$dimRedPlot.height,
                    width = input$dimRedPlot.width
                )
            )
        }
    })

    # When a double-click happens, check if there's a brush on the plot
    # If so, zoom to the brush bounds; if not, reset the zoom.
    observeEvent(input$dimReddblClick, {
        brush <- input$dimRedBrush
        if (!is.null(brush)) {
            ranges$x <- c(brush$xmin, brush$xmax)
            ranges$y <- c(brush$ymin, brush$ymax)

        } else {
            ranges$x <- NULL
            ranges$y <- NULL
        }
    })

    # * download DIM reduction plot & data -------------------------------------
    output$dimRedDownloadPlot <- downloadHandler(
        filename = function() {
            c("dimReduction.pdf")
        },
        content = function(file) {
            ggsave(
                file,
                plot = plotDimRed(
                    dimRedPlotData(), legendPos = input$dimRedPlot.legend,
                    colorPalette = input$colorPalleteDimRed,
                    transparent = input$dimRedDotAlpha,
                    textSize = input$dimRedPlot.textsize, font = input$font,
                    highlightTaxa = input$highlightDimRedTaxa,
                    dotZoom = input$dimRedPlot.dotzoom
                ) + coord_cartesian(
                    xlim = ranges$x, ylim = ranges$y, expand = TRUE),
                width = input$dimRedPlot.width * 0.056458333,
                height = input$dimRedPlot.height * 0.056458333,
                units = "cm", dpi = 300, device = "svg", limitsize = FALSE
            )
        }
    )

    output$dimRedDownloadData <- downloadHandler(
        filename = function() {
            c("dimRedData.RData")
        },
        content = function(fileName) {
            data4dimRed <- renameLabelsDimRed()
            dimRedCoord <- dimRedCoord()
            dimRedPlotData <- dimRedPlotData()
            save(data4dimRed, dimRedCoord, dimRedPlotData, file = fileName)
        }
    )

    # * create brushed DIM reduction table -------------------------------------
    brushedDimRedData <- reactive({
        # get list of selected gene(s)
        if (is.null(input$dimRedBrush))
            return()
        else {
            xmin <- input$dimRedBrush$xmin
            xmax <- input$dimRedBrush$xmax
            ymin <- input$dimRedBrush$ymin
            ymax <- input$dimRedBrush$ymax
            return(
                dimRedPlotData() %>%
                    filter(X >= xmin & X <= xmax & Y >= ymin & Y <= ymax)
            )
        }
    })

    output$dimRedTable.ui <- renderUI({
        if (input$dimRedType == "taxa") {
            list(
                hr(),
                h4("SELECTED TAXA"),
                DT::dataTableOutput("dimRedSpec.table"),
                downloadButton(
                    "downloadDimRedSpec.table", "Download table",
                    class = "butDL"
                ),
                hr(),
                h4("SELECTED SEED GENES"),
                DT::dataTableOutput("dimRedGenes.table"),
                downloadButton(
                    "downloadDimRedGenes.table", "Download table",
                    class = "butDL"
                )
            )
        } else {
            list(
                hr(),
                h4("SELECTED SEED GENES"),
                DT::dataTableOutput("dimRedGenes.table"),
                downloadButton(
                    "downloadDimRedGenes.table", "Download table",
                    class = "butDL"
                )
            )
        }
    })

    # ** DIM reduction selected taxa table -------------------------------------
    dimRedSelectedTaxa <- reactive({
        if (is.null(input$dimRedBrush$ymin)) {
            shinyjs::disable("addSpecDimRed")
            return()
        } else {
            lowestRank <- getLowestRank(getMainInput(), getTaxDBpath())
            # enable add taxa to cus. profile only when lowest rank is selected
            if (lowestRank == input$rankSelect) {
                shinyjs::enable("addSpecDimRed")
            } else shinyjs::disable("addSpecDimRed")
        }
        df <- as.data.frame(brushedDimRedData())
        if (nrow(df) > 0) {
            specDf <- df %>% dplyr::select(ncbiID, fullName, Label, Freq, n)
            colnames(specDf) <- c(
                "NCBI taxon ID", "Taxon name", "Label", "Number of genes", "Freq"
            )
            return(specDf)
        } else {
            shinyjs::disable("addSpecDimRed")
            return()
        }
    })

    output$dimRedSpec.table <- DT::renderDataTable(
        options = list(searching = TRUE, pageLength = 10
    ),{
        dimRedSelectedTaxa()
    })

    output$downloadDimRedSpec.table <- downloadHandler(
        filename = function() {
            c("dimRedSelectedTaxa.txt")
        },
        content = function(file) {
            dataOut <- dimRedSelectedTaxa()
            write.table(
                dataOut, file, sep = "\t", row.names = FALSE, quote = FALSE
            )
        }
    )

    # ** DIM reduction selected genes table ------------------------------------
    dimRedSelectedGenes <- reactive({
        if (is.null(input$dimRedBrush$ymin)){
            shinyjs::disable("addGeneDimRed")
            return()
        } else {
            if (input$addClusterCustomProfile == FALSE
                & input$addGeneAgeCustomProfile == FALSE
                & input$addCoreGeneCustomProfile == FALSE
                & input$addGCGenesCustomProfile == FALSE) {
                shinyjs::enable("addGeneDimRed")
            } else {
                shinyjs::disable("addGeneDimRed")
            }
        }
        df <- as.data.frame(brushedDimRedData())
        if (nrow(df) > 0) {
            removeDf <- df %>% dplyr::select(where(~ all(. == -1)))
            subDf <- df %>% dplyr::select(-c(colnames(removeDf), Label, Freq, X, Y, Z, n))
            if ("fullName" %in% colnames(subDf))
                subDf <- subDf %>% dplyr::select(-c("fullName"))
            if ("ncbiID" %in% colnames(df)) {
                meltedDf <- data.frame(melt(
                    as.data.table(subDf),
                    id.vars = "ncbiID", variable.name = "geneID"
                ))
            } else {
                meltedDf <- data.frame(melt(
                    as.data.table(subDf),
                    id.vars = "geneID", variable.name = "ncbiID"
                ))
            }
            geneDf <- meltedDf %>% filter(value >= 0) %>%
                select(geneID, ncbiID, value)
            geneDf$value <- round(geneDf$value,2)
            # get total number of taxa
            totalTaxa <- length(unique(geneDf$ncbiID))
            # rename "value" based on variable name(s)
            longDf <- getMainInput()
            if (ncol(longDf) >= 5) {
                colnames(geneDf)[colnames(geneDf) == "value"] <-
                    paste("Mean of", input$var1ID, "&", input$var2ID)
            } else if (ncol(longDf) == 4) {
                colnames(geneDf)[colnames(geneDf) == "value"] <- input$var1ID
            } else if (ncol(longDf) == 3) {
                geneDf <- geneDf %>% dplyr::select(geneID, ncbiID)
            }
            # check gene IDs containing character "X" at the beginning
            geneDf$geneID <- as.character(geneDf$geneID)
            if (all(!(unique(geneDf$geneID) %in% longDf$geneID))) {
                if (all(startsWith(unique(geneDf$geneID), "X"))){
                    geneDf$geneID <- sub("X","", geneDf$geneID, fixed = TRUE)
                }
            }
            # count number of taxa for each genes
            joinedIDs <- paste(geneDf$geneID, geneDf$ncbiID)
            t <- as.data.frame(do.call(rbind, strsplit(unique(joinedIDs), " ")))
            geneCountDf <- data.frame(table(t$V1))
            geneCountDf$Percentage <- round((geneCountDf$Freq / totalTaxa)*100, 1)
            colnames(geneCountDf) <- c("geneID", "Taxa count", "% taxa")
            geneDf <- merge(geneDf, geneCountDf, by = "geneID", all.x = TRUE)
            # return(geneDf)
            return(geneCountDf)
        } else {
            shinyjs::disable("addGeneDimRed")
            return()
        }
    })

    output$dimRedGenes.table <- DT::renderDataTable(
        options = list(searching = TRUE, pageLength = 10
    ),{
        dimRedSelectedGenes()
    })

    output$downloadDimRedGenes.table <- downloadHandler(
        filename = function() {
            c("dimRedSelectedGenes.txt")
        },
        content = function(file) {
            dataOut <- dimRedSelectedGenes()
            write.table(dataOut, file, sep = "\t", row.names = FALSE,
                        quote = FALSE)
        }
    )

    # * show DIM reduction hover info ------------------------------------------
    dimRedHoverTaxa <- reactive({
        req(input$dimRedHover)
        x <- input$dimRedHover$x
        y <- input$dimRedHover$y
        hoveredDf <- nearPoints(
            dimRedPlotData(), input$dimRedHover, xvar = "X", yvar = "Y"
        ) %>% dplyr::select(ncbiID, Label, Freq, fullName)
        if (nrow(hoveredDf) > 0) {
            hoveredDf <- hoveredDf %>% dplyr::select(fullName, Freq)
            colnames(hoveredDf) <- c("Taxon", "# genes")
            return(hoveredDf)
        } else return()
    })

    dimRedHoverGene <- reactive({
        req(input$dimRedHover)
        x <- input$dimRedHover$x
        y <- input$dimRedHover$y
        hoveredDf <- nearPoints(
            dimRedPlotData(), input$dimRedHover, xvar = "X", yvar = "Y"
        ) %>% dplyr::select(geneID, Freq)
        if (nrow(hoveredDf) > 0) {
            colnames(hoveredDf) <- c("Gene", "# taxa")
            return(hoveredDf)
        } else return()
    })

    output$dimRedHoverInfo <- renderTable({
        if(!is.null(input$dimRedHover)){
            if (input$dimRedType == "taxa") {
                dimRedHoverTaxa()
            } else dimRedHoverGene()
        }
    })

    # * check if genes are added anywhere else to the customized profile -------
    observe({
        if (input$addClusterCustomProfile == TRUE
            | input$addGeneAgeCustomProfile == TRUE
            | input$addCoreGeneCustomProfile == TRUE
            | input$addGCGenesCustomProfile == TRUE) {
            shinyjs::disable("addGeneDimRed")
        } else {
            shinyjs::enable("addGeneDimRed")
        }
    })

    output$addDimRedCustomProfileCheck.ui <- renderUI({
        if (input$addClusterCustomProfile == TRUE
            | input$addGeneAgeCustomProfile == TRUE
            | input$addCoreGeneCustomProfile == TRUE |
            input$addGCGenesCustomProfile == TRUE ) {
            HTML('<p><em>(Uncheck "Add to Customized profile" check box in
                 <strong>Profile clustering</strong> or
                 <strong>Gene age estimation</strong> or
                 <strong>Core genes finding</strong> or
                 <strong>Group comparison</strong>
                 &nbsp;to enable this function)</em></p>')
        }
    })

    # ============================== POINT INFO ================================

    # * get status of pointInfo for activating Detailed Plot button -----------
    output$pointInfoStatus <- reactive({
        if (input$tabs == "Main profile") {
            # info contains groupID,orthoID,supertaxon,mVar1,%spec,var2
            info <- mainpointInfo()
        } else if (input$tabs == "Customized profile") {
            info <- selectedpointInfo()
        } else info <- NULL
        return(is.null(info))
    })
    outputOptions(output, "pointInfoStatus", suspendWhenHidden = FALSE)

    # * show info into "point's info" box --------------------------------------
    output$pointInfo <- renderText({
        # GET INFO BASED ON CURRENT TAB
        if (input$tabs == "Main profile") {
            # info contains groupID,orthoID,supertaxon,mVar1,%spec,var2
            info <- mainpointInfo()
        } else if (input$tabs == "Customized profile") {
            info <- selectedpointInfo()
        } else return()

        req(info)
        orthoID <- info[[2]]
        if (length(info[[2]]) > 1) orthoID <- paste0(info[[2]][1], ",...")

        if (is.na(orthoID)) return()
        else {
            a <- toString(paste("Seed-ID:", info[[1]]))
            b <- toString(paste0(
                "Hit-ID: ", orthoID
            ))
            c <- ""
            if (input$var1ID != "") {
                c <- toString(paste(
                    input$var1AggregateBy, input$var1ID, ":", info[[5]]
                ))
            }
            d <- ""
            if (input$var2ID != "") {
                d <- toString(paste(
                    input$var2AggregateBy, input$var2ID, ":", info[[7]]
                ))
            }

            if (info[[10]] == "Y") {
                if (info[[3]] == 1) {
                    s <- toString(
                        paste0(info[[3]], " ortholog in ", info[[4]], ":")
                    )
                } else {
                    s <- toString(
                        paste0(info[[3]], " co-orthologs in ", info[[4]], ":")
                    )
                }
                e <- ""
            } else {
                s <- toString(paste0("Best ortholog in ", info[[4]], ":"))
                e <- toString(
                    paste0(
                        "% present taxa: ", info[[6]], " (", info[[8]], " out of ",
                        info[[9]],  ")", collapse = ""
                    )
                )
            }
            paste(a, s, b, c, d, e, sep = "\n")
        }
    })
    
    # * update highlight gene/taxa by selected point ---------------------------
    observe({
        input$highlightMain
        # GET INFO BASED ON CURRENT TAB
        if (input$tabs == "Main profile") {
            shinyjs::enable("highlightMain")
            # info contains groupID,orthoID,supertaxon,mVar1,%spec,var2
            info <- isolate(mainpointInfo())
        } else if (input$tabs == "Customized profile") {
            shinyjs::disable("highlightMain")
            info <- NULL
        }
        req(info)
        if (length(info) > 1) {
            updateSelectizeInput(
                session, "geneHighlight", selected = info[[1]][1]
            )
            updateSelectizeInput(
                session, "taxonHighlight", selected = info[[4]][1]
            )
        }
    })

    # ============================= DETAILED PLOT ==============================
    # * data for detailed plot -------------------------------------------------
    detailPlotDt <- reactive({
        req(v$doPlot)
        req(input$detailedBtn)
        # GET INFO BASED ON CURRENT TAB
        if (input$tabs == "Main profile") {
            # info contains groupID,orthoID,supertaxon,mVar1,%spec,var2
            info <- mainpointInfo()
        } else if (input$tabs == "Customized profile") {
            info <- selectedpointInfo()
        }
        req(info)
        withProgress(message = 'Getting data for detailed plot...', value=0.5, {
            ### get refspec name
            split <- strsplit(as.character(input$inSelect), "_")
            inSelect <- as.character(split[[1]][1])

            ### get info for present taxa in selected supertaxon (1)
            fullDf <- getFullData()
            ### filter data if needed
            if  (input$detailedFilter == TRUE) {
                fullDf <- filteredDataHeat()
                if (info[[4]] == inSelect) {
                    fullDf <- fullDf[
                        fullDf$var1 >= input$var1[1]
                        & fullDf$var1 <= input$var1[2],
                    ]
                    fullDf <- fullDf[
                        fullDf$var2 >= input$var2[1]
                        & fullDf$var2 <= input$var2[2],
                    ]
                }
                updateCheckboxInput(
                    session, "detailedRemoveNA", value = TRUE
                )
            }
            selTaxon <- info[[4]]
            selTaxon <- gsub("\\[","\\\\[",selTaxon)
            selTaxon <- gsub("\\]","\\\\]",selTaxon)
            plotTaxon <- unique(
                fullDf$supertaxon[grep(paste0("_",selTaxon,"$"), fullDf$supertaxon)]
            )
            plotGeneID <- info[[1]]
            selDf <- fullDf[fullDf$geneID == plotGeneID
                            & fullDf$supertaxon == plotTaxon, ]
            ### get all taxa of this supertaxon (2)
            allTaxaDf <- sortedtaxaList()
            allTaxaDf <- allTaxaDf[allTaxaDf$supertaxon == plotTaxon,
                                   c("abbrName", "fullName")]

            ### merge (1) and (2) together
            joinedDf <- merge(selDf, allTaxaDf, by = c("abbrName"), all.y =TRUE)
            joinedDf <- subset(
                joinedDf,
                select = c(
                    "abbrName", "fullName.y", "geneID", "orthoID", "var1","var2"
                )
            )
            names(joinedDf)[names(joinedDf) == "fullName.y"] <- "fullName"

            # replace var1/var2 as NA for all "NA orthologs"
            joinedDf$var1[is.na(joinedDf$orthoID)] <- NA
            joinedDf$var2[is.na(joinedDf$orthoID)] <- NA

            # remove NA orthologs if required
            if (input$detailedRemoveNA == TRUE) {
                joinedDf <- joinedDf[!is.na(joinedDf$orthoID), ]
            }

            ### return data for detailed plot
            return(joinedDf)
        })
    })

    # * render detailed plot ---------------------------------------------------
    pointInfoDetail <- callModule(
        createDetailedPlot, "detailedPlot",
        data = detailPlotDt,
        var1ID = reactive(input$var1ID),
        var2ID = reactive(input$var2ID),
        detailedText = reactive(input$detailedText),
        detailedHeight = reactive(input$detailedHeight),
        font = reactive(input$font)
    )

    # * render database links --------------------------------------------------
    parseProId <- function(protId, separator, seqIdFormat) {
        tmp <- ""
        if (separator == 1) {
            if (grepl("\\|", protId))
                tmp <- as.list(strsplit(protId, "\\|")[[1]])
        } else if (separator == 2) {
            if (grepl("@", protId))
                tmp <- as.list(strsplit(protId, "@")[[1]])
        } else if (separator == 3) {
            if (grepl("#", protId))
                tmp <- as.list(strsplit(protId, "#")[[1]])
        } else if (separator == 4) {
            if (grepl(";", protId))
                tmp <- as.list(strsplit(protId, ";")[[1]])
        }
        if (length(tmp) > 1) {
            if (seqIdFormat == 3) {
                protId <- tmp[[length(tmp)]]
            } else if (seqIdFormat == 4) {
                protId <- tmp[[1]]
            } else if (seqIdFormat == 1) {
                if (length(tmp) >= 3) {
                    protId <- tmp[[3]]
                } else {
                    warning("Wrong ID format was set!")
                    return(NULL)
                }
            }
        }
        return(protId)
    }

    output$dbLink.ui <- renderUI({
        info <- pointInfoDetail() # info = seedID, orthoID, var1, var2, ncbiID
        req(info)
        linkText <- ""
        # get seed ID
        seedId <- toString(info[1])
        separators <- c("|", "@", "#", ";")
        if (grepl(paste(separators, collapse = "|"), seedId)) {
            seedId <- parseProId(seedId, input$separator, input$seqIdFormat)
        }
        if (length(seedId) > 0) {
            if (input$seedSource == "ncbi") {
                linkText <- paste0(linkText, createDBlink(seedId, "NCBI"))
            } else if (input$seedSource == "uniprot") {
                linkText <- paste0(linkText, createDBlink(seedId, "UniProt"))
            } else if (input$seedSource == "orthodb") {
                linkText <- paste0(linkText, createDBlink(seedId, "OrthoDB", input$orthodbSeedVer))
            } else if (input$seedSource == "oma") {
                linkText <- paste0(linkText, createDBlink(seedId, "OMA"))
            }
        }
        # get ortho ID
        protId <- toString(info[2])
        if (grepl(paste(separators, collapse = "|"), protId)) {
            protId <- parseProId(protId, input$separator, input$seqIdFormat)
        }
        if (length(protId) > 0) {
            if (input$orthoSource == "ncbi") {
                linkText <- paste0(linkText, createDBlink(protId, "NCBI"))
            } else if (input$orthoSource == "uniprot") {
                linkText <- paste0(linkText, createDBlink(protId, "UniProt"))
            } else if (input$orthoSource == "orthodb") {
                linkText <- paste0(linkText,createDBlink(protId, "OrthoDB", "gene", input$orthodbOrthoVer))
            } else if (input$orthoSource == "oma") {
                linkText <- paste0(linkText, createDBlink(protId, "OMA", "gene"))
            }
        }
        # get taxon ID
        taxId <- gsub("ncbi", "", info[5])
        taxHierarchy <- PhyloProfile:::getTaxHierarchy(taxId, currentNCBIinfo)
        taxUrls <- paste(taxHierarchy[[1]]$link, collapse = ", ")
        taxUrls <- gsub("<p>", "", taxUrls)
        taxUrls <- gsub("</p>", "", taxUrls)
        linkText <- paste0(
            linkText, "<p><strong>NCBI taxonomy: </strong>", taxUrls, "</p>"
        )

        # render links
        linkText <- paste0(
            linkText,
            "<p><em><strong>Disclaimer:</strong> ",
            "External links are automatically generated and may point to ",
            "a wrong target. Please adapt the sequence ID format according ",
            "to your data in the Input and Settings tab (see <a ",
            "href=\"https://github.com/BIONF/PhyloProfile/wiki/FAQ",
            "#wrong-info-from-public-databases\" ",
            "target=\"_blank\">FAQ</a>)</em></p>"
        )
        HTML(linkText)
    })

    # * render FASTA sequence --------------------------------------------------
    output$fasta <- renderText({
        req(v$doPlot)
        info <- pointInfoDetail() # info = seedID, orthoID, var1

        req(info)
        seqID <- toString(info[2])

        if (input$demoData == "none") {
            filein <- input$mainInput
            inputType <- inputType()
            # get fata from oma
            if (inputType == "oma") {
                fastaOut <- getSelectedFastaOma(finalOmaDf(), seqID)
            }
            # get fasta from main input
            else if (inputType == "fasta") {
                seqIDMod <- paste0(info[1], "|", info[5], "|", info[2])
                fastaOut <- getFastaFromFasInput(
                    seqIDMod, file = filein$datapath
                )
            } else {
                # get from concaternated file
                if (input$inputType == "Concatenated fasta file") {
                    fastain <- input$concatFasta
                    fastaOut <- getFastaFromFile(seqID, fastain$datapath)
                }
                # get from folder
                else {
                    fastaOut <- getFastaFromFolder(
                        seqID,
                        input$path,
                        input$dirFormat,
                        input$fileExt,
                        input$idFormat
                    )
                }
            }
        } else {
            if (input$demoData == "preCalcDt") {
                if (!is.null(i_fastaInput))
                    if (!file.exists(i_fastaInput))
                        return("Fasta file not found!")
                    fastaOut <- getFastaFromFile(seqID, i_fastaInput)
            } else {
                # get fasta from demo online data
                fastaOut <- getFastaDemo(seqID, demoData = input$demoData)
            }
        }
        return(paste(fastaOut[1]))
    })

    # ======================== FEATURE ARCHITECTURE PLOT =======================
    # * check if seed and orthoID are specified and get domain file/path -------
    observe({
        infoTmp <- c()
        if (input$tabs == "Main profile") {
            infoTmp <- mainpointInfo()
        } else if (input$tabs == "Customized profile") {
            infoTmp <-selectedpointInfo()
        }
        if (!is.null(infoTmp)) {
            tmp <- getDomainFile()
        }
    })

    getDomainFile <- reactive({
        # get lowest rank
        # activate doDomainPlotMain if either working on the lowest rank
        # or only 1 ortholog present
        longDataframe <- getMainInput()
        req(longDataframe)
        req(input$rankSelect)
        lowestRank <- getLowestRank(longDataframe, getTaxDBpath())

        # get info from POINT INFO box
        info <- c()
        infoTmp <- c()
        if (input$tabs == "Main profile") {
            # info contains groupID,orthoID,supertaxon,mVar1,%spec,var2
            infoTmp <- mainpointInfo()
        } else if (input$tabs == "Customized profile") {
            infoTmp <- selectedpointInfo()
        }
        # only when no co-ortholog exists
        if (length(infoTmp[[2]]) == 1) {
            info <- c(infoTmp[[1]], as.character(infoTmp[[2]]))
        } else {
            # else, get info from detailed plot
            shinyBS::updateButton(session, "doDomainPlotMain", disabled = TRUE)
            if (!is.null(pointInfoDetail())) {
                info <- pointInfoDetail() # info = seedID, orthoID, var1
            }
        }

        if (is.null(info)) {
            shinyBS::updateButton(session, "doDomainPlot", disabled = TRUE)
            shinyBS::updateButton(session, "doDomainPlotMain", disabled = TRUE)
            return("noSelectHit")
        } else {
            if (
                input$demoData == "arthropoda" | input$demoData == "ampk-tor" |
                input$demoData == "preCalcDt"
            ) {
                shinyBS::updateButton(session, "doDomainPlot", disabled = FALSE)
                if (lowestRank == input$rankSelect || infoTmp[[8]][1] == 1)
                    shinyBS::updateButton(session, "doDomainPlotMain", disabled = FALSE)
                else
                    shinyBS::updateButton(session, "doDomainPlotMain", disabled = TRUE)
            } else {
                if (
                    input$mainInputType == "file" &
                    inputType() == "oma"
                ) {
                    shinyBS::updateButton(session, "doDomainPlot", disabled = FALSE)
                    if (lowestRank == input$rankSelect || infoTmp[[8]][1] == 1)
                        shinyBS::updateButton(session, "doDomainPlotMain", disabled = FALSE)
                    else
                        shinyBS::updateButton(session, "doDomainPlotMain", disabled = TRUE)
                } else {
                    if (input$annoLocation == "from file") {
                        inputDomain <- input$fileDomainInput
                        if (is.null(inputDomain)) {
                            shinyBS::updateButton(
                                session, "doDomainPlot", disabled = TRUE
                            )
                            shinyBS::updateButton(
                                session, "doDomainPlotMain", disabled = TRUE
                            )
                            return("noFileInput")
                        } else {
                            shinyBS::updateButton(
                                session, "doDomainPlot", disabled = FALSE
                            )
                            if (lowestRank == input$rankSelect || infoTmp[[8]][1] == 1)
                                shinyBS::updateButton(session, "doDomainPlotMain", disabled = FALSE)
                            else
                                shinyBS::updateButton(session, "doDomainPlotMain", disabled = TRUE)
                        }
                    } else {
                        domainDf <- parseDomainInput(
                            info[1], input$domainPath, "folder"
                        )
                        if (length(domainDf) == 1) {
                            if (domainDf == "noSelectHit" |
                                domainDf == "noFileInFolder") {
                                shinyBS::updateButton(
                                    session, "doDomainPlot", disabled = TRUE
                                )
                                shinyBS::updateButton(
                                    session, "doDomainPlotMain", disabled = TRUE
                                )
                                return(domainDf)
                            } else {
                                shinyBS::updateButton(
                                    session, "doDomainPlot", disabled = FALSE
                                )
                                if (lowestRank == input$rankSelect || infoTmp[[8]][1] == 1)
                                    shinyBS::updateButton(session, "doDomainPlotMain", disabled = FALSE)
                                else
                                    shinyBS::updateButton(session, "doDomainPlotMain", disabled = TRUE)
                            }
                        } else {
                            shinyBS::updateButton(
                                session, "doDomainPlot", disabled = FALSE
                            )
                            if (lowestRank == input$rankSelect || infoTmp[[8]][1] == 1)
                                shinyBS::updateButton(session, "doDomainPlotMain", disabled = FALSE)
                            else
                                shinyBS::updateButton(session, "doDomainPlotMain", disabled = TRUE)
                        }
                    }
                }
            }
        }
        return(info)
    })

    # * check domain file ------------------------------------------------------
    output$checkDomainFiles <- renderUI({
        fileDomain <- getDomainFile()
        if (length(fileDomain) == 1) {
            if (fileDomain == "noFileInput") {
                em("Domain file not provided!!")
            } else if (fileDomain == "noFileInFolder") {
                msg <- paste0(
                    "<p><em>Domain file not found!! </em></p>
                <p><em>Please make sure that file name has to be in this format:
                <strong>&lt;seedID&gt;.extension</strong>, where extension is
                limited to <strong>txt</strong>, <strong>csv</strong>,
                <strong>list</strong>, <strong>domains</strong> or
                <strong>architecture</strong>.</em></p>"
                )
                HTML(msg)
            } else if (fileDomain == "noSelectHit") {
                em("Please select one ortholog sequence!!")
            }
        }
    })

    # * render domain plot -----------------------------------------------------
    observeEvent(input$doDomainPlot, {
        callModule(
            createArchitecturePlot, "archiPlot",
            pointInfo = getDomainFile,
            domainInfo = getDomainInformation,
            currentNCBIinfo = reactive(currentNCBIinfo),
            font = reactive(input$font)
        )
    })
    observeEvent(input$doDomainPlotMain, {
        callModule(
            createArchitecturePlot, "archiPlotMain",
            pointInfo = getDomainFile,
            domainInfo = getDomainInformation,
            currentNCBIinfo = reactive(currentNCBIinfo),
            font = reactive(input$font)
        )
    })

    # ======================== FILTERED DATA DOWNLOADING =======================

    # * for main profile =======================================================
    mainFastaDownload <- reactive({
        downloadDf <- as.data.frame(downloadData())
        seqIDs <- downloadDf$orthoID

        if (input$demoData == "none") {
            filein <- input$mainInput
            inputType <- inputType()
            # get fata from oma
            if (inputType == "oma") {
                allOmaDf <- finalOmaDf()
                filteredDownloadDf <- as.data.frame(downloadData())
                filteredOmaDf <-
                    subset(allOmaDf,
                           allOmaDf$orthoID %in% filteredDownloadDf$orthoID &
                               allOmaDf$seed %in% filteredDownloadDf$geneID)
                mainFastaOut <- getAllFastaOma(filteredOmaDf)
            }
            # get fasta from main input
            else if (inputType == "fasta") {
                seqIDMod <- paste0(
                    as.character(downloadDf$geneID), "|",
                    as.character(downloadDf$ncbiID), "|",
                    as.character(downloadDf$orthoID)
                )
                mainFastaOut <- getFastaFromFasInput(
                    seqIDMod, file = filein$datapath
                )
            } else {
                # get from concaternated file
                if (input$inputType == "Concatenated fasta file") {
                    fastain <- input$concatFasta
                    mainFastaOut <- getFastaFromFile(seqIDs, fastain$datapath)
                }
                # get from folder
                else {
                    mainFastaOut <- getFastaFromFolder(
                        seqIDs,
                        input$path,
                        input$dirFormat,
                        input$fileExt,
                        input$idFormat
                    )
                }
            }
        } else {
            if (input$demoData == "preCalcDt") {
                if (!is.null(i_fastaInput))
                    mainFastaOut <- getFastaFromFile(seqIDs, i_fastaInput)
            } else {
                # get fasta from demo online data
                mainFastaOut <- getFastaDemo(seqIDs, demoData = input$demoData)
            }
        }
        return(mainFastaOut)
    })

    downloadData <- callModule(
        downloadFilteredMain,
        "filteredMainDownload",
        data = getFullData,
        taxaCount = getCountTaxa,
        fasta = mainFastaDownload,
        var1ID = reactive(input$var1ID),
        var2ID = reactive(input$var2ID),
        var1 = reactive(input$var1),
        var2 = reactive(input$var2),
        percent = reactive(input$percent)
    )

    # * for customized profile =================================================
    customizedFastaDownload <- reactive({
        downloadDf <- as.data.frame(downloadCustomData())
        seqIDs <- downloadDf$orthoID

        if (input$demoData == "none") {
            filein <- input$mainInput
            inputType <- inputType()
            # get fata from oma
            if (inputType == "oma") {
                allOmaDf <- finalOmaDf()
                filteredDownloadDf <- as.data.frame(downloadCustomData())
                filteredOmaDf <-
                    subset(allOmaDf,
                           allOmaDf$orthoID %in% filteredDownloadDf$orthoID &
                               allOmaDf$seed %in% filteredDownloadDf$geneID)
                fastaOutDf <- getAllFastaOma(filteredOmaDf)
            }
            # get fasta from main input
            else if (inputType == "fasta") {
                seqIDMod <- paste0(
                    as.character(downloadDf$geneID), "|",
                    as.character(downloadDf$ncbiID), "|",
                    as.character(downloadDf$orthoID)
                )
                fastaOutDf <- getFastaFromFasInput(
                    seqIDMod, file = filein$datapath
                )
            } else {
                # get from concaternated file
                if (input$inputType == "Concatenated fasta file") {
                    fastain <- input$concatFasta
                    fastaOutDf <- getFastaFromFile(seqIDs, fastain$datapath)
                }
                # get from folder
                else {
                    fastaOutDf <- getFastaFromFolder(
                        seqIDs,
                        input$path,
                        input$dirFormat,
                        input$fileExt,
                        input$idFormat
                    )
                }
            }
        } else {
            if (input$demoData == "preCalcDt") {
                if (!is.null(i_fastaInput))
                    fastaOutDf <- getFastaFromFile(seqIDs, i_fastaInput)
            } else {
                # get fasta from demo online data
                fastaOutDf <- getFastaDemo(seqIDs, demoData = input$demoData)
            }
        }
        return(fastaOutDf)
    })

    downloadCustomData <- callModule(
        downloadFilteredCustomized,
        "filteredCustomizedDownload",
        data = downloadData,
        fasta = customizedFastaDownload,
        inSeq = reactive(input$inSeq),
        inTaxa = reactive(input$inTaxa)
    )

    # ======================== PROCESSED DATA DOWNLOADING ======================
    # * description for plot settings downloading ------------------------------
    observe({
        desc = paste(
            "Here you can download the processed data in RDS format.",
            "These data allows for faster analysis, particularly useful when",
            "dealing with a large number of genes and taxa. To reuse the data,",
            "simply upload a folder containing these RDS files."
        )

        if (input$tabs == "Processed data") {
            shinyBS::createAlert(
                session, "descDownloadProcessedDataUI",
                "descDownloadProcessedData", content = desc, append = FALSE
            )
        }
    })

    # * get output path --------------------------------------------------------
    getProcDataPath <- reactive({
        shinyFiles::shinyDirChoose(
            input, "procDataOutDir", roots = homePath, session = session
        )
        settingPath <- shinyFiles::parseDirPath(homePath, input$procDataOutDir)
        return(replaceHomeCharacter(as.character(settingPath)))
    })

    output$procDataOutDir.ui <- renderUI({
        req(getProcDataPath())
        em(getProcDataPath())
    })

    # * do export processed data -----------------------------------------------
    observe({
        if (v$doPlot) shinyjs::enable("doDownloadProcData")
        else shinyjs::disable("doDownloadProcData")
    })

    downloadProcData <- function() {
        filein <- input$mainInput
        outDirName <- filein$name
        outDirName <- gsub("[^[:alnum:] ]", "_", outDirName)
        outFolder <- paste0(getProcDataPath(), "/", outDirName)
        if (dir.exists(outFolder)) {
            if (
                file.exists(paste0(outFolder, "/longDf.rds")) |
                file.exists(paste0(outFolder, "/sortedtaxaList.rds")) |
                file.exists(paste0(outFolder, "/preData.rds")) |
                file.exists(paste0(outFolder, "/fullData.rds"))
            ) {
                message("Other RDS files found! Please select a new output folder!")
                return()
            }
        } else dir.create(file.path(outFolder), showWarnings = TRUE)

        if (nrow(getFullData()) > 0) {
            saveRDS(getMainInput(), file = paste0(outFolder, "/longDf.rds"))
            saveRDS(sortedtaxaList(), file = paste0(outFolder, "/sortedtaxaList.rds"))
            saveRDS(preData(), file = paste0(outFolder, "/preData.rds"))
            saveRDS(getFullData(), file = paste0(outFolder, "/fullData.rds"))
            msg <- paste("Done! Data have been saved at", outFolder)
            message(msg)
            shinyjs::disable("doDownloadProcData")
        }
    }

    observeEvent(input$doDownloadProcData, {
        req(length(getProcDataPath()) > 0)
        withCallingHandlers({
            shinyjs::html("downloadProcDataStatus", "<p>Please wait...</p>")
            downloadProcData()
        },
        message = function(m) {
            shinyjs::html(
                id = "downloadProcDataStatus", html = m$message, add = TRUE
            )
        })
    })

    # ======================== PLOT SETTINGS DOWNLOADING =======================
    # * description for plot settings downloading -----------------------------
    observe({
        desc = paste(
            "Here you can download the plot settings (such as plot size, font",
            "size, filters, working taxonomy rank and the reference taxon,",
            "etc.) either as <em><strong>a list in YAML",
            "format</strong></em>, or as <em><strong>an Rscript</strong></em>",
            "that you can directly use to create the corresponding plot",
            "without the need of using the GUI of",
            "PhyloProfile."
        )

        if (input$tabs == "Export plot settings") {
            shinyBS::createAlert(
                session, "descExportSettingUI", "descExportSetting",
                content = desc, append = FALSE
            )
        }
    })

    # * render output file name -----------------------------------------------
    output$settingFile.ui <- renderUI({
        if (input$exportSetting == "list") {
            textInput(
                "settingFile",
                "",
                value = "plotSettings",
                width = "30%",
                placeholder = "Output file name"
            )
        } else {
            textInput(
                "settingFile",
                "",
                value = "preConfiguredPlot",
                width = "30%",
                placeholder = "Output file name"
            )
        }

    })

    # * do export plot settings -----------------------------------------------
    getPlotSettings <- function(){
        rawInput <- ""
        if (input$demoData == "preCalcDt") {
            rawInput <- createLongMatrix(i_mainInput)
        } else {
            rawInput <- input$mainInput
        }
        outputList <- list(
            "mainInput" = rawInput$datapath,
            "rank" = input$rankSelect,
            "refspec" = input$inSelect,
            "clusterProfile" = input$applyCluster,
            "profileTypeClustering" = input$profileType,
            "distMethodClustering" = input$distMethod,
            "clusterMethod" = input$clusterMethod,
            "taxonHighlight" = input$taxonHighlight,
            "geneHighlight" = input$geneHighlight,
            "var1AggregateBy" = input$var1AggregateBy,
            "var2AggregateBy" = input$var2AggregateBy,
            "percentCutoff" = input$percent,
            "var1Cutoff" = input$var1,
            "var2Cutoff" = input$var2,
            "var1Relation" = input$var1Relation,
            "var2Relation" = input$var2Relation,
            "taxDB" = getTaxDBpath()
        )
        outputList <- append(outputList, getParameterInputMain())
        return(outputList)
    }

    getSettingPath <- reactive({
        shinyFiles::shinyDirChoose(
            input, "settingDir", roots = homePath, session = session
        )
        settingPath <- shinyFiles::parseDirPath(homePath, input$settingDir)
        settingPathMod <- replaceHomeCharacter(as.character(settingPath))
        if (length(settingPathMod) > 0) {
            if (input$exportSetting == "list")
                settingPathMod <- paste0(
                    settingPathMod, "/", input$settingFile, ".yml"
                )
            else
                settingPathMod <- paste0(
                    settingPathMod, "/", input$settingFile #, ".rscript"
                )
        }
        return(settingPathMod)
    })

    output$settingDir.ui <- renderUI({
        req(getSettingPath())
        if (length(getSettingPath()) > 0) {
            outString <- getSettingPath()
            if (input$exportSetting == "list") em(outString)
            else
                em(
                    paste0(
                        paste0(outString, ".yml"), "; ",
                        paste0(outString, ".rscript")
                    )
                )
        }
    })

    observeEvent(input$doExportSetting, {
        req(length(getSettingPath()) > 0)
        fileCheck <- 0
        withCallingHandlers({
            shinyjs::html("exportSettingStatus", "")
            downloadPlotSettings(
                getPlotSettings(), getSettingPath(), input$exportSetting
            )
        },
        message = function(m) {
            shinyjs::html(
                id = "exportSettingStatus", html = m$message, add = TRUE
            )
        })
    })

    # ============================ ANALYSIS FUNCTIONS ==========================

    # * PROFILE CLUSTERING =====================================================
    # ** description for profile clustering function ---------------------------
    observe({
        desc = paste(
            "Cluster genes according to the distance of their phylogenetic
            profiles."
        )

        if (input$tabs == "Profiles clustering") {
            shinyBS::createAlert(
                session, "descClusteringUI", "descClustering",
                content = desc, append = FALSE
            )
        }
    })

    # ** check if genes are added anywhere else to the customized profile ------
    observe({
        if (input$addGeneAgeCustomProfile == TRUE
            | input$addCoreGeneCustomProfile == TRUE
            | input$addGCGenesCustomProfile == TRUE
            | input$addGeneDimRed == TRUE) {
            shinyjs::disable("addClusterCustomProfile")
        } else {
            shinyjs::enable("addClusterCustomProfile")
        }
    })

    output$addClusterCustomProfileCheck.ui <- renderUI({
        if (input$addGeneAgeCustomProfile == TRUE
            | input$addCoreGeneCustomProfile == TRUE
            | input$addGCGenesCustomProfile == TRUE
            | input$addGeneDimRed == TRUE) {
            HTML('<p><em>(Uncheck "Add to Customized profile" check box in
                 <strong>Gene age estimation</strong> or
                 <strong>Core genes finding</strong> or
                 <strong>Group comparison</strong> or
                 <strong>Dimension reduction (Selected genes)</strong>
                 &nbsp;to enable this function)</em></p>')
        }
    })

    # ** List of possible profile types ----------------------------------------
    output$selectProfileType <- renderUI({
        selectedType <- "binary"
        if (!is.null(i_profileType)) selectedType <- i_profileType
        if (input$colorByOrthoID == TRUE) selectedType <- "orthoID"
        if (input$var2ID != "") {
            radioButtons(
                "profileType",
                label = h5("Clustering profiles using"),
                choiceNames = list(
                    "Binary profile", input$var1ID,input$var2ID,"Ortholog IDs"),
                choiceValues = list("binary", "var1", "var2", "orthoID"),
                selected = selectedType,
                inline = FALSE)
        } else {
            if (selectedType == "var2") selectedType <- "binary"
            radioButtons(
                "profileType",
                label = h5("Clustering profiles using"),
                choiceNames =list("binary profile",input$var1ID,"Ortholog IDs"),
                choiceValues = list("binary", "var1", "orthoID"),
                selected = selectedType,
                inline = FALSE)
        }
    })

    # ** List of possible distance methods -------------------------------------
    output$selectDistMethod <- renderUI({
        if (is.null(input$profileType)) profileType <- i_profileType
        else profileType <- input$profileType

        if (profileType == "binary") {
            selectedMethod <- "euclidean"
            if (!is.null(i_distMethod)) selectedMethod <- i_distMethod
            selectInput(
                "distMethod",
                label = h5("Distance measure method:"),
                choices = list(
                    "euclidean" = "euclidean",
                    "maximum" = "maximum",
                    "manhattan" = "manhattan",
                    "canberra" = "canberra",
                    "binary" = "binary",
                    "pearson correlation coefficient" = "pearson",
                    "mutual information" = "mutualInformation",
                    "distance correlation" = "distanceCorrelation"
                ),
                selected = selectedMethod
            )
        } else {
            selectedMethod <- "mutualInformation"
            if (!is.null(i_distMethod)) selectedMethod <- i_distMethod
            selectInput(
                "distMethod",
                label = h5("Distance measure method:"),
                choices = list(
                    "mutual information" = "mutualInformation",
                    "distance correlation" = "distanceCorrelation"
                ),
                selected = selectedMethod
            )
        }
    })

    # ** create profiles for calculating distance matrix -----------------------
    getProfiles <- reactive({
        if (rtCheck) {
            checkpointc101 <- Sys.time()
            print(paste(
                "checkpoint-C101 - before get profiles for clustering", 
                checkpointc101
            ))
        }
        withProgress(message = 'Getting data for cluster...', value = 0.5, {
            req(dataHeat())
            if (is.null(input$profileType)) profileType <- i_profileType
            else profileType <- input$profileType
            if (input$keepOrder == TRUE) {
                fullDt <- getFullData()
                tmpDf <- calcPresSpec(fullDt, getCountTaxa())
                fullDt <- fullDt %>%
                    left_join(tmpDf, by = c("geneID", "supertaxon")) %>%
                    mutate(orthoID = ifelse(
                        presSpec == 0, NA, as.character(orthoID)
                    )) %>% mutate(
                        orthoID = factor(
                            orthoID, levels=unique(as.character(fullDt$orthoID))
                        )
                    )
                if (rtCheck) {
                    checkpointc101a <- Sys.time()
                    print(paste(
                        "checkpoint-C101a - start get profiles for clustering", 
                        checkpointc101a," --- ",checkpointc101a - checkpointc101
                    ))
                }
                profiles <- getDataClustering(
                    fullDt,
                    profileType,
                    input$var1AggregateBy,
                    input$var2AggregateBy
                )
            }
            else {
                if (rtCheck) {
                    checkpointc101a <- Sys.time()
                    print(paste(
                        "checkpoint-C101b - start get profiles for clustering", 
                        checkpointc101a," --- ",checkpointc101a - checkpointc101
                    ))
                }
                profiles <- getDataClustering(
                    dataHeat(),
                    profileType,
                    input$var1AggregateBy,
                    input$var2AggregateBy
                )
            }
            if (rtCheck) {
                checkpointc102 <- Sys.time()
                print(paste(
                    "checkpoint-C102 - get profiles for clustering done", 
                    checkpointc102, " --- ",  checkpointc102 - checkpointc101a
                ))   
            }
            return(profiles)
        })
    })

    # ** calculate distance matrix ---------------------------------------------
    getDistanceMatrixProfiles <- reactive({
        if (rtCheck) {
            checkpointc201 <- Sys.time()
            print(paste(
                "checkpoint-C201 - before calculate dist matrix", checkpointc201
            ))
        }
        withProgress(message = 'Calculating distance matrix...', value = 0.5, {
            req(getProfiles())
            if (is.null(input$distMethod))
                distMethod <- i_distMethod
            else distMethod <- input$distMethod
            distanceMatrix <- getDistanceMatrix(getProfiles(), distMethod)
            if (rtCheck) {
                checkpointc202 <- Sys.time()
                print(paste(
                    "checkpoint-C202 - calculate dist matrix done", 
                    checkpointc202, " --- ",  checkpointc202 - checkpointc201
                ))
            }
            return(distanceMatrix)
        })
    })

    # ** render cluster tree ---------------------------------------------------
    brushedClusterGene <- callModule(
        clusterProfile, "profileClustering",
        distanceMatrix = getDistanceMatrixProfiles,
        clusterMethod = reactive(input$clusterMethod),
        plotWidth = reactive(input$clusterPlot.width),
        plotHeight = reactive(input$clusterPlot.height)
    )

    # * DISTRIBUTION ANALYSIS ==================================================
    # ** description for distribution analysis function ------------------------
    observe({
        desc = paste(
            "Plot the distributions of the values incurred by the integrated
            information layers."
        )

        if (input$tabs == "Distribution analysis") {
            shinyBS::createAlert(
                session, "descDistributionUI", "descDistribution",
                content = desc, append = FALSE
            )
        }
    })

    # ** list of available variables for distribution plot ---------------------
    output$selected.distribution <- renderUI({
        if (nchar(input$var1ID) == 0 & nchar(input$var2ID) == 0) {
            varList <- "% present taxa"
        } else if (nchar(input$var1ID) == 0 & nchar(input$var2ID) > 0) {
            varList <- as.list(c(input$var2ID, "% present taxa"))
        } else if (nchar(input$var1ID) > 0 & nchar(input$var2ID) == 0) {
            varList <- as.list(c(input$var1ID, "% present taxa"))
        } else {
            varList <- as.list(c(input$var1ID, input$var2ID, "% present taxa"))
        }

        selectInput(
            "selectedDist", "Choose variable to plot:", varList, varList[1]
        )
    })

    # ** var1 / var2 distribution data -----------------------------------------
    distributionDf <- reactive({
        req(v$doPlot)
        withProgress(message = 'Getting data for analyzing...', value = 0.5, {
            splitDt <- createVariableDistributionData(
                getMainInput(), input$var1, input$var2
            )
            # filter data base on customized plot (if chosen)
            if (input$dataset.distribution == "Customized data") {
                req(input$inSeq)
                splitDt <- createVariableDistributionDataSubset(
                    getFullData(),
                    splitDt,
                    input$inSeq,
                    input$inTaxa
                )
            }
            # return dt
            return(splitDt)
        })
    })

    # ** render distribution plots ---------------------------------------------
    observe({
        req(v$doPlot)
        req(input$selectedDist)

        if (input$selectedDist == "% present taxa") {
            callModule(
                analyzeDistribution, "distPlot",
                data = reactive(
                    createPercentageDistributionData(
                        getMainInput(), input$rankSelect, getTaxDBpath()
                    )
                ),
                varID = reactive(input$selectedDist),
                varType = reactive("presSpec"),
                percent = reactive(input$percent),
                distTextSize = reactive(input$distTextSize),
                distWidth = reactive(input$distWidth)
            )
        } else {
            if (input$selectedDist == input$var1ID) {
                callModule(
                    analyzeDistribution, "distPlot",
                    data = distributionDf,
                    varID = reactive(input$selectedDist),
                    varType = reactive("var1"),
                    percent = reactive(input$percent),
                    distTextSize = reactive(input$distTextSize),
                    distWidth = reactive(input$distWidth)
                )
            } else if (input$selectedDist == input$var2ID) {
                callModule(
                    analyzeDistribution, "distPlot",
                    data = distributionDf,
                    varID = reactive(input$selectedDist),
                    varType = reactive("var2"),
                    percent = reactive(input$percent),
                    distTextSize = reactive(input$distTextSize),
                    distWidth = reactive(input$distWidth)
                )
            }
        }
    })

    # * GENE AGE ESTIMATION ====================================================
    # ** description for gene age estimation function --------------------------
    observe({
        desc = paste(
            "ESTIMATE THE EVOLUTIONARY AGE OF GENES from the phylogenetic
            profiles using an LCA algorithm. Specifically, the last common
            ancestor of the two most distantly related species displaying
            a given gene serves as the minimal gene age."
        )

        if (input$tabs == "Gene age estimation") {
            shinyBS::createAlert(
                session, "descGeneAgeUI", "descGeneAge",
                title = "", content = desc, append = FALSE
            )
        }
    })

    # ** check if genes are added anywhere else to the customized profile ------
    observe({
        if (input$addClusterCustomProfile == TRUE
            | input$addCoreGeneCustomProfile == TRUE
            | input$addGCGenesCustomProfile == TRUE
            | input$addGeneDimRed == TRUE) {
            shinyjs::disable("addGeneAgeCustomProfile")
        } else {
            shinyjs::enable("addGeneAgeCustomProfile")
        }
    })

    output$addGeneAgeCustomProfileCheck.ui <- renderUI({
        if (input$addClusterCustomProfile == TRUE
            | input$addCoreGeneCustomProfile == TRUE
            | input$addGCGenesCustomProfile == TRUE
            | input$addGeneDimRed == TRUE) {
            HTML('<p><em>(Uncheck "Add to Customized profile" check box in
           <strong>Profile clustering</strong> or
           <strong>Core genes finding</strong> or
           <strong>Group comparison</strong> or
           <strong>Dimension reduction (Selected genes)</strong>
           &nbsp;to enable this function)</em></p>')
        }
    })

    # ** reset geneAgeProtConfig --------------------------------------------
    observeEvent(input$resetGeneAgeProtConfig, {
        shinyjs::reset("geneAgeWidth")
        shinyjs::reset("geneAgeHeight")
        shinyjs::reset("geneAgeText")
    })

    # ** data for gene age estimation ------------------------------------------
    geneAgeDf <- reactive({
        req(v$doPlot)
        withProgress(message = 'Getting data for analyzing...', value = 0.5, {
            geneAgeDf <- estimateGeneAge(
                getFullData(),
                getCountTaxa(),
                toString(input$rankSelect),
                input$inSelect,
                input$var1, input$var2, input$percent, getTaxDBpath()
            )
            return(geneAgeDf)
        })
    })

    # ** render age distribution plot ------------------------------------------
    selectedgeneAge <- callModule(
        plotGeneAge, "geneAge",
        data = geneAgeDf,
        geneAgeWidth = reactive(input$geneAgeWidth),
        geneAgeHeight = reactive(input$geneAgeHeight),
        geneAgeText = reactive(input$geneAgeText),
        font = reactive(input$font)
    )

    # * CORE GENES IDENTIFICATION ==============================================
    # ** description for core gene identification function ---------------------
    observe({
        desc = paste(
            "IDENTIFY GENES THAT ARE SHARED AMONG SELECTED TAXA.",
            "You can set the minimal taxa that should be taken into
            account by using the \"Core taxa coverage\" cutoff.",
            "If you are working with a taxonomy level (e.g. Family)
            that is higher than the one in the input profile (e.g.
            Species), you can also identify a minimal fragtion of species
            that need to have an ortholog in each supertaxon with
            \"% of present taxa\" cutoff. WARNING: You should set the cutoffs
            before selecting taxa of interest!"
        )

        if (input$tabs == "Core gene identification") {
            shinyBS::createAlert(
                session, "descCoreGeneUI", "descCoreGene",
                title = "", content = desc, append = FALSE
            )
        }
    })

    # ** render list of available taxa -----------------------------------------
    output$taxaListCore.ui <- renderUI({
        if (
            input$demoData == "arthropoda" | input$demoData == "ampk-tor" |
            input$demoData == "preCalcDt" | input$mainInputType == "folder"
        ) {
            filein <- 1
        } else filein <- input$mainInput
        if (is.null(filein)) {
            return(selectInput("inTaxa", "Select taxa of interest:", "none"))
        }
        if (v$doPlot == FALSE) {
            return(selectInput("inTaxa", "Select taxa of interest:", "none"))
        } else {
            choice <- inputTaxonName()
            choice$fullName <- as.factor(choice$fullName)

            out <- as.list(levels(choice$fullName))
            out <- append("none", out)

            if (input$applyCoreTaxa == TRUE) {
                out <- coreTaxaName()
                return(selectInput(
                    "taxaCore",
                    "Select taxa of interest:",
                    out,
                    selected = out,
                    multiple = TRUE
                ))
            } else {
                return(selectInput(
                    "taxaCore",
                    "Select taxa of interest:",
                    out,
                    selected = out[1],
                    multiple = TRUE
                ))
            }
        }
    })

    # ** render popup for selecting group of taxa to find core genes -----------
    coreTaxaName <- callModule(
        selectTaxonRank,
        "selectTaxonRankCore",
        rankSelect = reactive(input$rankSelect),
        inputTaxonID = inputTaxonID,
        taxDB = getTaxDBpath
    )

    # ** check if genes are added anywhere else to the customized profile ------
    observe({
        if (input$addClusterCustomProfile == TRUE
            | input$addGeneAgeCustomProfile == TRUE
            | input$addGCGenesCustomProfile == TRUE
            | input$addGeneDimRed == TRUE) {
            shinyjs::disable("addCoreGeneCustomProfile")
        } else {
            shinyjs::enable("addCoreGeneCustomProfile")
        }
    })

    output$addCoreGeneCustomProfileCheck.ui <- renderUI({
        if (input$addClusterCustomProfile == TRUE
            | input$addGeneAgeCustomProfile == TRUE
            | input$addGCGenesCustomProfile == TRUE
            | input$addGeneDimRed == TRUE) {
            HTML('<p><em>(Uncheck "Add to Customized profile" check box in
           <strong>Profiles clustering</strong> or
           <strong>Gene age estimating</strong> or
           <strong>Group Comparioson</strong> or
           <strong>Dimension reduction (Selected genes)</strong>
           &nbsp;to enable this function)</em></p>')
        }
    })

    # ** render table contains list of core genes ------------------------------
    coreGeneDf <- callModule(
        identifyCoreGene,
        "coreGene",
        filteredData = getFullData,
        taxaCount = getCountTaxa,
        rankSelect = reactive(input$rankSelect),
        taxaCore = reactive(input$taxaCore),
        percentCore = reactive(input$percentCore),
        var1Cutoff = reactive(input$var1Core),
        var2Cutoff = reactive(input$var2Core),
        coreCoverage = reactive(input$coreCoverage),
        taxDB = getTaxDBpath
    )

    # ** download gene list from coreGene.table -------------------------------
    output$coreGeneTableDownload <- downloadHandler(
        filename = function() {
            c("coreGeneList.out")
        },
        content = function(file) {
            dataOut <- coreGeneDf()
            write.table(dataOut, file, sep = "\t", row.names = FALSE,
                        quote = FALSE)
        }
    )

    # * GROUP COMPARISON =======================================================
    # ** description for group comparison function -----------------------------
    observe({
        if (is.null(input$var1ID)) return()
        desc = paste("This function is used to COMPARE THE DISTRIBUTIONS of")
        if (input$var1ID == "") {
            desc = paste(desc, "two additional scores")
            shinyjs::disable("plotGC")
        } else if (input$var2ID == "") {
            desc = paste(desc, input$var1ID)
        } else {
            desc = paste(desc, input$var1ID, "and", input$var2ID)
        }
        desc = paste(
            desc,
            "between two taxon groups, an in- and an out-group. You can define
            the in-group below and all taxa not included in this are used as
            the out-group. The value distributions of the variables are then
            compared using statistical tests (Kolmogorov-Smirnov and
            Wilcoxon-Mann-Whitney) using the specified significant level.
            Genes that have a significantly different distribution are
            shown in the candidate gene list below."
        )

        if (input$tabs == "Group comparison") {
            shinyBS::createAlert(
                session, "descGCUI", "descGC", title = "",
                content = desc, append = FALSE
            )
        }
    })

    # ** reset configuration windows for GC plot config ------------------------
    observeEvent(input$resetConfigGC, {
        shinyjs::reset("xSizeGC")
        shinyjs::reset("ySizeGC")
        shinyjs::reset("titleSizeGC")
        shinyjs::reset("legendSizeGC")
        shinyjs::reset("widthVarGC")
        shinyjs::reset("heightVarGC")
        shinyjs::reset("legendGC")
        shinyjs::reset("widthFeatureGC")
        shinyjs::reset("heightFeatureGC")
    })

    observeEvent(input$applyConfigGC, {
        toggleModal(session, "gcPlotConfigBs", toggle = "close")
    })

    # ** check if genes are added anywhere else to the customized profile ------
    observe({
        if (input$addGeneAgeCustomProfile == TRUE
            | input$addCoreGeneCustomProfile == TRUE
            | input$addClusterCustomProfile == TRUE
            | input$addGeneDimRed == TRUE) {
            shinyjs::disable("addGCGenesCustomProfile")
        } else {
            shinyjs::enable("addGCGenesCustomProfile")
        }
    })

    output$addGCCustomProfileCheck <- renderUI({
        if (input$addGeneAgeCustomProfile == TRUE
            | input$addCoreGeneCustomProfile == TRUE
            | input$addClusterCustomProfile == TRUE
            | input$addGeneDimRed == TRUE) {
            HTML(
                '<p><em>(Uncheck "Add to Customized profile" check box in
                 <strong>Gene age estimation</strong> or
                <strong>Profile clustering</strong> or
                <strong>Core genes finding</strong> or
                <strong>Dimension reduction (Selected genes)</strong>
                &nbsp;to enable this function)</em></p>'
            )
        }
    })

    # ** render list of variables ----------------------------------------------
    output$variableGC <- renderUI({
        if (input$var1ID == "") variableList <- list("none" = "none")
        else if (input$var2ID == "")
            variableList <- list("1st Variable" = "var1")
        else
            variableList <- list(
                "1st Variable" = "var1", "2nd Variable" = "var2"
            )
        selectInput(
            inputId = "varNameGC",
            label = "Variable to compare:",
            choices = variableList,
            selected = "var1"
        )
    })

    # ** render list of all sequence IDs (same as customized profile) ----------
    output$listGenesGC <- renderUI({
        fileGC <- input$gcFile
        if (
            input$demoData == "arthropoda" | input$demoData == "ampk-tor" |
            input$demoData == "preCalcDt" | input$mainInputType == "folder"
        ) {
            filein <- 1
        } else filein <- input$mainInput

        if (v$doPlot == FALSE) {
            return(selectInput(
                "selectedGeneGC", "Sequence(s) of interest:", "none"
            ))
        } else {
            data <- as.data.frame(getFullData())
            data$geneID <- as.character(data$geneID)
            data$geneID <- as.factor(data$geneID)
            outAll <- as.list(levels(data$geneID))
            outAll <- append("all", outAll)

            if (is.null(fileGC)) {
                return(selectInput(
                    "selectedGeneGC", "Sequence(s) of interest:",
                    outAll,
                    selected = outAll[1],
                    multiple = TRUE,
                    selectize = FALSE
                ))
            } else {
                listGC <- read.table(file = fileGC$datapath, header = FALSE)
                out <- as.list(levels(listGC$V1))
                return(selectInput(
                    "selectedGeneGC", "Sequence(s) of interest:",
                    out,
                    selected = NULL,
                    multiple = FALSE,
                    selectize = FALSE
                ))
            }
        }
    })

    # ** render popup for selecting rank and return list of belonging taxa -----
    # ** (same as core gene identification)
    gcTaxaName <- callModule(
        selectTaxonRank,
        "selectTaxonRankGC",
        rankSelect = reactive(input$rankSelect),
        inputTaxonID = inputTaxonID,
        taxDB = getTaxDBpath
    )

    # ** check the validity of in-group/out-group taxa input file --------------
    inputTaxonGroupGC <- reactive({
        if (is.null(input$taxonGroupGC)) return()
        taxonGroupGCin <- input$taxonGroupGC
        uploadTaxonGC <- read.table(
            file = taxonGroupGCin$datapath,
            sep = "\t",
            header = FALSE,
            stringsAsFactors = FALSE
        )
        colnames(uploadTaxonGC) <- c("ncbiID", "type")
        return(uploadTaxonGC)
    })

    invalidTaxonGroupGC <- reactive({
        req(inputTaxonGroupGC())
        uploadTaxonGC <- inputTaxonGroupGC()
        # compare with input taxa IDs
        invalidID <- setdiff(uploadTaxonGC$ncbiID, inputTaxonID())
        if (length(invalidID) > 0)
            return(uploadTaxonGC[uploadTaxonGC$ncbiID %in% invalidID,])
    })

    output$checkTaxonGroupGC <- reactive({
        if (is.null(invalidTaxonGroupGC())) return(TRUE)
        else return(FALSE)
    })
    outputOptions(output, "checkTaxonGroupGC", suspendWhenHidden = FALSE)

    output$invalidTaxonGroupGC <- DT::renderDataTable({
        if (is.null(invalidTaxonGroupGC())) return()
        else return(invalidTaxonGroupGC())
    })

    # ** render list of taxa (and default in-group taxa are selected) ----------
    output$taxaListGC <- renderUI({
        if (
            input$demoData == "arthropoda" | input$demoData == "ampk-tor" |
            input$demoData == "preCalcDt" | input$mainInputType == "folder"
        ) {
            filein <- 1
        } else filein <- input$mainInput
        if (is.null(filein)) {
            return(selectInput("selectedInGroupGC", "In-group taxa:", "none"))
        }
        if (v$doPlot == FALSE) {
            return(selectInput("selectedInGroupGC", "In-group taxa:", "none"))
        } else {
            if (is.null(input$taxonGroupGC)) {
                choice <- inputTaxonName()
                choice$fullName <- as.factor(choice$fullName)
                out <- as.list(levels(choice$fullName))

                #' when the taxonomy rank was changed --------------------------
                if (input$applyTaxonGC == TRUE) {
                    out <- gcTaxaName()
                    selectInput(
                        "selectedInGroupGC", "In-group taxa:",
                        out,
                        selected = out,
                        multiple = TRUE,
                        selectize = FALSE
                    )
                }
                #' when the taxonomy is the same as the initially chosen one ---
                else {
                    # all input taxon IDs
                    inputTaxonID <- gsub("ncbi", "", inputTaxonID())
                    # get the next higher rank of the current working rank
                    ranks <- getTaxonomyRanks()
                    pos <- which(ranks == input$rankSelect) # pos in the list
                    higherRank <- ranks[pos + 1] # take the next higher rank
                    higherRankName <- as.character(higherRank[1])
                    # get ID of the selected reference taxon
                    nameList <- getNameList(taxDB = getTaxDBpath())
                    reference <- subset(
                        nameList, nameList$fullName == input$inSelect
                    )
                    # get the corresponding ID in the higher rank for the
                    # selected reference taxon
                    taxMatrix <- getTaxonomyMatrix(
                        getTaxDBpath(), TRUE, inputTaxonID()
                    )
                    higherRankID <- taxMatrix[
                        taxMatrix[, reference$rank] == reference$ncbiID,
                    ][,higherRankName][1]
                    # return selected in-group taxa
                    taxaHigherRank <- getSelectedTaxonNames(
                        inputTaxonID, input$rankSelect,
                        higherRankName, higherRankID, getTaxDBpath()
                    )
                    selectedSupertaxa <- getInputTaxaName(
                        input$rankSelect,
                        paste0("ncbi", taxaHigherRank$ncbiID),
                        getTaxDBpath()
                    )
                    selectInput(
                        "selectedInGroupGC", "Select inGroup taxa:",
                        out,
                        selected = unique(selectedSupertaxa$fullName),
                        multiple = TRUE,
                        selectize = FALSE
                    )
                }
            } else {
                if (is.null(invalidTaxonGroupGC())) {
                    return(
                        selectInput(
                            "selectedInGroupGC", "In-group taxa:", "From file"
                        )
                    )
                } else {
                    return(
                        selectInput(
                            "selectedInGroupGC", "In-group taxa:", "none"
                        )
                    )
                }
            }
        }
    })

    # ** get the ID list of in-group taxa --------------------------------------
    getInGroup <- reactive({
        if (is.null(input$selectedInGroupGC)) return()

        # list of selected in-group taxa names and working rank
        selectedTaxa <- input$selectedInGroupGC
        selectedRank <- input$rankSelect

        if (selectedTaxa[1] == "none") return()
        if (selectedTaxa[1] == "From file") {
            taxonGroupGC <- inputTaxonGroupGC()
            return(taxonGroupGC$ncbiID[taxonGroupGC$type == "in-group"])
        } else {
            # get IDs for selected in-group taxa
            nameList <- getNameList(taxDB = getTaxDBpath())
            selectedTaxaID <- nameList$ncbiID[
                nameList$fullName %in% selectedTaxa
                & nameList$rank == selectedRank]

            # get in-group IDs from raw input (regardless to the working rank)
            taxMatrix <- getTaxonomyMatrix(getTaxDBpath(), TRUE, inputTaxonID())

            inGroup <- as.character(
                taxMatrix$abbrName[
                    taxMatrix[, selectedRank] %in% selectedTaxaID]
            )

            if (length(inGroup) == 0) return()
            else {
                if (input$useCommonAncestor == TRUE) {
                    inGroupTMP <- getCommonAncestor(
                        inputTaxonID(), inGroup, getTaxDBpath()
                    )
                    inGroup <- inGroupTMP[[3]]$abbrName
                }
                return(as.character(inGroup))
            }
        }
    })

    # ** parameters for the plots in Group Comparison --------------------------
    plotParametersGC <- reactive({
        input$updateGC # for trigger changes
        inputData <- list(
            "xSize" = isolate(input$xSizeGC),
            "ySize" = isolate(input$ySizeGC),
            "angle" = isolate(input$angleGC),
            "legendPosition" = isolate(input$legendGC),
            "legendSize" = isolate(input$legendSizeGC),
            "titleSize" = isolate(input$titleSizeGC),
            "flipPlot" = isolate(input$xAxisGC),
            "mValue" = isolate(input$mValueGC),
            "widthVar" = isolate(input$widthVarGC),
            "heightVar" = isolate(input$heightVarGC),
            "widthFeature" = isolate(input$widthFeatureGC),
            "heightFeature" = isolate(input$heightFeatureGC),
            "inGroupName" = isolate(input$inGroupName),
            "outGroupName" = isolate(input$outGroupName)
        )
    })

    # ** data for group comparison ---------------------------------------------
    groupComparisonData <- reactive({
        req(getFullData())
        withProgress(message = 'Getting data for analyzing...', value = 0.5, {
            if (is.null(input$taxonGroupGC)) return(getFullData())
            else {
                taxonGroupGC <- inputTaxonGroupGC()
                dataFiltered <- getFullData()
                return(
                    dataFiltered[dataFiltered$ncbiID %in% taxonGroupGC$ncbiID,]
                )
            }
        })
    })

    # ** render plots for group comparison -------------------------------------
    candidateGenes <- callModule(
        groupComparison, "groupComparison",
        filteredDf = groupComparisonData,
        inGroup = getInGroup,
        variable = reactive(input$varNameGC),
        varName = reactive(c(input$var1ID, input$var2ID)),
        compareType = reactive(input$compareType),
        significanceLevel = reactive(input$significance),
        plotParameters = plotParametersGC,
        domainDf = getDomainInformation,
        doCompare = reactive(input$doCompare),
        doUpdate = reactive(input$updateGC),
        taxDB = getTaxDBpath
    )

    # * WORKING WITH NCBI TAXONOMY DATABASE ====================================
    observe({
        desc = paste(
            "<p><em>PhyloProfile</em> is provided with a set of pre-identified
            taxa (based on the Quest for Ortholog data set). The taxonomy
            information in <em>PhyloProfile</em> is stored in different files
            within the <code>PhyloProfile/PhyloProfile/data</code> folder (to
            check where <em>PhyloProfile</em> package is installed, im R
            Terminal type <code>find.package(\"PhyloProfile\")</code>). Two
            most important files are the <code>preProcessedTaxonomy.txt</code>
            and <code>taxonomyMatrix.txt</code>. The
            <code>preProcessedTaxonomy.txt</code> file stored a pre-processing
            NCBI taxonomy database. While the <code>taxonomyMatrix.txt</code>
            file is used for sorting the input taxa in the profile plot as well
            as dynamically changing the working systematic rank. For more info
            please refer to
            <a href=\"https://github.com/BIONF/PhyloProfile/wiki/PhyloProfile-and-the-NCBI-taxonomy-database\"
            target=\"_blank\">this post</a>.</p>
            <p>If your phylogenetic profiles contains taxa that are not part of
            that set, the new taxa will need to be parsed. Normally, if the new
            taxa can be found in the <code>preProcessedTaxonomy.txt</code>, you
            can easily parse the taxonomy info for those taxa by clicking on
            <strong>Parse taxnomoy info</strong> button in the Shiny App of
            <em>PhyloProfile</em>. In case your pre-processing NCBI taxonomy
            database is out-of-date and some of the new taxa are not in that
            old database, you will see the <strong>Add taxonomy info</strong>
            button instead. You have to either manually add those taxa using
            the <strong>Add taxonomy info</strong> button, or update the
            <code>preProcessedTaxonomy.txt</code> file by using the
            function <strong>\"<span style=\"text-decoration: underline;\">
            Update NCBI taxonomy DB</span>\"</strong>. This
            task will take some minutes depending on your internet
            connection. So, please be patient and wait until the process is
            done!</p>
            <p>Whenever a new taxon is added into the taxonomy data of
            <em>PhyloProfile</em>, the taxonomy files will be
            changed. If you encounter any troubles related to the taxonomy,
            such as error by parsing new taxa, not all input taxa can be found,
            the order of your taxa in the profile plot looks weird, etc., you
            should <strong>\"<span style=\"text-decoration: underline;\">
            Reset the NCBI taxonomy DB</span>\"</strong>.</p>
            <p>You can also <strong>\"<span style=\"text-decoration: underline;\">
            Export the taxonomy database</span>\"</strong> of your current
            working project for backing up.</p>
            <p>Last but not least, if you want to
            <strong>\"<span style=\"text-decoration: underline;\">
            Import any existing taxonomy files</span>\"</strong> from your old
            analysesto <em>PhyloProfile</em> by providing the folder that
            contains these files: <code>newTaxa.txt</code>,
            <code>idList.txt</code>, <code>rankList.txt</code>,
            <code>taxonNamesReduced.txt</code>,
            <code>taxonomyMatrix.txt</code> and
            <code>preCalcTree.nw</code> (optional).</p>"
        )

        if (input$tabs == "NCBI taxonomy data") {
            shinyBS::createAlert(
                session, "descNcbiTaxDbUI", "descNcbiTaxDb", title = "",
                content = desc, append = FALSE
            )
        }
    })

    # ** do update NCBI taxonomy database --------------------------------------
    observeEvent(input$doUpdateNcbi, {
        withCallingHandlers({
            shinyjs::html("updateNCBITaxStatus", "")
            updateNcbiTax()
        },
        message = function(m) {
            shinyjs::html(
                id = "updateNCBITaxStatus", html = m$message, add = TRUE
            )
        })
        shinyBS::updateButton(session, "doUpdateNcbi", disabled = TRUE)
    })

    # ** do reset taxonomy data ------------------------------------------------
    output$taxResetWarning.ui <- renderUI({
        defaultTaxDB <- system.file(
            "PhyloProfile", "data", package = "PhyloProfile", mustWork = TRUE
        )
        em(paste("NOTE: Only reset data in", defaultTaxDB))
    })

    observeEvent(input$doResetTax, {
        withCallingHandlers({
            shinyjs::html("resetTaxonomyDataStatus", "")
            resetTaxData()
        },
        message = function(m) {
            shinyjs::html(
                id = "resetTaxonomyDataStatus", html = m$message, add = TRUE
            )
        })
        shinyBS::updateButton(session, "doResetTax", disabled = TRUE)
    })

    # ** do export taxonomy data -----------------------------------------------
    output$taxExportWarning.ui <- renderUI({
        em(paste("from", getTaxDBpath()))
    })

    getTaxPathOut <- reactive({
        shinyFiles::shinyDirChoose(
            input, "taxDirOut", roots = homePath, session = session
        )
        taxPathOut <- shinyFiles::parseDirPath(homePath, input$taxDirOut)
        return(replaceHomeCharacter(as.character(taxPathOut)))
    })

    output$taxDirOut.ui <- renderUI({
        req(getTaxPathOut())
        if (length(getTaxPathOut()) > 0) {
            em(paste("to", getTaxPathOut()))
        }
    })

    observeEvent(input$doExportTax, {
        req(getTaxPathOut())
        withCallingHandlers({
            shinyjs::html("exportTaxonomyDataStatus", "")
            exportNcbiTax(getTaxPathOut(), getTaxDBpath())
        },
        message = function(m) {
            shinyjs::html(
                id = "exportTaxonomyDataStatus", html = m$message, add = TRUE
            )
        })
        shinyBS::updateButton(session, "doExportTax", disabled = TRUE)
    })

    # ** do import taxonomy data -----------------------------------------------
    output$taxImportWarning.ui <- renderUI({
        msg <- paste(
            "*** WARNING: This function will OVERWRITE the data in",
            getTaxDBpath(),
            "!!!"
        )
        strong(em(msg))
    })

    getInTaxPath <- reactive({
        shinyFiles::shinyDirChoose(
            input, "taxDir", roots = homePath, session = session
        )
        taxPath <- shinyFiles::parseDirPath(homePath, input$taxDir)
        return(replaceHomeCharacter(as.character(taxPath)))
    })

    output$taxDir.ui <- renderUI({
        req(getInTaxPath())
        if (length(getInTaxPath()) > 0) {
            outString <- getInTaxPath()
            em(outString)
        }
    })

    observeEvent(input$doImportTax, {
        req(getInTaxPath())
        withCallingHandlers({
            shinyjs::html("importTaxonomyDataStatus", "")
            importNcbiTax(getInTaxPath(), getTaxDBpath())
        },
        message = function(m) {
            shinyjs::html(
                id = "importTaxonomyDataStatus", html = m$message, add = TRUE
            )
        })
        shinyBS::updateButton(session, "doImportTax", disabled = TRUE)
    })
})
