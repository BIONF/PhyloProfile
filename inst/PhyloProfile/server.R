filterProfileData#' Import function files
sourceFiles = list.files(path = "R",
                          pattern = "*.R$",
                          full.names = TRUE)
lapply(sourceFiles, source, .GlobalEnv)

#' set size limit for input (9999mb)
options(
    shiny.maxRequestSize = 9999 * 1024 ^ 2, # size limit for input 9999mb
    shiny.usecairo = TRUE
)

#' MAIN SERVER ================================================================
shinyServer(function(input, output, session) {
    # Automatically stop a Shiny app when closing the browser tab
    # session$onSessionEnded(stopApp)
    session$allowReconnect(TRUE)

    # =========================== INITIAL CHECKING  ============================

    # * check for internet connection ------------------------------------------
    observe({
        if (hasInternet() == FALSE) {
            toggleState("demoData")
        }
    })

    output$noInternetMsg <- renderUI({
        if (hasInternet() == FALSE) {
            strong(em("Internet connection is required for using demo data!"),
                   style = "color:red")
        } else {
            return()
        }
    })

    # * check for the existence of taxonomy files ------------------------------
    observe({
        if (!file.exists(isolate("data/rankList.txt"))) {
            if (hasInternet() == TRUE) {
                fileUrl <- paste0(
                    "https://raw.githubusercontent.com/BIONF/",
                    "phyloprofile-data/master/rankList.txt"
                )
                ncol <- max(count.fields(fileUrl, sep = "\t"))
                df <- read.table(fileUrl,
                                 sep = "\t",
                                 quote = "",
                                 header = FALSE,
                                 fill = TRUE,
                                 na.strings = c("", "NA"),
                                 col.names = paste0("V", seq_len(ncol)))
                write.table(df, file = "data/rankList.txt",
                            col.names = FALSE,
                            row.names = FALSE,
                            quote = FALSE,
                            sep = "\t") #na = "",
            } else {
                file.create("data/rankList.txt")
            }
        }
    })

    observe({
        if (!file.exists(isolate("data/idList.txt"))) {
            if (hasInternet() == TRUE) {
                fileUrl <- paste0(
                    "https://raw.githubusercontent.com/BIONF/",
                    "phyloprofile-data/master/idList.txt"
                )
                ncol <- max(count.fields(fileUrl, comment.char = "",
                                         sep = "\t"))
                df <- read.table(fileUrl,
                                 sep = "\t",
                                 header = FALSE,
                                 fill = TRUE,
                                 comment.char = "",
                                 na.strings = c("", "NA"),
                                 col.names = paste0("V", seq_len(ncol)))
                write.table(df, file = "data/idList.txt",
                            col.names = FALSE,
                            row.names = FALSE,
                            quote = FALSE,
                            sep = "\t") #na = "",
            } else {
                file.create("data/idList.txt")
            }
        }
    })

    observe({
        if (!file.exists(isolate("data/taxonNamesReduced.txt"))) {
            if (hasInternet() == TRUE) {
                fileUrl <- paste0(
                    "https://raw.githubusercontent.com/BIONF/",
                    "phyloprofile-data/master/taxonNamesReduced.txt"
                )
                ncol <- max(count.fields(fileUrl, sep = "\t"))
                df <- read.table(fileUrl,
                                 sep = "\t",
                                 quote = "",
                                 header = FALSE,
                                 fill = TRUE,
                                 na.strings = c("", "NA"),
                                 col.names = paste0("V", seq_len(ncol)))
                write.table(df, file = "data/taxonNamesReduced.txt",
                            col.names = FALSE,
                            row.names = FALSE,
                            quote = FALSE,
                            sep = "\t")
            } else {
                system("cp data/newTaxa.txt data/taxonNamesReduced.txt")
            }
        }
    })

    observe({
        if (!file.exists(isolate("data/taxonomyMatrix.txt"))) {
            if (hasInternet() == TRUE) {
                fileUrl <- paste0(
                    "https://raw.githubusercontent.com/BIONF/",
                    "phyloprofile-data/master/taxonomyMatrix.txt"
                )
                ncol <- max(count.fields(fileUrl, sep = "\t"))
                df <- read.table(fileUrl,
                                 sep = "\t",
                                 quote = "",
                                 header = FALSE,
                                 fill = TRUE,
                                 na.strings = c("", "NA"),
                                 col.names = paste0("V", seq_len(ncol)))

                write.table(df, file = "data/taxonomyMatrix.txt",
                            col.names = FALSE,
                            row.names = FALSE,
                            quote = FALSE,
                            sep = "\t")
            } else {
                file.create("data/taxonomyMatrix.txt")
            }
        }
    })

    observe({
        if (!file.exists(isolate("data/taxonNamesFull.txt"))) {
            if (hasInternet() == TRUE) {
                fileUrl <- paste0(
                    "https://raw.githubusercontent.com/BIONF/",
                    "phyloprofile-data/master/taxonNamesFull.txt"
                )
                ncol <- max(count.fields(fileUrl, sep = "\t"))
                df <- read.table(fileUrl,
                                 sep = "\t",
                                 quote = "",
                                 header = FALSE,
                                 fill = TRUE,
                                 na.strings = c("", "NA"),
                                 col.names = paste0("V", seq_len(ncol)))
                write.table(df, file = "data/taxonNamesFull.txt",
                            col.names = FALSE,
                            row.names = FALSE,
                            quote = FALSE,
                            sep = "\t")
            } else {
                system("cp data/newTaxa.txt data/taxonNamesFull.txt")
            }
        }
    })

    # ======================== INPUT & SETTINGS TAB ============================
    # * Render input message ---------------------------------------------------
    observe({
        filein <- input$mainInput
        if (is.null(filein) & input$demoData == "none") {
            msg <- paste0(
                "Please <em>upload an input file</em> or
        <em>select a demo data</em><br /> to begin!
        To learn more about the <em>input data</em>, please visit
        <span style=\"text-decoration: underline;\">
        <a href=\"https://github.com/BIONF/PhyloProfile/wiki/Input-Data\">
        <span style=\"color: #ff0000;\">our wiki</span></span></a>."
            )
            createAlert(session, "inputMsgUI", "inputMsg", title = "",
                        content = msg,
                        append = FALSE)
        } else {
            closeAlert(session, "inputMsg")
        }
    })

    # * check the validity of input file and render inputCheck.ui -------------
    output$inputCheck.ui <- renderUI({
        filein <- input$mainInput
        if (is.null(filein)) return()
        inputType <- checkInputValidity(filein$datapath)

        if (inputType[1] == "noGeneID") {
            updateButton(session, "do", disabled = TRUE)
            HTML(
                "<font color=\"red\"><em><strong>ERROR: Unsupported input
                format.<a
                href=\"https://github.com/BIONF/PhyloProfile/wiki/Input-Data\"
                target=\"_blank\">Click here for more
                info</a></em></strong></font>"
            )
        } else if (inputType[1] == "emptyCell") {
            updateButton(session, "do", disabled = TRUE)
            em(strong("ERROR: Rows have unequal length",
                      style = "color:red"))
        }
        else if (inputType[1] == "moreCol") {
            updateButton(session, "do", disabled = TRUE)
            em(strong("ERROR: More columns than column names",
                      style = "color:red"))
        }
        else {
            validType = c("xml", "fasta", "wide", "long", "oma")
            if (!(inputType[1] %in% validType)) {
                updateButton(session, "do", disabled = TRUE)
                invalidOma <- paste(inputType, collapse = "; ")
                msg <- paste0("ERROR: Invalid IDs found! ", invalidOma)
                em(strong(msg,
                          style = "color:red"))
            } else {
                updateButton(session, "do", disabled = FALSE)
                return()
            }
        }
    })

    # * render download link for Demo online files -----------------------------
    output$mainInputFile.ui <- renderUI({
        if (input$demoData == "lca-micros") {
            url <- paste0(
                "https://raw.githubusercontent.com/BIONF/",
                "phyloprofile-data/master/demo/test.main.long"
            )
            strong(a("Download demo input file",
                     href = url,
                     target = "_blank"))
        } else if (input$demoData == "ampk-tor") {
            url <- paste0(
                "https://raw.githubusercontent.com/BIONF/phyloprofile-data/",
                "master/expTestData/ampk-tor/ampk-tor.phyloprofile"
            )
            strong(a("Download demo input file",
                     href = url,
                     target = "_blank"))
        } else {
            fileInput("mainInput", h5("Upload input file:"))
        }
    })

    output$domainInputFile.ui <- renderUI({
        if (input$demoData == "lca-micros") {
            url <- paste0(
                "https://raw.githubusercontent.com/BIONF/phyloprofile-data/",
                "master/demo/domainFiles"
            )
            strong(a("Download demo domain files",
                     href = url,
                     target = "_blank"))
        } else if (input$demoData == "ampk-tor") {
            url <- paste0(
                "https://raw.githubusercontent.com/BIONF/phyloprofile-data/",
                "master/expTestData/ampk-tor/ampk-tor.domains_F"
            )
            strong(a("Download demo domain file",
                     href = url,
                     target = "_blank"))
        } else {
            if (input$annoLocation == "from file") {
                fileInput("fileDomainInput", "")
            } else {
                textInput("domainPath", "", "")
            }
        }
    })

    output$downloadFastaDemo.ui <- renderUI({
        if (input$demoData == "lca-micros") {
            url <- paste0(
                "https://raw.githubusercontent.com/BIONF/phyloprofile-data/",
                "master/demo/fastaFile/concatenatedSeq.fa"
            )
            strong(a("Download demo fasta file",
                     href = url,
                     target = "_blank"))
        } else if (input$demoData == "ampk-tor") {
            url <- paste0(
                "https://raw.githubusercontent.com/BIONF/phyloprofile-data/",
                "master/expTestData/ampk-tor/ampk-tor.extended.fa"
            )
            strong(a("Download demo fasta file",
                     href = url,
                     target = "_blank"))
        }
    })

    # * render description for Demo data ---------------------------------------
    output$demoDataDescribe <- renderUI({
        if (input$demoData == "none") {
            return()
        } else if (input$demoData == "ampk-tor") {
            url <- paste0(
                "https://github.com/BIONF/phyloprofile-data/blob/master/",
                "expTestData/ampk-tor/README.md"
            )
            em(a("Data description",
                 href = url,
                 target = "_blank"))
        } else {
            url <- paste0(
                "https://github.com/BIONF/phyloprofile-data/blob/master/",
                "demo/README.md"
            )
            em(a("Data description",
                 href = url,
                 target = "_blank"))
        }
    })

    # * check OMA input --------------------------------------------------------
    output$checkOmaInput <- reactive({
        filein <- input$mainInput
        if (is.null(filein)) return()
        inputType <- checkInputValidity(filein$datapath)
        inputType == "oma"
    })
    outputOptions(output, "checkOmaInput", suspendWhenHidden = FALSE)

    # * download OMA data after parsing ----------------------------------------
    output$downloadFilesOma <- downloadHandler(
        filenname <- function() {
            "omaDataToPhyloprofileInput.zip"
        },
        content <- function(file) {
            write.table(getMainInput(), "phyloprofile.txt",
                        sep = "\t",
                        row.names = FALSE,
                        col.names = TRUE,
                        quote = FALSE)

            write.table(getAllFastaOma(finalOmaDf()), "fasta.txt",
                        sep = "\t",
                        row.names = FALSE,
                        col.names = FALSE,
                        quote = FALSE)

            write.table(getDomainInformation(), "domain.txt",
                        sep = "\t",
                        row.names = FALSE,
                        col.names = FALSE,
                        quote = FALSE)

            zip(zipfile = file,
                files = c("phyloprofile.txt", "domain.txt", "fasta.txt"))
        },
        contentType = "application/zip"
    )

    # * close OMA parsing popup windows -------------------------------------
    observeEvent(input$getDataOma, {
        toggleModal(session, "getOmaDataWindows", toggle = "close")
        updateButton(session, "getDataOma", disabled = TRUE)
        toggleState("mainInput")
        toggleState("fileDomainInput")
        toggleState("fastaUpload")
    })

    # * render textinput for Variable 1 & 2 ------------------------------------
    output$var1ID.ui <- renderUI({
        longDataframe <- getMainInput()
        if (is.null(longDataframe)) {
            textInput("var1ID",
                      h5("1st variable:"),
                      value = "Variable 1",
                      width = "100%",
                      placeholder = "Name of first variable")
        } else {
            textInput("var1ID", h5("1st variable:"),
                      value = colnames(longDataframe)[4],
                      width = "100%",
                      placeholder = "Name of first variable")
        }
    })

    output$var2ID.ui <- renderUI({
        longDataframe <- getMainInput()
        if (is.null(longDataframe)) {
            textInput("var2ID",
                      h5("2st variable:"),
                      value = "Variable 2",
                      width = "100%",
                      placeholder = "Name of second variable")
        } else {
            textInput("var2ID", h5("2st variable:"),
                      value = colnames(longDataframe)[5],
                      width = "100%",
                      placeholder = "Name of second variable")
        }
    })

    # * render 2. variable relationship according to demo data -----------------
    output$var2Relation.ui <- renderUI({
        if (input$demoData == "ampk-tor") {
            selectInput("var2Relation", label = h5("Relationship:"),
                        choices = list("Prot-Prot" = "protein",
                                       "Prot-Spec" = "species"),
                        selected = "protein",
                        width = 130)
        } else {
            selectInput("var2Relation", label = h5("Relationship:"),
                        choices = list("Prot-Prot" = "protein",
                                       "Prot-Spec" = "species"),
                        selected = "species",
                        width = 130)
        }
    })

    # * check the existance of the input concatenate fasta file ----------------
    output$concatFasta.existCheck <- renderUI({
        if (is.null(input$concatFasta)) return()
        else {
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
                        HTML('<p><span style="color: #0000ff;">
                 <strong>Please click CLOSE to comfirm!</strong></span></p>')
                    } else {
                        helpText("is not a fasta file!!")
                    }
                }
            }
        }
    })

    # * check the validity of input tree file and render checkNewick.ui --------
    checkNewickID <- reactive({
        req(input$inputTree)
        req(input$mainInput)

        filein <- input$inputTree
        tree <- read.table(
            file = filein$datapath,
            header = FALSE,
            check.names = FALSE,
            comment.char = "",
            fill = FALSE
        )

        checkNewick <- checkNewick(tree, inputTaxonID())
        if (checkNewick == 0) {
            updateButton(session, "do", disabled = FALSE)
        }
        return(checkNewick)
    })

    output$checkNewick.ui <- renderUI({
        checkNewick <- checkNewickID()
        if (checkNewick == 1) {
            updateButton(session, "do", disabled = TRUE)
            HTML("<p><em><span style=\"color: #ff0000;\"><strong>
           ERROR: Parenthesis(-es) missing!</strong></span></em></p>")
        } else if (checkNewick == 2) {
            updateButton(session, "do", disabled = TRUE)
            HTML("<p><em><span style=\"color: #ff0000;\"><strong>
           ERROR: Comma(s) missing!</strong></span></em></p>")
        } else if (checkNewick == 3) {
            updateButton(session, "do", disabled = TRUE)
            HTML("<p><em><span style=\"color: #ff0000;\"><strong>
           ERROR: Tree contains singleton!</strong></span></em></p>")
        } else if (checkNewick == 0) {
            return()
        } else {
            updateButton(session, "do", disabled = TRUE)
            strong(em(paste0(checkNewick, " not exist in main input file!")),
                   style = "color:red")
        }
    })

    # * reset profile plot colors ----------------------------------------------
    observeEvent(input$defaultColorVar2, {
        shinyjs::reset("lowColorVar2")
        shinyjs::reset("highColorVar2")
    })

    observeEvent(input$defaultColorVar1, {
        shinyjs::reset("lowColorVar1")
        shinyjs::reset("highColorVar1")
    })

    observeEvent(input$defaultColorPara, {
        shinyjs::reset("paraColor")
    })

    # * render list of taxonomy ranks ------------------------------------------
    output$rankSelect <- renderUI({
        if (input$demoData == "lca-micros") {
            selectInput("rankSelect", label = "",
                        choices = getTaxonomyRanks(),
                        selected = "phylum")
        } else if (input$demoData == "ampk-tor") {
            selectInput("rankSelect", label = "",
                        choices = getTaxonomyRanks(),
                        selected = "species")
        } else {
            selectInput("rankSelect", label = "",
                        choices = getTaxonomyRanks(),
                        selected = "species")
        }
    })

    # * render list of (super)taxa ---------------------------------------------
    output$select <- renderUI({
        choice <- inputTaxonName()
        choice$fullName <- as.factor(choice$fullName)

        if (input$demoData == "lca-micros") {
            hellemDf <- data.frame("name" = c("Encephalitozoon hellem",
                                              "Encephalitozoon hellem",
                                              "Encephalitozoon",
                                              "Unikaryonidae",
                                              "Apansporoblastina",
                                              "Apansporoblastina",
                                              "Microsporidia",
                                              "Fungi",
                                              "Eukaryota"),
                                   "rank" = c("strain",
                                              "species",
                                              "genus",
                                              "family",
                                              "order",
                                              "class",
                                              "phylum",
                                              "kingdom",
                                              "superkingdom"))
            # rankSelect <- input$rankSelect
            # rankName <- substr(rankSelect,
            #                    4,
            #                    nchar(rankSelect))
            rankName <- input$rankSelect

            selectInput("inSelect", "",
                        as.list(levels(choice$fullName)),
                        hellemDf$name[hellemDf$rank == rankName])
        } else if (input$demoData == "ampk-tor") {
            humanDf <- data.frame("name" = c("Homo sapiens",
                                             "Homo sapiens",
                                             "Homo",
                                             "Hominidae",
                                             "Primates",
                                             "Mammalia",
                                             "Chordata",
                                             "Metazoa",
                                             "Eukaryota"),
                                  "rank" = c("strain",
                                             "species",
                                             "genus",
                                             "family",
                                             "order",
                                             "class",
                                             "phylum",
                                             "kingdom",
                                             "superkingdom"))
            # rankSelect <- input$rankSelect
            # rankName <- substr(rankSelect, 4, nchar(rankSelect))
            rankName <- input$rankSelect

            selectInput("inSelect", "",
                        as.list(levels(choice$fullName)),
                        humanDf$name[humanDf$rank == rankName])
        } else {
            selectInput("inSelect", "",
                        as.list(levels(choice$fullName)),
                        levels(choice$fullName)[1])
        }
    })

    # * enable "PLOT" button ---------------------------------------------------
    observeEvent(input$rankSelect,  ({
        if (input$rankSelect == "") {
            updateButton(session, "do", disabled = TRUE)
        } else {
            unkTaxa <- unkTaxa()
            if (length(unkTaxa) == 0) {
                updateButton(session, "do", disabled = FALSE)
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
            toggleState("mainInput")
            toggleState("geneListSelected")
            toggleState("demoData")
        }
    })

    # * update var2AggregateBy to mean if using demo lca-micros data ---------
    observe({
        if (input$demoData == "lca-micros") {
            ### update var2AggregateBy to mean
            updateSelectInput(session, "var2AggregateBy",
                              choices = list("Max" = "max",
                                             "Min" = "min",
                                             "Mean" = "mean",
                                             "Median" = "median"),
                              selected = "mean")
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
        createSliderCutoff(
            "var1cus",
            paste(input$var1ID, "cutoff:"),
            input$var1[1], input$var1[2], input$var1ID
        )
    })

    output$var2Filter.ui <- renderUI({
        createSliderCutoff(
            "var2cus",
            paste(input$var2ID, "cutoff:"),
            input$var2[1], input$var2[2], input$var2ID
        )
    })

    output$percentFilter.ui <- renderUI({
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
        longDataframe <- getMainInput()
        req(longDataframe)

        if (is.null(longDataframe)) {
            inputTaxa <- c("NA")
        } else {
            inputTaxa <- levels(longDataframe$ncbiID)
        }

        if (inputTaxa[1] == "NA") {
            return()
        } else {
            inputTaxa <- unlist(strsplit(inputTaxa, split = "\t"))
            if (inputTaxa[1] == "geneID") {
                # remove "geneID" element from vector inputTaxa
                inputTaxa <- inputTaxa[-1]
            }

            if (!file.exists(isolate("data/rankList.txt"))) {
                return(inputTaxa)
            } else {
                info <- file.info("data/rankList.txt")
                if (info$size == 0) {
                    return(inputTaxa)
                } else {
                    rankListFile <- paste0(getwd(), "/data/rankList.txt")
                    allTaxa <- as.factor(
                        unlist(
                            data.table::fread(
                                file = rankListFile, select = 1
                            )
                        )
                    )

                    # list of unknown taxa
                    unkTaxa <- inputTaxa[!(inputTaxa %in% allTaxa)]
                    if (identical(unkTaxa, character(0))) return()

                    # get non-ncbi taxa
                    unkTaxa <- data.frame("TaxonID" = unkTaxa)
                    unkTaxa$id <- substring(unkTaxa$TaxonID, 5)
                    unkTaxa$Source <- "ncbi"

                    nameFullFile <- paste0(getwd(), "/data/taxonNamesFull.txt")
                    ncbiTaxa <- as.factor(
                        unlist(
                            data.table::fread(
                                file = nameFullFile, select = 1
                            )
                        )
                    )

                    ncbiID <- levels(ncbiTaxa)
                    maxNCBI <- max(sort(as.numeric(ncbiID[ncbiID != "ncbiID"])))

                    if (nrow(unkTaxa[!(unkTaxa$id %in% ncbiTaxa),]) > 0) {
                        unkTaxa <- unkTaxa[!(unkTaxa$id %in% ncbiTaxa),]$id
                        unkTaxa[unkTaxa$id %in% unkTaxa,]$Source <- "unknown"
                        if (any(unkTaxa < maxNCBI)) {
                            unkTaxa[unkTaxa$id %in% unkTaxa &
                                    unkTaxa$id < maxNCBI,]$Source <- "invalid"
                        }
                    }

                    newTaxaFile <- paste0(getwd(), "/data/newTaxa.txt")
                    newTaxa <- as.factor(
                        unlist(
                            data.table::fread(
                                file = newTaxaFile, select = 1
                            )
                        )
                    )

                    if (nrow(unkTaxa[unkTaxa$id %in% newTaxa,]) > 0) {
                        unkTaxa[unkTaxa$id %in% newTaxa,]$Source <- "new"
                    }

                    # check for invalid newly generated IDs in newTaxa.txt file
                    if (length(newTaxa) > 1) {
                        newTaxaList <- levels(newTaxa)
                        newTaxaList <- as.integer(
                            newTaxaList[newTaxaList != "ncbiID"]
                        )

                        if (min(newTaxaList) < maxNCBI) {
                            invalidList <- as.data.frame(
                                newTaxaList[newTaxaList < maxNCBI]
                            )
                            colnames(invalidList) <- c("id")
                            invalidList$Source <- "newTaxa.txt"
                            invalidList$TaxonID <- "invalid"
                            unkTaxa <- rbind(invalidList, unkTaxa)
                        }
                    }

                    # return list of unkTaxa
                    return(unkTaxa)
                }
            }
        }
        # return input taxa
        return(inputTaxa)
    })

    # * check the status of unkTaxa --------------------------------------------
    output$unkTaxaStatus <- reactive({
        unkTaxa <- unkTaxa()
        if (length(unkTaxa) > 0) {
            if ("invalid" %in% unkTaxa$TaxonID) return("invalid")
            if ("unknown" %in% unkTaxa$Source) return("unknown")
            else return("ncbi")
        } else {
            return(0)
        }
    })
    outputOptions(output, "unkTaxaStatus", suspendWhenHidden = FALSE)

    # * render list of unkTaxa -------------------------------------------------
    output$unkTaxaFull <-
        renderDataTable(options = list(searching = FALSE, pageLength = 10),{
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
            write.table(dataOut, file,
                        sep = "\t",
                        row.names = FALSE,
                        quote = FALSE)
        }
    )

    # * update the form for adding new taxa ------------------------------------
    newTaxa <- reactiveValues()
    newTaxa$Df <- data.frame("ncbiID" = numeric(),
                             "fullName" = character(),
                             "rank" = character(),
                             "parentID" = numeric(),
                             stringsAsFactors = FALSE)
    newIndex <- reactiveValues()
    newIndex$value <- 1

    observeEvent(input$newAdd, {
        newTaxa$Df[newIndex$value, ] <- c(input$newID,
                                          input$newName,
                                          input$newRank,
                                          input$newParent)
        newIndex$value <- newIndex$value + 1
        updateTextInput(session, "newID",
                        value = as.numeric(input$newID) + 1)
        updateTextInput(session, "newName", value = "")
        updateTextInput(session, "newRank", value = "norank")
        updateTextInput(session, "newParent", value = "")
        shinyjs::enable("newDone")
    })

    # * get info for new taxa from uploaded file -------------------------------
    newTaxaFromFile <- reactive({
        filein <- input$newTaxaFile
        if (is.null(filein)) return()
        tmpDf <- read.table(
            file = filein$datapath,
            sep = "\t",
            header = TRUE,
            check.names = FALSE,
            comment.char = ""
        )
        if (ncol(tmpDf) != 4) {
            createAlert(session, "wrongNewTaxa", "wrongNewTaxaMsg",
                        content = "Wrong format. Please check your file!",
                        append = FALSE)
            shinyjs::disable("newDone")
            return()
        } else {
            createAlert(session, "wrongNewTaxa", "wrongNewTaxaMsg",
                        content = "Click Finish adding to continue!",
                        append = FALSE)
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
        write.table(newTaxa$Df, "data/newTaxa.txt",
                    sep = "\t",
                    eol = "\n",
                    row.names = FALSE,
                    quote = FALSE)
    })

    # * check if data is loaded and "parse" button is clicked and confirmed ----
    v1 <- reactiveValues(parse = FALSE)
    observeEvent(input$butParse, {
        toggleModal(session, "parseConfirm", toggle = "close")
        v1$parse <- input$butParse
        updateButton(session, "butParse", disabled = TRUE)
        toggleState("newTaxaAsk")
        toggleState("mainInput")
    })

    # * create rankList.txt, idList.txt, ---------------------------------------
    invalidID <- reactive({
        invalidID <- data.frame("id" = as.character(),
                                "type" = as.character(),
                                stringsAsFactors = FALSE)

        filein <- input$mainInput
        if (is.null(filein)) return()
        inputType <- checkInputValidity(filein$datapath)

        if (inputType == "xml" |
            inputType == "long" |
            inputType == "wide" |
            inputType == "fasta" |
            inputType == "oma") {
            inputDf <- read.table(
                file = filein$datapath,
                sep = "\t",
                header = TRUE,
                check.names = FALSE,
                comment.char = ""
            )

            if (v1$parse == FALSE) return()
            else {
                # get list of taxa need to be parsed (taxa mising taxonomy info)
                if (v1$parse == TRUE) {
                    unkTaxaDf <- unkTaxa()
                    unkTaxa <- as.character(substring(unkTaxaDf$TaxonID, 5))
                    titleline <- c("geneID", unkTaxa)
                }
                # invalidIDtmp <- list()

                ncbiTaxonInfo <- fread("data/taxonNamesFull.txt")
                newTaxaFromFile <- fread("data/newTaxa.txt",
                                         colClasses = c("ncbiID" = "character"))

                ## join all ncbi taxa and new taxa together
                allTaxonInfo <- rbind(newTaxaFromFile, ncbiTaxonInfo)

                ## check missing ids
                if (any(!(unkTaxa %in% allTaxonInfo$ncbiID))) {
                    invalidMissing <-
                        unkTaxa[!(unkTaxa %in% allTaxonInfo$ncbiID)]
                    invalidIDTmp <- data.frame(
                        "id" = invalidMissing,
                        "type" = rep("missing", length(invalidMissing))
                    )
                    invalidID <- rbind(invalidID, invalidIDTmp)
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
                    invalidID <- rbind(invalidID, invalidIDTmp)

                    newTaxaFromFile <- newTaxaFromFile[
                        !(newTaxaFromFile$ncbiID %in% ncbiTaxonInfo$ncbiID),
                    ]
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
                    invalidID <- rbind(invalidID, invalidIDTmp)
                }

                if (nrow(invalidID) > 0) {
                    return(invalidID)
                }

                ## parse taxonomy info
                withProgress(
                    message = "Parsing new taxa...",
                    value = 0, {
                        taxonomyInfo <- getIDsRank(
                            titleline[2:length(titleline)], allTaxonInfo
                        )
                        rankList <- as.data.frame(taxonomyInfo[2])
                        idList <- as.data.frame(taxonomyInfo[1])
                        reducedInfoList <- as.data.frame(taxonomyInfo[3])
                    }
                )

                withProgress(message = "Generating taxonomy file...",
                             value = 0, {
                    # open existing files
                    # (idList, rankList and taxonNamesReduced.txt)
                    ncol <- max(count.fields("data/rankList.txt", sep = "\t"))
                    oldIDList <- read.table(
                        "data/idList.txt",
                        sep = "\t",
                        header = FALSE,
                        check.names = FALSE,
                        comment.char = "",
                        fill = TRUE,
                        stringsAsFactors = TRUE,
                        na.strings = c("", "NA"),
                        col.names = paste0("X", seq_len(ncol))
                    )

                    oldRankList <- read.table(
                        "data/rankList.txt",
                        sep = "\t",
                        header = FALSE,
                        check.names = FALSE,
                        comment.char = "",
                        fill = TRUE,
                        stringsAsFactors = TRUE,
                        na.strings = c("", "NA"),
                        col.names = paste0("X", seq_len(ncol))
                    )

                    oldNameList <- read.table(
                        "data/taxonNamesReduced.txt",
                        sep = "\t",
                        header = TRUE,
                        check.names = FALSE,
                        comment.char = "",
                        fill = TRUE,
                        stringsAsFactors = TRUE
                    )

                    # and append new info into those files
                    newIDList <- rbind.fill(oldIDList, idList)
                    newRankList <- rbind.fill(oldRankList, rankList)
                    newNameList <- rbind.fill(oldNameList, reducedInfoList)

                    # write output files
                    # (idList, rankList and taxonNamesReduced)
                    write.table(newIDList[!duplicated(newIDList), ],
                                file  = "data/idList.txt",
                                col.names = FALSE,
                                row.names = FALSE,
                                quote = FALSE,
                                sep = "\t")
                    write.table(newRankList[!duplicated(newRankList), ],
                                file = "data/rankList.txt",
                                col.names = FALSE,
                                row.names = FALSE,
                                quote = FALSE,
                                sep = "\t")
                    write.table(newNameList[!duplicated(newNameList), ],
                                file = "data/taxonNamesReduced.txt",
                                col.names = TRUE,
                                row.names = FALSE,
                                quote = FALSE,
                                sep = "\t")

                    # create taxonomy matrix (taxonomyMatrix.txt)
                    taxMatrix <- taxonomyTableCreator(
                        "data/idList.txt", "data/rankList.txt")
                    write.table(taxMatrix,
                                file = "data/taxonomyMatrix.txt",
                                sep = "\t",
                                eol = "\n",
                                row.names = FALSE,
                                quote = FALSE)
                })
            }
        }
        return()
    })

    # * output invalid NCBI ID -------------------------------------------------
    output$invalidID.output <- renderTable({
        if (is.null(invalidID())) return()
        else {
            outDf <- invalidID()
            colnames(outDf) <- c("Invalid ID(s)", "Type")
            return(outDf)
        }
    })

    # * download list of invalidID ---------------------------------------------
    output$invalidID.download <- downloadHandler(
        filename = function() {
            c("invalidIDs.txt")
        },
        content = function(file) {
            dataOut <- invalidID()
            colnames(dataOut) <- c("Invalid ID(s)", "Type")
            write.table(dataOut, file,
                        sep = "\t",
                        row.names = FALSE,
                        quote = FALSE)
        }
    )

    # * render final msg after taxon parsing -----------------------------------
    output$endParsingMsg <- renderUI({
        if (is.null(invalidID())) {
            strong(h4("PLEASE RELOAD THIS TOOL WHEN FINISHED!!!"),
                   style = "color:red")
        } else {
            HTML('<p><strong><span style="color: #e12525;"> SOME INVALID TAXON
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
           <p>For IDs with type of <em><span style="color: #0000ff;">"missing"
           </span></em>, please check the validity of them&nbsp;in
           <a href="https://www.ncbi.nlm.nih.gov/taxonomy" target="_blank"
           rel="noopener"> NCBI taxonomy database</a>!</p>')
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
        if (is.null(filein) & input$demoData == "none") {
            v$doPlot <- FALSE
            updateButton(session, "do", disabled = TRUE)
        }
    })

    # * check if "no ordering gene IDs" has been checked -----------------------
    output$applyClusterCheck.ui <- renderUI({
        if (input$ordering == FALSE) {
            HTML('<p><em>(Check "Ordering sequence IDs" check box in
           <strong>Input & settings tab</strong>&nbsp;to enable this function)
           </em></p>')
        }
    })

    # * to enable clustering ---------------------------------------------------
    observe({
        if (input$ordering == FALSE) {
            shinyjs::disable("applyCluster")
        } else {
            shinyjs::enable("applyCluster")
        }
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
        return(data.frame(rbindlist(tmp, use.names = TRUE)))
    }

    finalOmaDf <- reactive({
        filein <- input$mainInput
        if (is.null(filein)) return()
        inputType <- checkInputValidity(filein$datapath)

        if (inputType == "oma") {
            if (input$getDataOma[1] == 0) return()
            omaIDs <- read.table(
                file = filein$datapath,
                sep = "\t",
                header = FALSE,
                check.names = FALSE,
                comment.char = ""
            )
            omaIDs[,1] <- as.character(omaIDs[,1])
            finalOmaDf <- getOmaBrowser(omaIDs[,1],
                                            input$selectedOmaType)
            return(finalOmaDf)
        } else {
            return()
        }
    })

    # * convert main input file in any format into long format dataframe -------
    getMainInput <- reactive({
        if (input$demoData == "lca-micros") {
            longDataframe <- createLongMatrix("lca-micros")
        } else if (input$demoData == "ampk-tor") {
            longDataframe <- createLongMatrix("ampk-tor")
        } else {
            filein <- input$mainInput
            if (is.null(filein)) return()
            inputType <- checkInputValidity(filein$datapath)
            if (inputType == "oma") {
                if (input$getDataOma[1] == 0) return()
                longDataframe <- createProfileFromOma(finalOmaDf())
                longDataframe <- as.data.frame(unclass(longDataframe))
            } else {
                longDataframe <- createLongMatrix(filein$datapath)
            }
        }
        return(longDataframe)
    })

    # * parse domain info into data frame --------------------------------------
    getDomainInformation <- reactive({
        if (input$demoData == "none") {
            filein <- input$mainInput
            inputType <- checkInputValidity(filein$datapath)
        } else {
            inputType <- "demo"
        }

        if (inputType == "oma") {
            domainDf <- getAllDomainsOma(finalOmaDf())
        } else {
            message("Getting the domains...")
            mainInput <- getMainInput()

            if (inputType == "demo") {
                domainDf <- parseDomainInput(
                    unlist(mainInput$geneID),
                    input$demoData,
                    "demo"
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

    # * get ID list of input taxa from main input ------------------------------
    inputTaxonID <- reactive({
        if (input$demoData == "lca-micros" |
            input$demoData == "ampk-tor" |
            length(unkTaxa()) == 0) {
            longDataframe <- getMainInput()
            inputTaxa <- getInputTaxaID(longDataframe)
        } else {
            return()
        }
    })

    # * get NAME list of all (super)taxa ---------------------------------------
    inputTaxonName <- reactive({
        filein <- input$mainInput
        if (is.null(filein) & input$demoData == "none") return()

        if (length(unkTaxa()) > 0) return()
        # rankSelect <- input$rankSelect
        # if (rankSelect == "") return()

        # get rank name from rankSelect
        # rankName <- substr(rankSelect, 4, nchar(rankSelect))
        rankName <- input$rankSelect
        if (rankName == "") return()
        inputTaxaName <- getInputTaxaName(rankName, inputTaxonID())

        return(inputTaxaName)
    })

    # * sort taxonomy data of input taxa ---------------------------------------
    sortedtaxaList <- reactive({
        if (v$doPlot == FALSE) return()

        # get selected rank
        # rankSelect <- input$rankSelect
        # rankName <- substr(rankSelect, 4, nchar(rankSelect))
        rankName <- input$rankSelect

        # get input taxonomy tree
        inputTaxaTree <- NULL
        treeIn <- input$inputTree
        if (!is.null(treeIn)) inputTaxaTree <- read.tree(file = treeIn$datapath)

        # sort taxonomy matrix based on selected refTaxon
        sortedOut <- sortInputTaxa(
            taxonIDs = inputTaxonID(),
            taxonNames = inputTaxonName(),
            rankName = rankName,
            refTaxon = input$inSelect,
            taxaTree = inputTaxaTree
        )
        # return
        return(sortedOut)
    })

    # * get subset data (default: first 30 genes) for plotting -----------------
    preData <- reactive({
        # isolate start and end gene index
        input$updateBtn

        if (input$autoUpdate == TRUE) {
            startIndex <- input$stIndex
            endIndex <- input$endIndex
        } else {
            startIndex <- isolate(input$stIndex)
            endIndex <- isolate(input$endIndex)
        }

        # get list of gene of interest (from a separated file)
        listGene <- list()
        if (is.na(endIndex)) endIndex <- 30

        # get list of genes based on selected gene index
        if (input$geneListSelected == "from file") {
            listIn <- input$list
            if (!is.null(listIn)) {
                list <- read.table(file = listIn$datapath, header = FALSE)
                listGeneOri <- list$V1
                if (startIndex <= length(listGeneOri)) {
                    listGene <-
                        listGeneOri[listGeneOri[startIndex:endIndex]]
                } else {
                    listGene <- listGeneOri
                }
            }
        }

        longDataframe <- getMainInput()
        if (is.null(longDataframe)) return()

        if (is.null(longDataframe)) {
            data <- data.frame("geneID" = character(),
                               "ncbiID" = character(),
                               "orthoID" = character(),
                               "var1" = character(),
                               "var2" = character(),
                               stringsAsFactors = FALSE)
        } else {
            longDataframe <- unsortID(longDataframe, input$ordering)

            if (length(listGene) >= 1) {
                data <- longDataframe[longDataframe$geneID %in% listGene, ]
            } else {
                subsetID <-
                    levels(longDataframe$geneID)[startIndex:endIndex]
                data <- longDataframe[longDataframe$geneID %in% subsetID, ]
            }

            if (ncol(data) < 5) {
                for (i in seq_len(5 - ncol(data))) {
                    data[paste0("newVar", i)] <- 1
                }
            }

            colnames(data) <- c("geneID", "ncbiID", "orthoID", "var1", "var2")
        }
        # return preData
        return(data)
    })

    # * creating main dataframe for subset taxa (in species/strain level) ------
    # * get (super)taxa names (1)
    # * calculate percentage of presence (2),
    # * max/min/mean/median VAR1 (3) and VAR2 (4)
    getDataFiltered <- reactive({
        if (is.null(preData())) return()

        fullMdData <- parseInfoProfile(
            inputDf = preData(),
            sortedInputTaxa = sortedtaxaList(),
            var1AggregateBy = input$var1AggregateBy,
            var2AggregateBy = input$var2AggregateBy
        )

        return(fullMdData)
    })

    # * reduce data from lowest level to supertaxon (e.g. phylum) --------------
    # * This data set contain only supertaxa
    # * and their value (%present, mVar1 & mVar2) for each gene
    dataSupertaxa <- reactive({
        fullMdData <- getDataFiltered()
        superDfExt <- reduceProfile(fullMdData)
        return(superDfExt)
    })

    # * heatmap data input -----------------------------------------------------
    dataHeat <- reactive({
        {
            input$plotCustom
            input$updateBtn
        }

        # get all cutoffs
        if (input$autoUpdate == TRUE) {
            percentCutoff <- input$percent
            coorthologCutoffMax <- input$coortholog
            var1Cutoff <- input$var1
            var2Cutoff <- input$var2
        } else {
            percentCutoff <- isolate(input$percent)
            coorthologCutoffMax <- isolate(input$coortholog)
            var1Cutoff <- isolate(input$var1)
            var2Cutoff <- isolate(input$var2)
        }

        # check input file
        filein <- input$mainInput
        if (input$demoData == "lca-micros" | input$demoData == "ampk-tor") {
            filein <- 1
        }
        if (is.null(filein)) return()

        # get selected supertaxon name
        split <- strsplit(as.character(input$inSelect), "_")
        inSelect <- as.character(split[[1]][1])

        # get gene categories
        inputCatDt <- NULL
        if (input$colorByGroup == TRUE) {
            # get gene category
            geneCategoryFile <- input$geneCategory
            if (!is.null(geneCategoryFile)) {
                inputCatDt <-  read.table(
                    file = geneCategoryFile$datapath,
                    sep = "\t",
                    header = TRUE,
                    check.names = FALSE,
                    comment.char = "",
                    fill = TRUE
                )
                colnames(inputCatDt) <- c("geneID","group")
            } else {
                inputCatDt <- NULL
            }
        }

        # create data for heatmap plotting
        dataHeat <- filterProfileData(
            superTaxonData = dataSupertaxa(),
            refTaxon = inSelect,
            percentCutoff,
            coorthologCutoffMax,
            var1Cutoff,
            var2Cutoff,
            var1Relation = input$var1Relation,
            var2Relation = input$var2Relation,
            groupByCat = input$colorByGroup,
            catDt = inputCatDt
        )

        return(dataHeat)
    })

    # * clustered heatmap data -------------------------------------------------
    clusteredDataHeat <- reactive({
        dataHeat <- dataHeat()
        dat <- getProfiles()

        # do clustering based on distance matrix
        row.order <- hclust(getDistanceMatrixProfiles(),
                            method = input$clusterMethod)$order

        # re-order distance matrix accoring to clustering
        datNew <- dat[row.order, ] #col.order

        # return clustered gene ID list
        clusteredGeneIDs <- as.factor(row.names(datNew))

        # sort original data according to clusteredGeneIDs
        dataHeat$geneID <- factor(dataHeat$geneID,
                                  levels = clusteredGeneIDs)

        dataHeat <- dataHeat[!is.na(dataHeat$geneID),]
        return(dataHeat)
    })

    # =========================== MAIN PROFILE TAB =============================

    # * get total number of genes ----------------------------------------------
    output$totalGeneNumber.ui <- renderUI({
        geneList <- preData()
        geneList$geneID <- as.factor(geneList$geneID)
        out <- as.list(levels(geneList$geneID))

        listIn <- input$list
        if (!is.null(listIn)) {
            list <- read.table(file = listIn$datapath, header = FALSE)
            out <- as.list(list$V1)
        }

        if (length(out) > 0) {
            strong(paste0("Total number of genes:  ", length(out)))
        }
    })

    # * get list of taxa for highlighting --------------------------------------
    output$highlightTaxonUI <- renderUI({
        choice <- inputTaxonName()
        choice$fullName <- as.factor(choice$fullName)

        out <- as.list(levels(choice$fullName))
        out <- append("none", out)

        selectInput("taxonHighlight", "Select (super)taxon to highlight:",
                    out, selected = out[1])
    })

    # * get list of genes for highlighting -------------------------------------
    output$highlightGeneUI <- renderUI({
        geneList <- dataHeat()
        geneList$geneID <- as.factor(geneList$geneID)

        out <- as.list(levels(geneList$geneID))
        out <- append("none", out)

        selectInput("geneHighlight", "Highlight:", out, selected = out[1])
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
        if (input$autoUpdate == TRUE) {
            inputPara <- list(
                "xAxis" = input$xAxis,
                "var1ID" = input$var1ID,
                "var2ID"  = input$var2ID,
                "lowColorVar1" =  input$lowColorVar1,
                "highColorVar1" = input$highColorVar1,
                "lowColorVar2" = input$lowColorVar2,
                "highColorVar2" = input$highColorVar2,
                "paraColor" = input$paraColor,
                "xSize" = input$xSize,
                "ySize" = input$ySize,
                "legendSize" = input$legendSize,
                "mainLegend" = input$mainLegend,
                "dotZoom" = input$dotZoom,
                "xAngle" = input$xAngle,
                "guideline" = 1,
                "width" = input$width,
                "height" = input$height,
                "colorByGroup" = input$colorByGroup
            )
        } else {
            inputPara <- isolate(
                list(
                    "xAxis" = input$xAxis,
                    "var1ID" = input$var1ID,
                    "var2ID"  = input$var2ID,
                    "lowColorVar1" =  input$lowColorVar1,
                    "highColorVar1" = input$highColorVar1,
                    "lowColorVar2" = input$lowColorVar2,
                    "highColorVar2" = input$highColorVar2,
                    "paraColor" = input$paraColor,
                    "xSize" = input$xSize,
                    "ySize" = input$ySize,
                    "legendSize" = input$legendSize,
                    "mainLegend" = input$mainLegend,
                    "dotZoom" = input$dotZoom,
                    "xAngle" = input$xAngle,
                    "guideline" = 1,
                    "width" = input$width,
                    "height" = input$height,
                    "colorByGroup" = input$colorByGroup
                )
            )
        }
        return(inputPara)
    })

    # * render dot size to dotSizeInfo ---------------------------------------
    output$dotSizeInfo <- renderUI({
        if (v$doPlot == FALSE) return()

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
        typeProfile = reactive("mainProfile")
    )

    # ======================== CUSTOMIZED PROFILE TAB ==========================

    # * get list of all sequence IDs for customized profile -----
    output$geneIn <- renderUI({
        filein <- input$mainInput
        fileCustom <- input$customFile

        if (input$demoData == "lca-micros" | input$demoData == "ampk-tor") {
            filein <- 1
        }

        if (is.null(filein) & is.null(fileCustom)) {
            return(selectInput("inSeq", "", "all"))
        }
        if (v$doPlot == FALSE) {
            return(selectInput("inSeq", "", "all"))
        } else {
            # full list
            data <- as.data.frame(getDataFiltered())
            data$geneID <- as.character(data$geneID)
            data$geneID <- as.factor(data$geneID)
            outAll <- as.list(levels(data$geneID))
            outAll <- append("all", outAll)
            out <- list()
            if (input$addGeneAgeCustomProfile == TRUE) {
                out <- as.list(selectedgeneAge())
            } else if (input$addClusterCustomProfile == TRUE) {
                out <- as.list(brushedClusterGene())
            } else if (input$addCoreGeneCustomProfile == TRUE) {
                out <- as.list(coreGeneDf())
            } else if (input$addGCGenesCustomProfile == TRUE) {
                out <- as.list(geneListGC())
            } else {
                if (!is.null(fileCustom)) {
                    customList <- read.table(
                        file = fileCustom$datapath, header = FALSE
                    )

                    customList$V1 <- as.factor(customList$V1)
                    out <- as.list(levels(customList$V1))
                }
            }

            if (length(out) > 0) {
                createSelectGene("inSeq", out, out)
            } else {
                createSelectGene("inSeq", outAll, outAll[1])
            }
        }
    })

    # * render popup for selecting taxon rank and return list of subset taxa ---
    cusTaxaName <- callModule(
        selectTaxonRank,
        "selectTaxonRank",
        rankSelect = reactive(input$rankSelect),
        inputTaxonID = inputTaxonID
    )

    # * get list of all taxa for customized profile ----------------------------
    output$taxaIn <- renderUI({
        filein <- input$mainInput
        if (input$demoData == "lca-micros" | input$demoData == "ampk-tor") {
            filein <- 1
        }

        if (is.null(filein)) return(selectInput("inTaxa", "", "all"))
        if (v$doPlot == FALSE) return(selectInput("inTaxa", "", "all"))
        else {
            choice <- inputTaxonName()
            choice$fullName <- as.factor(choice$fullName)

            out <- as.list(levels(choice$fullName))
            out <- append("all", out)
            if (input$applyCusTaxa == TRUE) {
                out <- cusTaxaName()
                selectInput("inTaxa", "",
                            out,
                            selected = out,
                            multiple = TRUE,
                            selectize = FALSE)
            } else {
                selectInput("inTaxa", "",
                            out,
                            selected = out[1],
                            multiple = TRUE,
                            selectize = FALSE)
            }
        }
    })

    # * check if all genes and all species are selected ------------------------
    output$sameProfile <- reactive({
        if (v$doPlot == FALSE) return(FALSE)
        if (length(input$inSeq[1]) == 0) return(FALSE)
        else {
            if (input$inSeq[1] == "all" & input$inTaxa[1] == "all") {
                return(TRUE)
            }
        }
    })
    outputOptions(output, "sameProfile", suspendWhenHidden = FALSE)

    # * change value of autoUpdateSelected checkbox --------------------------
    observe({
        updateCheckboxInput(session,
                            "autoUpdateSelected", value = input$autoUpdate)
        shinyjs::disable("autoUpdateSelected")
        if (input$autoUpdate == TRUE) {
            shinyjs::disable("plotCustom")
        } else {
            shinyjs::enable("plotCustom")
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
        if (input$autoUpdateSelected == TRUE) {
            inputPara <- list(
                "xAxis" = input$xAxisSelected,
                "var1ID" = input$var1ID,
                "var2ID"  = input$var2ID,
                "lowColorVar1" =  input$lowColorVar1,
                "highColorVar1" = input$highColorVar1,
                "lowColorVar2" = input$lowColorVar2,
                "highColorVar2" = input$highColorVar2,
                "paraColor" = input$paraColor,
                "xSize" = input$xSizeSelect,
                "ySize" = input$ySizeSelect,
                "legendSize" = input$legendSizeSelect,
                "mainLegend" = input$selectedLegend,
                "dotZoom" = input$dotZoomSelect,
                "xAngle" = input$xAngleSelect,
                "guideline" = 0,
                "width" = input$selectedWidth,
                "height" = input$selectedHeight,
                "colorByGroup" = input$colorByGroup
            )
        } else {
            inputPara <- isolate(
                list(
                    "xAxis" = input$xAxisSelected,
                    "var1ID" = input$var1ID,
                    "var2ID"  = input$var2ID,
                    "lowColorVar1" =  input$lowColorVar1,
                    "highColorVar1" = input$highColorVar1,
                    "lowColorVar2" = input$lowColorVar2,
                    "highColorVar2" = input$highColorVar2,
                    "paraColor" = input$paraColor,
                    "xSize" = input$xSizeSelect,
                    "ySize" = input$ySizeSelect,
                    "legendSize" = input$legendSizeSelect,
                    "mainLegend" = input$selectedLegend,
                    "dotZoom" = input$dotZoomSelect,
                    "xAngle" = input$xAngleSelect,
                    "guideline" = 0,
                    "width" = input$selectedWidth,
                    "height" = input$selectedHeight,
                    "colorByGroup" = input$colorByGroup
                )
            )
        }
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
        typeProfile = reactive("customizedProfile")
    )

    # ============================== POINT INFO ================================

    # * get status of pointInfo for activating Detailed Plot button -----------
    output$pointInfoStatus <- reactive({
        if (input$tabs == "Main profile") {
            # info = groupID,orthoID,supertaxon,mVar1,%spec,var2
            info <- mainpointInfo()
        } else if (input$tabs == "Customized profile") {
            info <- selectedpointInfo()
        } else {
            info <- NULL
        }
        is.null(info)
    })
    outputOptions(output, "pointInfoStatus", suspendWhenHidden = FALSE)

    # * show info into "point's info" box --------------------------------------
    output$pointInfo <- renderText({
        # GET INFO BASED ON CURRENT TAB
        if (input$tabs == "Main profile") {
            # info = groupID,orthoID,supertaxon,mVar1,%spec,var2
            info <- mainpointInfo()
        } else if (input$tabs == "Customized profile") {
            info <- selectedpointInfo()
        } else {
            return()
        }

        if (is.null(info)) return()
        else {
            orthoID <- info[2]

            if (is.na(orthoID)) return()
            else {
                # if (orthoID=="NA") {orthoID <- info[2]}
                ## print output
                a <- toString(paste("Seed-ID:", info[1]))
                b <- toString(paste0("Hit-ID: ",
                                     orthoID,
                                     " (",
                                     substr(info[3], 6, nchar(info[3])),
                                     ")"))
                c <- ""
                if (input$var1ID != "") {
                    c <- toString(paste(input$var1AggregateBy,
                                        input$var1ID,
                                        ":",
                                        info[4]))
                }
                d <- ""
                if (input$var2ID != "") {
                    d <- toString(paste(input$var2AggregateBy,
                                        input$var2ID,
                                        ":",
                                        info[6]))
                }
                e <- toString(paste("% present taxa:", info[5]))
                paste(a, b, c, d, e, sep = "\n")
            }
        }
    })

    # ============================= DETAILED PLOT ==============================

    # * data for detailed plot -------------------------------------------------
    detailPlotDt <- reactive({
        if (v$doPlot == FALSE) return()

        # GET INFO BASED ON CURRENT TAB
        if (input$tabs == "Main profile") {
            # info = groupID,orthoID,supertaxon,mVar1,%spec,var2
            info <- mainpointInfo()
        } else if (input$tabs == "Customized profile") {
            info <- selectedpointInfo()
        }

        if (is.null(info)) return()
        else {
            ### get info for present taxa in selected supertaxon (1)
            plotTaxon <- info[3]
            plotGeneID <- info[1]
            fullDf <- getDataFiltered()
            selDf <- as.data.frame(fullDf[fullDf$geneID == plotGeneID
                                          & fullDf$supertaxon == plotTaxon, ])

            ### get all taxa of this supertaxon (2)
            allTaxaDf <- sortedtaxaList()
            allTaxaDf <- allTaxaDf[allTaxaDf$supertaxon == plotTaxon, ]
            allTaxaDf <- subset(allTaxaDf, select = c("abbrName", "fullName"))

            ### merge (1) and (2) together
            joinedDf <- merge(selDf, allTaxaDf,
                              by = c("abbrName"),
                              all.y = TRUE)
            joinedDf <- subset(joinedDf,
                               select = c("abbrName",
                                          "fullName.y",
                                          "geneID",
                                          "orthoID",
                                          "var1",
                                          "var2"))
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
        }
    })

    # * render detailed plot ---------------------------------------------------

    pointInfoDetail <- callModule(
        createDetailedPlot, "detailedPlot",
        data = detailPlotDt,
        var1ID = reactive(input$var1ID),
        var2ID = reactive(input$var2ID),
        detailedText = reactive(input$detailedText),
        detailedHeight = reactive(input$detailedHeight)
    )

    # * render FASTA sequence --------------------------------------------------
    output$fasta <- renderText({
        if (v$doPlot == FALSE) return()

        info <- pointInfoDetail() # info = seedID, orthoID, var1

        if (is.null(info)) return()
        else {
            data <- getDataFiltered()

            seqID <- toString(info[2])
            groupID <- toString(info[1])
            ncbiID <- gsub("ncbi", "", toString(info[5]))

            if (input$demoData == "none") {
                filein <- input$mainInput
                inputType <- checkInputValidity(
                    filein$datapath
                )
            } else {
                inputType <- "demo"
            }

            if (inputType == "oma") {
                fastaOut <- getSelectedFastaOma(finalOmaDf(), seqID)
            } else {
                seqDf <- data.frame("geneID" = groupID,
                                    "orthoID" = seqID,
                                    "ncbiID" = ncbiID)

                filein <- input$mainInput
                fastain <- input$concatFasta
                fastaOut <- getFastaSeqs(
                    seqDf, filein$datapath, input$demoData,
                    input$inputType, fastain$datapath,
                    input$path,
                    input$dirFormat,
                    input$fileExt,
                    input$idFormat
                )
            }

            return(paste(fastaOut[1]))
        }
    })

    # ======================== FEATURE ARCHITECTURE PLOT =======================

    # * get domain file/path ---------------------------------------------------
    checkDomainFile <- reactive({
        # click info
        info <- pointInfoDetail() # info = seedID, orthoID, var1
        group <- as.character(info[1])
        ortho <- as.character(info[2])
        var1 <- as.character(info[3])

        if (is.null(info)) {
            updateButton(session, "doDomainPlot", disabled = TRUE)
            return("noSelectHit")
        } else {
            if (input$demoData == "lca-micros" |
                input$demoData == "ampk-tor") {
                updateButton(session, "doDomainPlot", disabled = FALSE)
            } else {
                if (input$annoLocation == "from file") {
                    inputDomain <- input$fileDomainInput
                    if (is.null(inputDomain)) {
                        updateButton(session, "doDomainPlot",
                                     disabled = TRUE)
                        return("noFileInput")
                    } else {
                        updateButton(session, "doDomainPlot",
                                     disabled = FALSE)
                    }
                } else {
                    domainDf <- parseDomainInput(
                        info[1],
                        input$domainPath,
                        "folder"
                    )
                    if (length(domainDf) == 1) {
                        if (domainDf == "noSelectHit" |
                            domainDf == "noFileInFolder") {
                            updateButton(session, "doDomainPlot",
                                         disabled = TRUE)
                            return(domainDf)
                        } else {
                            updateButton(session, "doDomainPlot",
                                         disabled = FALSE)
                        }
                    } else {
                        updateButton(session, "doDomainPlot",
                                     disabled = FALSE)
                    }
                }
            }
        }

        return("correct")
    })

    # * check domain file ------------------------------------------------------
    output$checkDomainFiles <- renderUI({
        fileDomain <- checkDomainFile()
        if (fileDomain == "noFileInput") {
            em("Domain file not provided!!")
        } else if (fileDomain == "noFileInFolder") {
            msg <- paste0(
                "<p><em>Domain file not found!! </em></p>
        <p><em>Please make sure that file name has to be in this format:
        <strong>&lt;seedID&gt;.extension</strong>, where extension is limited
        to <strong>txt</strong>, <strong>csv</strong>, <strong>list</strong>,
        <strong>domains</strong> or <strong>architecture</strong>.
        </em></p>"
            )
            HTML(msg)
        } else if (fileDomain == "noSelectHit") {
            em("Please select one ortholog sequence!!")
        }
    })

    # * render domain plot -----------------------------------------------------
    observeEvent(input$doDomainPlot, {
        callModule(
            createArchitecturePlot, "archiPlot",
            pointInfo = pointInfoDetail,
            domainInfo = getDomainInformation,
            labelArchiSize = reactive(input$labelArchiSize),
            titleArchiSize = reactive(input$titleArchiSize),
            archiHeight = reactive(input$archiHeight),
            archiWidth = reactive(input$archiWidth)
        )
    })

    # ======================== FILTERED DATA DOWNLOADING =======================

    # * for main profile =======================================================
    mainFastaDownload <- reactive({
        if (input$demoData == "none") {
            filein <- input$mainInput
            inputType <- checkInputValidity(filein$datapath)
        } else {
            inputType <- "demo"
        }

        if (inputType == "oma") {
            allOmaDf <- finalOmaDf()
            filteredDownloadDf <- as.data.frame(downloadData())
            filteredOmaDf <-
                subset(allOmaDf,
                       allOmaDf$orthoID %in% filteredDownloadDf$orthoID &
                           allOmaDf$seed %in% filteredDownloadDf$geneID)
            mainFastaOut <- getAllFastaOma(filteredOmaDf)
        } else {
            mainFastaOut <- getFastaSeqs(
                as.data.frame(downloadData()),
                input$mainInput, input$demoData,
                input$inputType, input$concatFasta,
                input$path,
                input$dirFormat,
                input$fileExt,
                input$idFormat
            )
        }

        return(mainFastaOut)
    })

    downloadData <- callModule(
        downloadFilteredMain,
        "filteredMainDownload",
        data = getDataFiltered,
        fasta = mainFastaDownload,
        var1ID = reactive(input$var1ID),
        var2ID = reactive(input$var2ID),
        var1 = reactive(input$var1),
        var2 = reactive(input$var2),
        percent = reactive(input$percent)
    )

    # * for customized profile =================================================
    customizedFastaDownload <- reactive({
        if (input$demoData == "none") {
            filein <- input$mainInput
            inputType <- checkInputValidity(filein$datapath)
        } else {
            inputType <- "demo"
        }

        if (inputType == "oma") {
            allOmaDf <- finalOmaDf()
            filteredDownloadDf <- as.data.frame(downloadCustomData())
            filteredOmaDf <-
                subset(allOmaDf,
                       allOmaDf$orthoID %in% filteredDownloadDf$orthoID &
                           allOmaDf$seed %in% filteredDownloadDf$geneID)
            fastaOutDf <- getAllFastaOma(filteredOmaDf)
        } else {
            fastaOutDf <- getFastaSeqs(
                as.data.frame(downloadCustomData()),
                input$mainInput, input$demoData,
                input$inputType, input$concatFasta,
                input$path,
                input$dirFormat,
                input$fileExt,
                input$idFormat
            )
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

    # ============================ ANALYSIS FUNCTIONS ==========================

    # * PROFILE CLUSTERING =====================================================
    # ** description for profile clustering function ---------------------------
    observe({
        desc = paste("Cluster genes according to the distance of their
                 phylogenetic profiles.")

        if (input$tabs == "Profiles clustering") {
            createAlert(session, "descClusteringUI", "descClustering",
                        content = desc, append = FALSE)
        }
    })

    # ** check if genes are added anywhere else to the customized profile ------
    observe({
        if (input$addGeneAgeCustomProfile == TRUE
            | input$addCoreGeneCustomProfile == TRUE
            | input$addGCGenesCustomProfile == TRUE) {
            shinyjs::disable("addClusterCustomProfile")
        } else {
            shinyjs::enable("addClusterCustomProfile")
        }
    })

    output$addClusterCustomProfileCheck.ui <- renderUI({
        if (input$addGeneAgeCustomProfile == TRUE
            | input$addCoreGeneCustomProfile == TRUE |
            input$addGCGenesCustomProfile == TRUE ) {
            HTML('<p><em>(Uncheck "Add to Customized profile" check box in
           <strong>Gene age estimation</strong> or
           <strong>Core genes finding</strong> or
           <strong>Group comparison</strong>
           &nbsp;to enable this function)</em></p>')
        }
    })

    # ** List of possible profile types ----------------------------------------
    output$selectProfileType <- renderUI({
        variable1 <- paste0("profile using ", input$var1ID)
        if (input$var2ID != "") {
            variable2 <- paste0("profile using ", input$var2ID)
            radioButtons(
                "profileType",
                label = h5("Select the profile type"),
                choiceNames = list(
                    "binary profile",
                    variable1,
                    variable2),
                choiceValues = list(
                    "binary", "var1", "var2"
                ),
                selected = "binary",
                inline = FALSE)
        }
        else {
            radioButtons(
                "profileType",
                label = h5("Select the profile type"),
                choiceNames = list(
                    "binary profile",
                    variable1),
                choiceValues = list(
                    "binary", "var1"
                ),
                selected = "binary",
                inline = FALSE)
        }
    })

    # ** List of possible distance methods -------------------------------------
    output$selectDistMethod <- renderUI({
        if (is.null(input$profileType)) return()

        if (input$profileType == "binary") {
            selectInput(
                "distMethod",
                label = h5("Distance measure method:"),
                choices = list("euclidean" = "euclidean",
                               "maximum" = "maximum",
                               "manhattan" = "manhattan",
                               "canberra" = "canberra",
                               "binary" = "binary",
                               "pearson correlation coefficient" = "pearson",
                               "mutual information" = "mutualInformation",
                               "distance correlation" = "distanceCorrelation"
                ),
                selected = "euclidean"
            )
        } else {
            selectInput(
                "distMethod",
                label = h5("Distance measure method:"),
                choices = list("mutual information" = "mutualInformation",
                               "distance correlation" = "distanceCorrelation"
                ),
                selected = "mutualInformation"
            )
        }
    })

    # ** Distance matrix -------------------------------------------------------
    getDistanceMatrixProfiles <- reactive({
        if (is.null(input$distMethod)) return()
        profiles <- getProfiles()
        distanceMatrix <- getDistanceMatrix(
            profiles, input$distMethod
        )
        return(distanceMatrix)
    })

    # ** Phylogenetic profiles -------------------------------------------------
    getProfiles <- reactive({
        dataHeat <- dataHeat()
        req(dataHeat)
        if (is.null(input$distMethod)) return()

        profiles <- getDataClustering(dataHeat,
                                        input$profileType,
                                        input$var1AggregateBy,
                                        input$var2AggregateBy)
        return(profiles)
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
        desc = paste("Plot the distributions of the values incurred by the
                 integrated information layers.")

        if (input$tabs == "Distribution analysis") {
            createAlert(session, "descDistributionUI", "descDistribution",
                        content = desc, append = FALSE)
        }
    })

    # ** list of available variables for distribution plot ---------------------
    output$selected.distribution <- renderUI({
        if (nchar(input$var1ID) == 0 & nchar(input$var2ID) == 0) {
            varList <- "% present taxa"
        } else if (nchar(input$var1ID) == 0 & nchar(input$var2ID) > 0) {
            varList <- as.list(c(input$var2ID, "% present taxa"))
        } else if (nchar(input$var1ID) > 0 & nchar(input$var2ID) == 0) {
            varList <- as.list(c(input$var1ID,
                                 "% present taxa"))
        } else {
            varList <- as.list(c(input$var1ID,
                                 input$var2ID,
                                 "% present taxa"))
        }

        selectInput("selectedDist",
                    "Choose variable to plot:",
                    varList,
                    varList[1])
    })

    # ** var1 / var2 distribution data -----------------------------------------
    distributionDf <- reactive({
        if (v$doPlot == FALSE) return()
        dataOrig <- getMainInput()

        splitDt <- createVariableDistributionData(
            dataOrig,
            input$var1[1],
            input$var1[2],
            input$var2[1],
            input$var1[2]
        )

        # filter data base on customized plot (if chosen)
        if (input$dataset.distribution == "Customized data") {
            splitDt <- createVariableDistributionDataSubset(
                getDataFiltered(),
                splitDt,
                input$inSeq,
                input$inTaxa
            )
        }

        # return dt
        return(splitDt)
    })

    # ** calculate % present species in supertaxa ------------------------------
    presSpecAllDt <- reactive({
        # open main input file
        mdData <- getMainInput()
        # get list of sorted input taxa
        sortedTaxa <- sortedtaxaList()
        # calculate % present species for the whole data
        # rankName <- substr(input$rankSelect, 4, nchar(input$rankSelect))
        rankName <- input$rankSelect
        finalPresSpecDt <- createPercentageDistributionData(mdData,
                                                               rankName)

        return(finalPresSpecDt)
    })

    # ** render distribution plots ---------------------------------------------
    observe({
        if (v$doPlot == FALSE) return()

        if (is.null(input$selectedDist)) {
            return()
        } else {
            if (input$selectedDist == "% present taxa") {
                callModule(
                    analyzeDistribution, "distPlot",
                    data = presSpecAllDt,
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
            createAlert(session, "descGeneAgeUI", "descGeneAge",
                        title = "", content = desc, append = FALSE)
        }
    })

    # ** check if genes are added anywhere else to the customized profile ------
    observe({
        if (input$addClusterCustomProfile == TRUE
            | input$addCoreGeneCustomProfile == TRUE
            | input$addGCGenesCustomProfile == TRUE ) {
            shinyjs::disable("addGeneAgeCustomProfile")
        } else {
            shinyjs::enable("addGeneAgeCustomProfile")
        }
    })

    output$addGeneAgeCustomProfileCheck.ui <- renderUI({
        if (input$addClusterCustomProfile == TRUE
            | input$addCoreGeneCustomProfile == TRUE
            | input$addGCGenesCustomProfile == TRUE) {
            HTML('<p><em>(Uncheck "Add to Customized profile" check box in
           <strong>Profile clustering</strong> or
           <strong>Core genes finding</strong> or
           <strong>Group comparison</strong>
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
        if (v$doPlot == FALSE) return()
        geneAgeDf <- estimateGeneAge(
            getDataFiltered(),
            toString(input$rankSelect),
            input$inSelect,
            input$var1, input$var2, input$percent
        )
        return(geneAgeDf)
    })

    # ** render age distribution plot ------------------------------------------
    selectedgeneAge <- callModule(
        plotGeneAge, "geneAge",
        data = geneAgeDf,
        geneAgeWidth = reactive(input$geneAgeWidth),
        geneAgeHeight = reactive(input$geneAgeHeight),
        geneAgeText = reactive(input$geneAgeText)
    )

    # * CORE GENES IDENTIFICATION ==============================================
    # ** description for core gene identification function ---------------------
    observe({
        desc = paste("IDENTIFY GENES THAT ARE SHARED AMONG SELECTED TAXA.",
                     "You can set the minimal taxa that should be taken into
                 account by using the \"Core taxa coverage\" cutoff.",
                     "If you are working with a taxonomy level (e.g. Family)
                 that is higher than the one in the input profile (e.g.
                 Species), you can also identify a minimal fragtion of species
                 that need to have an ortholog in each supertaxon with
                 \"% of present taxa\" cutoff.")

        if (input$tabs == "Core gene identification") {
            createAlert(session, "descCoreGeneUI", "descCoreGene",
                        title = "", content = desc, append = FALSE)
        }
    })

    # ** render list of available taxa -----------------------------------------
    output$taxaListCore.ui <- renderUI({
        filein <- input$mainInput
        if (input$demoData == "lca-micros" | input$demoData == "ampk-tor") {
            filein <- 1
        }
        if (is.null(filein)) {
            return(selectInput("inTaxa",
                               "Select taxa of interest:",
                               "none"))
        }
        if (v$doPlot == FALSE) {
            return(selectInput("inTaxa",
                               "Select taxa of interest:",
                               "none"))
        } else {
            choice <- inputTaxonName()
            choice$fullName <- as.factor(choice$fullName)

            out <- as.list(levels(choice$fullName))
            out <- append("none", out)

            if (input$applyCoreTaxa == TRUE) {
                out <- coreTaxaName()
                selectInput("taxaCore",
                            "Select taxa of interest:",
                            out,
                            selected = out,
                            multiple = TRUE)
            } else {
                selectInput("taxaCore",
                            "Select taxa of interest:",
                            out,
                            selected = out[1],
                            multiple = TRUE)
            }
        }
    })

    # ** render popup for selecting group of taxa to find core genes -----------
    coreTaxaName <- callModule(
        selectTaxonRank,
        "selectTaxonRankCore",
        rankSelect = reactive(input$rankSelect),
        inputTaxonID = inputTaxonID
    )

    # ** check if genes are added anywhere else to the customized profile ------
    observe({
        if (input$addClusterCustomProfile == TRUE
            | input$addGeneAgeCustomProfile == TRUE
            | input$addGCGenesCustomProfile == TRUE) {
            shinyjs::disable("addCoreGeneCustomProfile")
        } else {
            shinyjs::enable("addCoreGeneCustomProfile")
        }
    })

    output$addCoreGeneCustomProfileCheck.ui <- renderUI({
        if (input$addClusterCustomProfile == TRUE
            | input$addGeneAgeCustomProfile == TRUE
            | input$addGCGenesCustomProfile == TRUE) {
            HTML('<p><em>(Uncheck "Add to Customized profile" check box in
           <strong>Profiles clustering</strong> or
           <strong>Gene age estimating</strong> or
           <strong>Group Comparioson</strong>
           &nbsp;to enable this function)</em></p>')
        }
    })

    # ** render table contains list of core genes ------------------------------
    coreGeneDf <- callModule(
        identifyCoreGene,
        "coreGene",
        filteredData = getDataFiltered,
        rankSelect = reactive(input$rankSelect),
        taxaCore = reactive(input$taxaCore),
        percentCore = reactive(input$percentCore),
        var1Cutoff = reactive(input$var1Core),
        var2Cutoff = reactive(input$var2Core),
        coreCoverage = reactive(input$coreCoverage)
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
        desc = paste(desc, "between two taxon groups, an in- and
                 an out-group. You can define the in-group below and all taxa
                 not included in this are used as the out-group. The value
                 distributions of the variables are then compared using
                 statistical tests (Kolmogorov-Smirnov and
                 Wilcoxon-Mann-Whitney) using the specified significant level.
                 Genes that have a significantly different distribution are
                 shown in the candidate gene list below.")

        if (input$tabs == "Group comparison") {
            createAlert(session, "descGCUI", "descGC", title = "",
                        content = desc, append = FALSE)
        }
    })

    # ** list of all available genes -------------------------------------------
    output$listGenesGC <- renderUI({
        filein <- input$mainInput

        fileGC <- input$gcFile
        if (input$demoData == "lca-micros" | input$demoData == "ampk-tor") {
            filein <- 1
        }

        if (v$doPlot == FALSE) {
            return(selectInput(
                "listSelectedGenesGC", "Select sequence(s):", "none"
            ))
        }

        if (is.null(filein) & is.null(fileGC)) {
            return(selectInput(
                "listSelectedGenesGC", "Select sequence(s):", "none"
            ))
        } else {
            # full list
            data <- as.data.frame(getDataFiltered())
            data$geneID <- as.factor(data$geneID)
            outAll <- as.list(levels(data$geneID))
            outAll <- append("all", outAll)

            if (is.null(fileGC)) {
                selectInput("listSelectedGenesGC", "Select sequence(s):",
                            outAll,
                            selected = outAll[1],
                            multiple = TRUE,
                            selectize = FALSE)
            } else {
                listGC <- read.table(file = fileGC$datapath, header = FALSE)
                listGC$V1 <- as.factor(listGC$V1)
                out <- as.list(levels(listGC$V1))
                selectInput("listSelectedGenesGC", "Select sequence(s):",
                            out,
                            selected = NULL,
                            multiple = FALSE,
                            selectize = FALSE)
            }
        }
    })

    # ** popup for selecting taxon rank and return list of belonging taxa ------
    gcTaxaName <- callModule(
        selectTaxonRank,
        "selectTaxonRankGC",
        rankSelect = reactive(input$rankSelect),
        inputTaxonID = inputTaxonID
    )

    # ** list of available taxa (for selecting as inGroup) --------------------
    output$taxaListGC <- renderUI({
        filein <- input$mainInput
        if (input$demoData == "lca-micros" | input$demoData == "ampk-tor") {
            filein <- 1
        }
        if (is.null(filein)) {
            return(selectInput("inTaxa", "Select inGroup taxa:", "none"))
        }
        if (v$doPlot == FALSE) {
            return(selectInput("inTaxa", "Select inGroup taxa:", "none"))
        } else {
            choice <- inputTaxonName()
            choice$fullName <- as.factor(choice$fullName)

            out <- as.list(levels(choice$fullName))

            #' when the taxonomy rank was changed ------------------------------
            if (input$applyTaxaGC == TRUE) {
                out <- gcTaxaName()
                selectInput("selectedInGroupGC", "Select inGroup taxa:",
                            out,
                            selected = out,
                            multiple = TRUE,
                            selectize = FALSE)
            }
            #' when the taxonomy is the same as the initially chosen one -------
            else {
                #' check for the rank of the rank in the input
                ranks <- getTaxonomyRanks()
                pos <- which(ranks == input$rankSelect) # position in the list
                higherRank <- ranks[pos + 1] # take the next higher rank
                higherRankName <- as.character(higherRank[1])

                nameList <- getNameList() # get the taxon names
                dt <- getTaxonomyMatrix(FALSE, NULL) # get the taxa

                #' get the info for the reference protein from the namelist
                reference <- subset(
                    nameList, nameList$fullName == input$inSelect
                )

                #' get the id for every rank for the reference protein
                rankName <- input$rankSelect
                referenceDt <- dt[dt[, rankName] == reference$ncbiID, ]

                #' save the next higher rank
                referenceHigherRank <- referenceDt[higherRankName]
                referenceHigherRank <-
                    referenceHigherRank[!duplicated(referenceHigherRank), ]

                #' get all the taxa with the same id in the next higher rank
                selectedTaxaDt <-
                    subset(
                        dt, dt[, higherRankName] %in% referenceHigherRank
                    )
                selectedTaxaDt <-
                    selectedTaxaDt[!duplicated(selectedTaxaDt[rankName]), ]

                #' get list with all ids with referenceHigherRank as parent
                selectedTaxaIDs <- selectedTaxaDt[rankName]
                if (length(selectedTaxaIDs[[1]]) >= 1) {
                    selectedTaxaIDs <- selectedTaxaIDs[[1]]
                }

                selectedTaxa <- subset(nameList, nameList$rank == rankName)
                selectedTaxa <-
                    subset(
                        selectedTaxa,
                        selectedTaxa$ncbiID %in% selectedTaxaIDs
                    )

                defaultSelect <- selectedTaxa$fullName

                selectInput("selectedInGroupGC", "Select inGroup taxa:",
                            out,
                            selected = defaultSelect,
                            multiple = TRUE,
                            selectize = FALSE)
            }
        }
    })

    # ** buttons to choose the variable ----------------------------------------
    output$variableButtonGC <- renderUI({
        radioButtons(
            inputId = "varNameGC",
            label = "Select variable(s) to compare:",
            choices = list(input$var1ID, input$var2ID, "Both"),
            selected = input$var1ID,
            inline = FALSE
        )
    })

    # ** slider to set significance level --------------------------------------
    output$significance.ui <- renderUI({
        msg <- paste0(
            "P-value cut-off of the statistic test"
        )
        popify(
            sliderInput(
                "significance",
                paste("Significance level:"),
                min = 0,
                max = 1,
                step = 0.005,
                value = c(0.05),
                width = 200
            ),
            "",
            msg
        )
    })

    # ** reset plots config ----------------------------------------------------
    observeEvent(input$resetConfigGC, {
        shinyjs::reset("xSizeGC")
        shinyjs::reset("ySizeGC")
        shinyjs::reset("angleGC")
        shinyjs::reset("legendSizeGC")
    })

    observeEvent(input$applyConfigGC, {
        toggleModal(session, "gcPlotConfigBs", toggle = "close")
    })

    # ** check if genes are added anywhere else to the customized profile ------
    observe({
        if (input$addGeneAgeCustomProfile == TRUE |
            input$addCoreGeneCustomProfile == TRUE |
            input$addClusterCustomProfile == TRUE) {
            shinyjs::disable("addGCGenesCustomProfile")
        } else {
            shinyjs::enable("addGCGenesCustomProfile")
        }
    })

    output$addGCCustomProfileCheck <- renderUI({
        if (input$addGeneAgeCustomProfile == TRUE |
            input$addCoreGeneCustomProfile == TRUE |
            input$addClusterCustomProfile == TRUE) {
            HTML('<p><em>(Uncheck "Add to Customized profile" check box in
           <strong>Gene age estimation</strong> or
           <strong>Profile clustering</strong> or
           <strong>Core genes finding</strong>&nbsp;to enable this function)
           </em></p>')
        }
    })

    # ** parameters for the plots in Group Comparison --------------------------
    getParameterInputGC <- reactive({
        inputData <- list(
            "showPValue" = input$showPValue,
            "highlightSignificant" = input$highlightSignificant,
            "significance" = input$significance,
            "var1ID" = input$var1ID,
            "var2ID" = input$var2ID,
            "xSizeGC" = input$xSizeGC,
            "ySizeGC" = input$ySizeGC,
            "interestingFeatures" = input$interestingFeatures,
            "angleGC" = input$angleGC,
            "legendGC" = input$legendGC,
            "legendSizeGC" = input$legendSizeGC,
            "pValuesSize" = input$pValuesSizeGC
        )
    })

    # ** render plots for group comparison -------------------------------------
    geneListGC <- callModule(
        groupComparison, "groupComparison",
        selectedInGroup = reactive(input$selectedInGroupGC),
        selectedGenesList = reactive(input$listSelectedGenesGC),
        mainRank = reactive(input$rankSelect),
        selectedVariable = reactive(input$varNameGC),
        useCommonAncestor = reactive(input$useCommonAncestor),
        referenceTaxon = reactive(input$inSelect),
        ncbiIDList = inputTaxonID,
        filteredData = getDataFiltered,
        rightFormatFeatures = reactive(input$rightFormatFeatures),
        domainInformation = getDomainInformation,
        plot = reactive(input$plotGC),
        parameter = getParameterInputGC,
        selectedPoint = reactive(input$showPointGC)
    )
})
