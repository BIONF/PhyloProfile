#' Protein domain architecture plot
#' @param pointInfo() info of clicked point
#' (from reactive fn "pointInfoDetail")
#' @param domainInfo() domain information
#' (from reactive fn "getDomainInformation")
#' @param labelArchiSize lable size (from input$labelArchiSize)
#' @param titleArchiSize title size (from input$titleArchiSize)
#' @param archiHeight plot height (from input$archiHeight)
#' @param archiWidth plot width (from input$archiWidth)
#' @param currentNCBIinfo dataframe of the pre-processed NCBI taxonomy data
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

source("R/functions.R")

createArchitecturePlotUI <- function(id) {
    ns <- NS(id)
    tagList(
        fluidRow(
            column(
                3,
                radioButtons(
                    ns("resolveOverlap"),
                    "Merge non-overlapped features",
                    choices = c("Yes","No"), selected = "Yes",
                    inline = TRUE
                ),
                checkboxGroupInput(
                    ns("namePostion"),
                    "Display feature names",
                    choices = c(
                        "On the plot" = "plot",
                        "As a legend" = "legend",
                        "On the y-axis" = "axis"
                    ),
                    selected = c("plot","axis")
                )
            ),
            column(
                3,
                selectInput(
                    ns("feature"),
                    "Exclude features",
                    choices = c(
                        "flps","seg","coils","signalp","tmhmm",
                        "smart","pfam",
                        "without E-value" = "noEvalue",
                        "without Bit-score" = "noBitscore"
                    ),
                    multiple = TRUE
                ),
                checkboxInput(
                    ns("featureOpt"), "Other feature options", value = FALSE
                ),
                checkboxInput(
                    ns("plotConfig"), "Plot configuration", value = FALSE
                )
            ),
            column(
                3,
                strong("Show information"),
                checkboxGroupInput(
                    ns("showWeight"),
                    "",
                    choices = "Weight"
                ),
                checkboxGroupInput(
                    ns("showScore"),
                    "",
                    choices = c(
                        "E-value", "Bit-score"
                    )
                ),
                uiOutput(ns("filterEvalue.ui")),
                uiOutput(ns("filterBitscore.ui"))
            ),
            column(
                3,
                selectInput(
                    ns("linearizationBy"),
                    "Linearizing architecture using",
                    choices = c(
                        "None" = "none",
                        "Best E-value" = "evalue",
                        "Best Bit-score" = "bitscore",
                        "Paths" = "path"
                    )
                )
            )
        ),
        br(),
        fluidRow(
            conditionalPanel(
                condition = {sprintf("input['%s'] == 1", ns("featureOpt"))},
                column(
                    3,
                    radioButtons(
                        ns("nameType"),"Type of feature names", inline = TRUE,
                        choices = c("Labels","Texts"), selected = "Labels"
                    )

                ),
                column(
                    4,
                    conditionalPanel(
                        condition = {
                            sprintf("input['%s'] == 'Labels'", ns("nameType"))
                        },
                        radioButtons(
                            ns("labelPos"),"Label position", inline = TRUE,
                            choices = c("Above","Inside","Below"),
                            selected = "Above"
                        )
                    ),
                    conditionalPanel(
                        condition = {
                            sprintf("input['%s'] == 'Texts'", ns("nameType"))
                        },
                        colourpicker::colourInput(
                            ns("nameColor"),
                            "Feature name color",
                            value = "#000000"
                        )
                    )
                ),
                column(
                    5,
                    selectInput(
                        ns("excludeNames"),
                        "Exclude feature names of",
                        choices = c(
                            "flps","seg","coils","signalp","tmhmm",
                            "smart","pfam"
                        ),
                        selected = c("tmhmm","signalp","seg","coils"),
                        multiple = TRUE
                    )
                ),
                column(
                    3,
                    radioButtons(
                        ns("featureClassSort"),
                        "Sort feature classes by shared features",
                        choices = c("Yes","No"), selected = "Yes",
                        inline = TRUE
                    )
                ),
                column(
                    5,
                    conditionalPanel(
                        condition = {
                            sprintf("input['%s']=='No'",ns("featureClassSort"))
                        },
                        selectInput(
                            ns("featureClassOrder"),
                            "Feature class order",
                            choices = c(
                                "pfam", "smart", "tmhmm", "coils", "signalp",
                                "seg", "flps"
                            ),
                            selected = c(
                                "pfam", "smart", "tmhmm", "coils", "signalp",
                                "seg", "flps"
                            ),
                            multiple = TRUE
                        )
                    )
                )
            )
        ),
        br(),
        fluidRow(
            conditionalPanel(
                condition = {sprintf("input['%s'] == 1", ns("plotConfig"))},
                column(
                    3,
                    createPlotSize(ns("archiHeight"),"Plot height(px)",400,200),
                    createPlotSize(ns("archiWidth"),"Plot width (px)", 800, 200)
                ),
                column(
                    3,
                    createTextSize(
                        ns("titleArchiSize"), "Title/Seq ID size (px)", 14, 200
                    ),
                    createTextSize(
                        ns("labelArchiSize"), "Axis label size(px)", 12, 200
                    ),
                    createTextSize(
                        ns("legendArchiSize"), "Legend font size(px)", 12, 200
                    )
                ),
                column(
                    6,
                    column(
                        6,
                        createTextSize(
                            ns("segmentSize"),"Feature segment size (mm)",5,200
                        )
                    ),
                    column(
                        6,
                        createTextSize(
                            ns("nameSize"), "Feature ID size (mm)", 3, 200
                        )
                    ),
                    column(
                        12,
                        sliderInput(
                            ns("firstDist"),
                            "Distance between plot title and the 1st feature",
                            min = 0, max = 5, value = 0.5, step = 0.1, width=400
                        )
                    )
                ),
                column(
                    6,
                    radioButtons(
                        ns("colorType"),"Color feature instances", inline=TRUE,
                        choices = c("Shared","Unique","All","Feature type"),
                        selected = "All"
                    ),
                    checkboxInput(
                        ns("ignoreInstanceNo"), "Ignore number of instances",
                        value = FALSE
                    ),
                ),
                column(
                    6,
                    selectInput(
                        ns("colorPallete"),
                        "Color pallete",
                        choices = c(
                            "Paired", "Set1", "Set2", "Set3", "Accent", "Dark2"
                        ),
                        selected = "Paired"
                    )
                )
            )
        ),
        hr(),
        uiOutput(ns("archiPlot.ui")),
        verbatimTextOutput(ns("hover_info")),
        br(),
        downloadButton(ns("archiDownload"), "Download plot", class = "butDL"),
        hr(),
        tableOutput(ns("linkTable")),
        checkboxInput(
            ns("showDomainTable"), "Show detailed feature table", value = FALSE
        ),
        conditionalPanel(
            condition = {sprintf("input['%s'] == 1", ns("showDomainTable"))},
            DT::dataTableOutput(ns("domainTable"))
        )
    )
}

createArchitecturePlot <- function(
    input, output, session, pointInfo, domainInfo, currentNCBIinfo, font="sans"
){
    # * update excludeNames if no feature type on the y-axis ===================
    observe({
        req(input$namePostion)
        if (
            !("axis" %in% input$namePostion) & !("legend" %in% input$namePostion)
        ) {
            updateSelectInput(
                session, "excludeNames",
                "Exclude feature names of",
                choices = c("seg","coils","signalp","tmhmm","smart","pfam")
            )
        } else if (
            "axis" %in% input$namePostion | "legend" %in% input$namePostion
        ) {
            updateSelectInput(
                session, "excludeNames",
                "Exclude feature names of",
                choices = c(
                    "flps","seg","coils","signalp","tmhmm","smart","pfam"
                ),
                selected = c("tmhmm","signalp","seg","coils")
            )
        }
    })

    # * render e-value / bitscore filter =======================================
    output$filterEvalue.ui <- renderUI({
        req(domainInfo())
        df <- domainInfo()
        maxEvalue = 1
        if ("evalue" %in% colnames(df))
            maxEvalue <- format(
                max(df$evalue[!is.na(df$evalue)]), scientific = TRUE, digits = 2
            )
        if ("E-value" %in% input$showScore) {
            numericInput(
                session$ns("minEvalue"), "Filter E-value:",
                min = 0,
                max = maxEvalue,
                value = format(0.00001, scientific = TRUE, digits = 2)
            )
        }
    })

    output$filterBitscore.ui <- renderUI({
        req(domainInfo())
        df <- domainInfo()
        if ("Bit-score" %in% input$showScore) {
            numericInput(
                session$ns("minBitscore"), "Filter Bit-score:",
                min = min(df$bitscore[!is.na(df$bitscore)]),
                max = 9999,
                value = min(df$bitscore[!is.na(df$bitscore)])
            )
        }
    })

    # * filter domain features =================================================
    filterDomainDf <- reactive({
        outDf <- domainInfo()
        if (is.null(outDf)) stop("Domain info is NULL!")
        if (is.null(nrow(outDf))) stop("Domain info is NULL!")

        # get subset domainDf containing only domains for seed and ortholog
        group <- as.character(pointInfo()[1])
        ortho <- as.character(pointInfo()[2])
        # get sub dataframe based on selected groupID and orthoID
        group <- gsub("\\|", ":", group)
        ortho <- gsub("\\|", ":", ortho)
        grepID <- paste(group, "#", ortho, sep = "")
        outDf <- outDf[outDf$seedID == grepID, ]

        # filter domain df by features
        if (nrow(outDf) == 0) stop("Domain info is NULL!")

        outDf[c("feature_type","feature_id")] <-
            stringr::str_split_fixed(outDf$feature, '_', 2)
        outDf <- outDf[!(outDf$feature_type %in% input$feature),]

        # filter filters without e-value and/or bitscore
        if ("evalue" %in% colnames(outDf)) {
            if ("noEvalue" %in% input$feature)
                outDf <- outDf[!is.na(outDf$evalue),]
            if ("noBitscore" %in% input$feature)
                outDf <- outDf[!is.na(outDf$bitscore),]
        }

        # modify feature IDs
        outDf$feature_id_mod <- outDf$feature_id
        outDf$feature_id_mod <- gsub("SINGLE", "LCR", outDf$feature_id_mod)
        outDf$feature_id_mod[outDf$feature_type == "coils"] <- "Coils"
        outDf$feature_id_mod[outDf$feature_type == "seg"] <- "LCR"
        outDf$feature_id_mod[outDf$feature_type == "tmhmm"] <- "TM"
        # exclude features IDs
        if (!is.null(input$excludeNames)) {
            outDf$feature_id_mod[outDf$feature_type %in% input$excludeNames]<-NA
        }

        # enable/disable option for showing evalue/bitscore
        if ("evalue" %in% colnames(outDf)) {
            shinyjs::enable("showScore")
        } else {
            shinyjs::disable("showScore")
        }

        # enable/disable option for showing weight
        if ("weight" %in% colnames(outDf)) {
            shinyjs::enable("showWeight")
        } else {
            shinyjs::disable("showWeight")
        }

        # Filter data by e-value, bit-score and feature path
        if ("evalue" %in% colnames(outDf)) {
            # filter by e-value and/or bit-score
            if ("E-value" %in% input$showScore) {
                minEvalue <- format(input$minEvalue, scientific = FALSE)
                naOutDf <- outDf[is.na(outDf$evalue),]
                outDf <- outDf[!is.na(outDf$evalue) & outDf$evalue <= minEvalue,]
                outDf <- rbind(outDf,naOutDf)
            }
            if ("Bit-score" %in% input$showScore) {
                naOutDf <- outDf[is.na(outDf$bitscore),]
                outDf <- outDf[
                    !is.na(outDf$bitscore) & outDf$bitscore >= input$minBitscore,
                ]
                outDf <- rbind(outDf,naOutDf)
            }
            # get only best instances
            if ("evalue" %in% input$linearizationBy) {
                linearizedDfs <- lapply(
                    levels(as.factor(outDf$orthoID)),
                    function(orthoID) {
                        return(linearizeArchitecture(outDf, orthoID, "evalue"))
                    }
                )
                outDf <- do.call(rbind, linearizedDfs)
            }
            if ("bitscore" %in% input$linearizationBy) {
                linearizedDfs <- lapply(
                    levels(as.factor(outDf$orthoID)),
                    function(orthoID) {
                        return(linearizeArchitecture(outDf, orthoID,"bitscore"))
                    }
                )
                outDf <- do.call(rbind, linearizedDfs)
            }
            if ("path" %in% input$linearizationBy) {
                outDf <- outDf %>% dplyr::group_by(feature) %>%
                    dplyr::filter(path == "Y")
            }
            # Format e-values
            outDf$evalue[!is.na(outDf$evalue)] <-
                format(
                    outDf$evalue[!is.na(outDf$evalue)], scientific = TRUE,
                    digits = 2
                )
        }
        return(outDf[!is.na(outDf$seedID),])
    })

    # * render plot ============================================================
    output$archiPlot <- renderPlot({
        if (is.null(nrow(filterDomainDf()))) stop("Domain info is NULL!")
        # remove user specified features (from input$featureList)
        df <- filterDomainDf()
        df <- df[!(df$feature_id %in% input$featureList),]
        g <- createArchiPlot(
            pointInfo(), df,
            input$labelArchiSize, input$titleArchiSize, input$legendArchiSize,
            input$showScore, input$showWeight, input$namePostion, 
            input$firstDist,input$nameType, input$nameSize, input$segmentSize, 
            input$nameColor, input$labelPos, input$colorType, 
            input$ignoreInstanceNo, currentNCBIinfo(), input$featureClassSort, 
            input$featureClassOrder, input$colorPallete, input$resolveOverlap, 
            font()
        )
        if (is.character(g) && any(g == "No domain info available!")) {
            msgPlot()
        } else {
            suppressWarnings(grid::grid.draw(g))
        }
    })

    output$archiPlot.ui <- renderUI({
        ns <- session$ns
        if (is.null(nrow(filterDomainDf()))) {
            msg <- paste0(
                "<p><em>Wrong domain file has been uploaded!
        Please check the correct format in
        <a href=\"https://github.com/BIONF/PhyloProfile/wiki/",
                "Input-Data#ortholog-annotations-eg-domains\"
        target=\"_blank\" rel=\"noopener\">our PhyloProfile wiki</a>.</em></p>"
            )
            HTML(msg)
        } else {
            shinycssloaders::withSpinner(
                plotOutput(
                    ns("archiPlot"),
                    height = input$archiHeight,
                    width = input$archiWidth,
                    click = ns("archiClick")
                )
            )
        }
    })

    output$archiDownload <- downloadHandler(
        filename = function() {
            c("domains.svg")
        },
        content = function(file) {
            # remove user specified features (from input$featureList)
            df <- filterDomainDf()
            df <- df[!(df$feature_id %in% input$featureList),]
            g <- createArchiPlot(
                pointInfo(), df,
                input$labelArchiSize,input$titleArchiSize,input$legendArchiSize, 
                input$showScore, input$showWeight, input$namePostion, 
                input$firstDist, input$nameType, input$nameSize, 
                input$segmentSize, input$nameColor, input$labelPos, 
                input$colorType, input$ignoreInstanceNo, currentNCBIinfo(),
                input$featureClassSort, input$featureClassOrder,
                input$colorPallete, input$resolveOverlap, font()
            )
            # save plot to file
            suppressWarnings(ggsave(
                file, plot = g,
                width = input$archiWidth * 0.035,
                height = input$archiHeight * 0.035,
                units = "cm", dpi = 300, device = "svg", limitsize = FALSE
            ))
        }
    )

    # output$selectedDomain <- renderText({
    #     if (is.null(input$archiClick$y)) return("No domain selected!")
    #     y <- input$archiClick$y
    #     # paste(y, round(y), convertY(unit(y, "npc"), "px"))
    #
    # })

    output$featureList.ui <- renderUI({
        ns <- session$ns
        allFeats <- getAllFeatures(pointInfo(), filterDomainDf())
        selectizeInput(
            ns("featureList"), "Exclude individual feature(s)", multiple = TRUE,
            choices = allFeats,
            options=list(placeholder = 'None')
        )
    })

    output$linkTable <- renderTable({
        if (is.null(nrow(filterDomainDf()))) return("No domain info available!")
        features <- getDomainLink(pointInfo(), filterDomainDf())
        features <- features[!duplicated(features),]
    }, sanitize.text.function = function(x) x)

    output$domainTable <- DT::renderDataTable({
        req(filterDomainDf())
        createDomainInfoTable(filterDomainDf())
    }, rownames= FALSE)
}

#' plot error message
#' @return error message in a ggplot object
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
msgPlot <- function() {
    msg <- paste(
        "No information about domain architecture!",
        "Please check:","if you uploaded the correct domain file/folder; or ",
        "if the selected genes (seed & ortholog) do exist in the uploaded file",
        "(please search for the corresponding seedID and hitID)",
        sep = "\n"
    )
    x <- c(1,2,3,4,5)
    y <- c(1,2,3,4,5)
    g <- ggplot(data.frame(x, y), aes(x,y)) +
        geom_point(color = "white") +
        annotate(
            "text", label = msg, x = 3.5, y = 0.5, size = 5, colour = "red"
        ) +
        theme(axis.line = element_blank(), axis.text = element_blank(),
              axis.ticks = element_blank(), axis.title = element_blank(),
              panel.background = element_blank(),
              panel.border = element_blank(),
              panel.grid = element_blank(),
              plot.background = element_blank()) +
        ylim(0,1)
    return(g)
}

#' Get list of all features
#' @return A list of all features
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
getAllFeatures <- function(info, domainDf) {
    group <- as.character(info[1])
    ortho <- as.character(info[2])
    # get sub dataframe based on selected groupID and orthoID
    group <- gsub("\\|", ":", group)
    ortho <- gsub("\\|", ":", ortho)
    grepID <- paste(group, "#", ortho, sep = "")
    subdomainDf <- domainDf[grep(grepID, domainDf$seedID), ]
    orthoID <- feature <- NULL
    if (nrow(subdomainDf) < 1) return(paste0("No domain info available!"))
    else {
        # ortho & seed domains df
        orthoDf <- subdomainDf[subdomainDf$orthoID == ortho,]
        seedDf <- subdomainDf[subdomainDf$orthoID != ortho,]
        feature <- c(
            levels(as.factor(orthoDf$feature_id)),
            levels(as.factor(seedDf$feature_id))
        )
    }
    return(levels(as.factor(feature)))
}

#' get pfam and smart domain info (domain name, acc, profile HMM,...)
#' @return dataframe for each type of database
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
getDomainInfo <- function(info, domainDf, type) {
    group <- as.character(info[1])
    ortho <- as.character(info[2])
    # get sub dataframe based on selected groupID and orthoID
    group <- gsub("\\|", ":", group)
    ortho <- gsub("\\|", ":", ortho)
    grepID <- paste(group, "#", ortho, sep = "")
    domainDf <- domainDf[grep(grepID, domainDf$seedID),]
    domainDf <- domainDf[domainDf$feature_type == type,]
    domainInfoDf <- data.frame()
    if (nrow(domainDf) > 0) {
        # get all features of for this pair proteins
        feature <- getAllFeatures(info, domainDf)
        # filter domain info
        if ("acc" %in% colnames(domainDf)) {
            domainInfoDf <- domainDf[
                , c("orthoID", "feature_id", "acc", "evalue", "bitscore")
            ]
            domainInfoDf$evalue <- format(domainInfoDf$evalue, scientific=TRUE)
        } else {
            domainInfoDf <- domainDf[, c("orthoID", "feature_id")]
            domainInfoDf <- domainInfoDf[!(duplicated(domainInfoDf)),]
            domainInfoDf$acc <- rep(NA, nrow(domainInfoDf))
            domainInfoDf$evalue <- rep(NA, nrow(domainInfoDf))
            domainInfoDf$bitscore <- rep(NA, nrow(domainInfoDf))
        }
    }
    return(domainInfoDf)
}

#' get pfam and smart domain links
#' @return dataframe with domain IDs and their database links
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
getDomainLink <- function(info, domainDf) {
    featurePfam <- getDomainInfo(info, domainDf, "pfam")
    pfamDf <- createLinkTable(featurePfam, "pfam")
    if (nrow(pfamDf) > 0) pfamDf$db <- "PFAM"
    featureSmart <- getDomainInfo(info, domainDf, "smart")
    smartDf <- createLinkTable(featureSmart, "smart")
    if (length(smartDf) > 0) smartDf$db <- "SMART"
    featDf <- rbind(pfamDf, smartDf)
    if (nrow(featDf) < 1) stop("No PFAM/SMART domain found!")
    featDf <- subset(featDf, select = c(feature_id, db, link))
    colnames(featDf) <- c("Feature ID", "Source", "URL")
    return(featDf)
}

#' plot error message
#' @param featureDf dataframe contains 2 columns feature_id and acc
#' @param featureType pfam or smart
#' @return dataframe contains 2 columns feature_id and link
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
createLinkTable <- function(featureDf, featureType) {
    featureDf <- featureDf[!(duplicated(featureDf)),]
    if (nrow(featureDf) > 0) {
        if (featureType == "pfam") {
            featureDf$link[is.na(featureDf$acc)] <- paste0(
                "<a href='http://pfam-legacy.xfam.org/family/",
                featureDf$feature_id,
                "' target='_blank'>", "PFAM", "</a>"
            )
            featureDf$link[!(is.na(featureDf$acc))] <- paste0(
                "<a href='https://www.ebi.ac.uk/interpro/entry/pfam/",
                featureDf$acc,
                "' target='_blank'>", "INTERPRO", "</a>"
            )
        } else {
            featureDf$link <- paste0(
                "<a href='http://smart.embl-heidelberg.de/smart/",
                "do_annotation.pl?BLAST=DUMMY&DOMAIN=",
                # featureDf$feature_id, "' target='_blank'>",
                gsub("_smart", "",featureDf$feature_id), "' target='_blank'>",
                "SMART", "</a>"
            )
        }
    }
    return(featureDf)
}

#' plot error message
#' @param domainDf dataframe contains feature information
#' @return dataframe as subset of domainDf
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

createDomainInfoTable <- function(domainDf) {
    if (is.null(domainDf)) stop("No information!")

    selectedCols <- c("orthoID", "length", "feature", "start", "end")
    if (ncol(domainDf) <= 11) {
        if (!("length" %in% colnames(domainDf)))
            selectedCols <- selectedCols[!selectedCols == "length"]
        outDf <- subset(domainDf, select = selectedCols)
        if ("length" %in% colnames(domainDf)) {
            colnames(outDf) <- c(
                "orthoID", "Length", "Feature", "Start", "End"
            )
        } else {
            colnames(outDf) <- c(
                "orthoID", "Feature", "Start", "End"
            )
        }
    } else if (ncol(domainDf) == 17 & "evalue" %in% colnames(domainDf)) {
        selectedCols <- c(
            selectedCols, "evalue", "bitscore", "pStart", "pEnd", "pLen"
        )
        outDf <- subset(domainDf, select = selectedCols)
        colnames(outDf) <- c(
            "orthoID", "Length", "Feature", "Start", "End", "E-value",
            "Bit-score", "pHMM start", "pHMM end", "pHMM length"
        )
    } else {
        return(paste("Wrong number of columns:", ncol(domainDf)))
    }
    outDf$Feature <- sub("_", " ", outDf$Feature)
    outDf$`Gene ID` <- unlist(lapply(
        strsplit(outDf$orthoID, split = ":"), function (x) return(x[3])
    ))
    outDf <- subset(outDf, select = -c(orthoID))
    outDf <- outDf %>% dplyr::select(`Gene ID`, everything())
    return(outDf)
}
