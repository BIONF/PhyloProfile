#' Group Comparison
#'
#' @export
#' @param selectedInGroup selected In-group (input$selectedInGroupGC)
#' @param selectedGenesList list of genes to calculate the plots for
#'  (input$listSelectedGenesGC)
#' @param mainRank in input and settings selected taxonomy rank
#'  (input$rankSelect)
#' @param selectedVariable variable(s) to claculate the plots for
#'  (input$varNameGC)
#' @param useCommonAncestor boolean if the next common anchestor should be
#'  used (input$useCommonAncestor)
#' @param referenceTaxon selected taxon name (input$inSelect)
#' @param ncbiIDList list of ncbi ids (from reactive fn "inputTaxonID")
#' @param filteredData full processed main data
#' (from reactive fn "getDataFiltered")
#' @param rightFormatFeatures boolean if the features have the right format
#' (input$rightFormatFeatures)
#' @param domainInformation dataframe holding the domain input
#' (from reactive fn "getDomainInformation)
#' @param plot information if the plots should be generated (input$plotGC)
#' @param parameter list of parameters needed to generate the plots
#' (from reactive fn "getParameterInputGC")
#' @return list of candidate genes
#' @author Carla Mölbert {carla.moelbert@gmx.de}

source("R/functions.R")

groupComparisonUI <- function(id){
    ns <- NS(id)
    fluidPage(
        sidebarPanel(
            shinycssloaders::withSpinner(uiOutput(ns("candidateGenes"))),
            bsPopover(
                "candidateGenes",
                "",
                "Select gene to show the plots",
                "right"
            ),
            shinycssloaders::withSpinner(uiOutput(ns("featuresOfInterestUI"))),
            bsPopover(
                "featuresOfInterestUI",
                "",
                "This function is only use full if the features are
                saved in the right format: featuretype_featurename"
            ),
            sliderInput(
                ns("domainsThreshold"),
                "Persentage of proteins in which the domain needs to have an
                instance",
                min = 0 , max = 100, value = 0,
                step = 1, round = FALSE
            ),

            downloadButton(ns("downloadPlots"), "Download plots"),
            width = 3
        ),
        mainPanel(
            tags$style(
                HTML("#plotsUI { height:650px; overflow-y:scroll}")
            ),
            shinycssloaders::withSpinner(uiOutput(ns("plotsUI"))),
            width = 9
        )
    )
}

groupComparison <- function(input, output, session,
                            selectedInGroup,
                            selectedGenesList,
                            mainRank,
                            selectedVariable,
                            useCommonAncestor,
                            referenceTaxon,
                            ncbiIDList,
                            filteredData,
                            rightFormatFeatures,
                            domainInformation,
                            plot,
                            parameter,
                            selectedPoint){
    # Dataframe for the significant Genes ======================================
    #' contains geneID, inGroup, outGroup, pvalues, features, databases,
    #' rank, var
    candidateGenes <- reactiveValues(plots = NULL)

    # List with all candidate genes ============================================
    output$candidateGenes <- renderUI({
        ns <- session$ns
        plot()

        isolate({
            candidateGenes$plots <- {
                getSignificantGenes(
                    selectedInGroup(),
                    selectedGenesList(),
                    mainRank(),
                    selectedVariable(),
                    useCommonAncestor(),
                    referenceTaxon(),
                    parameter(),
                    ncbiIDList(),
                    filteredData(),
                    rightFormatFeatures(),
                    domainInformation()
                )
            }

            if (is.data.frame(candidateGenes$plots)) {
                significantGenes <- candidateGenes$plots
                x <- as.vector(significantGenes$geneID)
                choices <- c("all", x)

                selectInput(ns("selectedGene"), "Candidate gene(s):",
                            choices,
                            selected = choices[2],
                            multiple = FALSE)

            } else{
                selectInput(ns("selectedGene"), "Candidate gene(s):",
                            NULL,
                            selected = NULL,
                            multiple = FALSE)
            }
        })
    })

    # Output of the plots for the selected gene(s) =============================
    output$plotsUI <- renderUI({
        if (is.character(candidateGenes$plots)) return(candidateGenes$plots)
        getPlots()
    })

    # List with possible features for the selected gene ========================
    output$featuresOfInterestUI <- renderUI({
        ns <- session$ns
        input$selectedGene
        isolate({
            gene <- input$selectedGene
            if (!rightFormatFeatures()) {
                selectInput(ns("interestingFeatures"),
                            "Feature type(s) of interest:",
                            NULL,
                            selected = NULL,
                            multiple = TRUE,
                            selectize = FALSE)
            } else if (is.null(gene)) {
                selectInput(ns("interestingFeatures"),
                            "Feature type(s) of interest:",
                            NULL,
                            selected = NULL,
                            multiple = TRUE,
                            selectize = FALSE)
            } else if (gene == "") {
                selectInput(ns("interestingFeatures"),
                            "Feature type(s) of interest:",
                            NULL,
                            selected = NULL,
                            multiple = TRUE,
                            selectize = FALSE)
            } else{
                significantGenes <- candidateGenes$plots
                choices <- c("all")
                if (gene == "all") {
                    for (currentGene in significantGenes$geneID) {
                        subsetCurrentGene <-
                            subset(significantGenes,
                                   significantGenes$geneID == currentGene)
                        choices <- append(
                            choices, unlist(subsetCurrentGene$databases)
                        )
                    }
                    # show each database only once
                    choices <- choices[!duplicated(choices)]
                } else {

                    subsetGene <- subset(significantGenes,
                                          significantGenes$geneID == gene)

                    choices <- append(choices, unlist(subsetGene$databases))

                }
                selectInput(ns("interestingFeatures"),
                            "Feature type(s) of interest:",
                            choices,
                            selected = choices[1],
                            multiple = TRUE,
                            selectize = FALSE)
            }
        })
    })

    # download file with the shown plots =======================================
    output$downloadPlots <- downloadHandler(
        filename = "plotSignificantGenes.zip",
        content = function(file){
            genes <- input$selectedGene
            significantGenes <- candidateGenes$plots
            if ("all" %in% genes) {
                genes <- significantGenes$geneID
            }

            fs <- c()
            #tmpdir <- tempdir()
            setwd(tempdir())

            for (gene in genes) {
                path <- paste(gene, ".pdf", sep = "")
                fs <- c(fs, path)
                pdf(path)
                getPlotsToDownload(
                    gene, parameter(), input$interestingFeatures, 
                    significantGenes,input$domainsThreshold, selectedPoint()
                )
                dev.off()
            }
            zip(zipfile = file, files = fs)
        },
        contentType = "application/zip"
    )

    #' observer for the download functions
    observe({
        if (is.null(selectedInGroup())
            | length(selectedGenesList()) == 0) {
            shinyjs::disable("downloadPlots")
        } else if (plot() == FALSE) {
            shinyjs::disable("downloadPlots")
        } else if (input$selectedGene == "") {
            shinyjs::disable("downloadPlots")
        } else {
            shinyjs::enable("downloadPlots")
        }
    })

    # Deciding which plots will be shown =======================================
    getPlots <- reactive({
        input$interestingFeatures
        input$domainsThreshold
        gene <- as.character(input$selectedGene)
        plot()
        if (is.null(candidateGenes$plots)) return()
        significantGenes <- candidateGenes$plots
        if (gene == "all") {
            plotOutputList <- getPlotOutputList(significantGenes,
                                                parameter(),
                                                input$interestingFeatures,
                                                input$domainsThreshold,
                                                selectedPoint())
        }else{
            geneInfo <- {
                significantGenes[significantGenes$geneID == gene, ]
            }
            if (nrow(geneInfo) == 0) return()
            plotOutputList <- getPlotOutputList(geneInfo,
                                                     parameter(),
                                                     input$interestingFeatures,
                                                     input$domainsThreshold,
                                                     selectedPoint())
        }
        #' List with all plots that will be shown
        return(plotOutputList)
    })

    # List of genes for the customized profile =================================
    geneList <- reactive({
        if (!is.null(candidateGenes$plots)) {
            significantGenes <- candidateGenes$plots
            return(significantGenes$geneID)
        }
    })
    return(geneList)
}

# FUNCTIONS ===================================================================

#' Generate the list with all plots -------------------------------------------
#' @export
#' @param genes list of genes
#' @param parameters contains "showPValue","highlightSignificant",
#' "significance", "var1ID", "var2ID", "xSizeGC", "ySizeGC",
#' "interestingFeatures", "angleGC", "legendGC", "legendSizeGC"
#' @param interestingFeatures list of databases to take the features from
#' @return list with all plots
#' @author Carla Mölbert (carla.moelbert@gmx.de)
getPlotOutputList <- function(genes, parameters, interestingFeatures,
                                 domainsThreshold, selectedPoint) {
    if (is.null(genes)) return()
    # Insert plot output objects the list
    plotOutputList <- lapply(1:nrow(genes), function(i) {
        plotname <- paste(genes[i, 1])
        plotOutputObject <- renderPlot(getMultiplot(genes[i, ], parameters,
                                                    interestingFeatures,
                                                    domainsThreshold,
                                                    selectedPoint),
                                         height = 650, width = 700)
    })
    do.call(tagList, plotOutputList) # needed to display properly.
    return(plotOutputList)
}

#' Put the plots for one spicific gene in one multiplot -----------------------
#' @export
#' @param geneInfo contains "geneID",  "inGroup",  "outGroup", "pvalues",
#' "features",  "databases", "var", "rank"
#' @param parameters contains "showPValue","highlightSignificant",
#' "significance", "var1ID", "var2ID", "xSizeGC", "ySizeGC",
#' "interestingFeatures", "angleGC", "legendGC", "legendSizeGC"
#' @param interestingFeatures list of databases to take the features from
#' @return grid arrange with the plots that should be shown for this gene
#' @author Carla Mölbert (carla.moelbert@gmx.de)
getMultiplot <- function(
    geneInfo, parameters, interestingFeatures, domainsThreshold, selectedPoint
){
    #' Sorting the information to the selected gene ----------------------------
    gene <- as.character(geneInfo$geneID)
    inGroup <- as.data.frame(geneInfo$inGroup)
    outGroup <- as.data.frame(geneInfo$outGroup)
    features <- as.data.frame(geneInfo$features)

    var <- geneInfo$var

    #' Get the barplot ---------------------------------------------------------
    barplot <-  getBarplotGC(gene,
                               inGroup,
                               outGroup,
                               features,
                               parameters,
                               interestingFeatures,
                               domainsThreshold)

    if (is.null(barplot)) {
        barplot <- textGrob("The selected domains are not found in the gene")
    }

    #' Get the boxplots  for two variables  ------------------------------------
    if (var == "Both") {
        pValue1 <- round(geneInfo$pvalues1, 2)
        pValue2 <- round(geneInfo$pvalues2, 2)

        #' Check if the pValues should be printed
        if (parameters$showPValue == TRUE) {
            infoPValue1 <- paste("P-value:", pValue1, sep = " ")
            infoPValue2 <- paste("P-value:", pValue2, sep = " ")
        } else {
            infoPValue1 <- " "
            infoPValue2 <- " "
        }

        #' Get information about the plot colour
        if (parameters$highlightSignificant == TRUE) {
            if (is.na(pValue1)) colour1 <- "grey"
            else if (pValue1 < parameters$significance) {
                colour1 <- "indianred2"
            } else colour1 <- "grey"

            if (is.na(pValue2)) colour2 <- "grey"
            else if (pValue2 < parameters$significance) {
                colour2 <- "indianred2"
            } else colour2 <- "grey"
        } else {
            colour1 <- "grey"
            colour2 <- "grey"
        }

        #' Generate the boxplots
        boxplot1 <- getBoxplotGC(inGroup,
                                   outGroup,
                                   parameters$var1ID,
                                   gene,
                                   colour1,
                                   infoPValue1, parameters,
                                   selectedPoint)

        boxplot2 <- getBoxplotGC(inGroup,
                                   outGroup,
                                   parameters$var2ID,
                                   gene,
                                   colour2,
                                   infoPValue2, parameters,
                                   selectedPoint)

        plots <- grid.arrange(textGrob(gene),
                              arrangeGrob(boxplot1, boxplot2, ncol = 2),
                              barplot,
                              heights = c(0.02, 0.45, 0.458), ncol = 1)
    }
    #' get the boxplot if one varibale is selected  ----------------------------
    else {
        p <- round(geneInfo$pvalue, 2)

        #' Check if the pValues should be printed
        if (parameters$showPValue == TRUE) {
            info <- paste("P-value:", p)
        } else {
            info <- " "
        }

        #' Generate the plot
        boxplot <- getBoxplotGC(inGroup,
                                  outGroup,
                                  var,
                                  gene,
                                  "grey",
                                  info,
                                  parameters,
                                  selectedPoint)

        plots <- grid.arrange(textGrob(gene),
                              boxplot,
                              barplot,
                              heights = c(0.02, 0.45, 0.458), ncol = 1)
    }

    #' return the plots --------------------------------------------------------
    return(plots)
}

#' Create a Boxplot -----------------------------------------------------------
#' @export
#' @param inGroupDf contains "supertaxon", "geneID", "ncbiID", "orthoID",
#'  "var1", "var2", "paralog", "abbrName", "taxonID", "fullname",
#'  "supertaxonID", "rank", "category", "presSpec", "mVar1", "mVar2"
#' @param outGroupDf contains "supertaxon", "geneID", "ncbiID", "orthoID",
#'  "var1", "var2", "paralog", "abbrName", "taxonID", "fullname",
#' "supertaxonID", "rank", "category", "presSpec", "mVar1", "mVar2"
#' @param var variable to consider in the boxplot
#' @param  gene gene for which the plot is generated
#' @param  colour colour of the boxes
#' @param  info info about the p-values
#' @param parameters contains "showPValue","highlightSignificant",
#' "significance", "var1ID", "var2ID", "xSizeGC", "ySizeGC",
#' "interestingFeatures", "angleGC", "legendGC", "legendSizeGC"
#' @return boxplot
#' @author Carla Mölbert (carla.moelbert@gmx.de)
getBoxplotGC <- function(inGroupDf,
                           outGroupDf,
                           var,
                           gene,
                           colour,
                           info,
                           parameters,
                           selectedPoint){

    #' pre-processing the data for the boxplot ---------------------------------
    if (var == parameters$var1ID) {
        inGroup <- inGroupDf$var1
        outGroup <- outGroupDf$var1
    } else if (var == parameters$var2ID) {
        inGroup <- inGroupDf$var2
        outGroup <- outGroupDf$var2
    }

    inGroup <- inGroup[!is.na(inGroup)]
    outGroup <- outGroup[!is.na(outGroup)]

    lengthInGroup <- length(inGroup)
    lengthOutGroup <- length(outGroup)

    inGroup <- as.data.frame(inGroup)
    names(inGroup)[1] <- paste("values")
    inGroup$group <- "inGroup"

    outGroup <- as.data.frame(outGroup)
    names(outGroup)[1] <- paste("values")
    outGroup$group <- "Out-Group"

    dataBoxplot <- rbind(inGroup, outGroup)
    dataBoxplot <- dataBoxplot[complete.cases(dataBoxplot), ]

    names <- c(paste("In-Group \n n=", lengthInGroup, sep = ""),
               paste("Out-Group \n n=", lengthOutGroup, sep = ""))

    #' Generate the boxplot ----------------------------------------------------
    boxplotGC <- ggplot(dataBoxplot, aes(group, values)) +
        geom_violin(position = position_dodge(), scale = "width",
                    fill = colour) +
        labs(x = "", y = var, caption = paste(info), colour = "") +
        scale_x_discrete(labels = names) +
        theme_minimal() +
        stat_summary(aes(colour = selectedPoint),fun.y = selectedPoint,
                     geom = "point", size = 3, show.legend = TRUE)

    boxplotGC <- boxplotGC +
        theme(axis.text.x = element_text(size = parameters$xSizeGC,
                                         hjust = 1),
              axis.text.y = element_text(size = parameters$ySizeGC),
              axis.title.y = element_text(size = parameters$ySizeGC),
              plot.caption = element_text(size = parameters$pValuesSize),
              legend.position = "bottom",
              legend.text = element_text(size = parameters$legendSizeGC ),
              legend.title = element_text(size = parameters$legendSizeGC)) +
        scale_color_manual("", values = c("green"))

    #' return the boxplot ------------------------------------------------------
    return(boxplotGC)
}


#' Create a Barplot -----------------------------------------------------------
#' @export
#' @param  selectedGene gene for which the plot is generated
#' @param inGroupDf contains "supertaxon", "geneID", "ncbiID", "orthoID",
#'  "var1", "var2", "paralog", "abbrName", "taxonID", "fullname",
#'  "supertaxonID", "rank", "category", "presSpec", "mVar1", "mVar2"
#' @param outGroupDf contains "supertaxon", "geneID", "ncbiID", "orthoID",
#'  "var1", "var2", "paralog", "abbrName", "taxonID", "fullname",
#' "supertaxonID", "rank", "category", "presSpec", "mVar1", "mVar2"
#' @param features contains "seedID",  "orthoID", "feature", "start",   "end"
#' @param parameters contains "showPValue","highlightSignificant",
#' "significance", "var1ID", "var2ID", "xSizeGC", "ySizeGC",
#' "interestingFeatures", "angleGC", "legendGC", "legendSizeGC"
#' @param interestingFeatures list of databases for which the features should
#' be included
#' @return barplot
#' @author Carla Mölbert (carla.moelbert@gmx.de)
getBarplotGC <- function(selectedGene,
                           inGroup, outGroup,
                           features, parameters,
                           interestingFeatures,
                           domainsThreshold){

    #' information about the features ------------------------------------------
    features$feature <- as.character(features$feature)

    #' part in inGroup and out-group  -----------------------------------------
    inGroup$orthoID <- gsub("\\|", ":", inGroup$orthoID)
    inGroupDomainDf  <-  {
        subset(features, features$orthoID %in% inGroup$orthoID)
    }
    outGroup$orthoID <- gsub("\\|", ":", outGroup$orthoID)
    outGroupDomainDf <- {
        subset(features, features$orthoID %in% outGroup$orthoID)
    }

    #' get the dataframes for in- and out-group --------------------------------
    if (nrow(inGroupDomainDf) == 0) {
        dataIn <- NULL
    } else {
        feature <- unique(inGroupDomainDf$feature)
        dataIn <- as.data.frame(feature)
        dataIn$amount <- 0
        dataIn$proteins <- 0
    }

    if (nrow(outGroupDomainDf) == 0) {
        dataOut <- NULL
    } else {
        feature <- unique(outGroupDomainDf$feature)
        dataOut <- as.data.frame(feature)
        dataOut$amount <- 0
        dataOut$proteins <- 0
    }

    seeds <- features$seedID
    seeds <- unique(seeds)

    #' Get the values for the boxplot ------------------ -----------------------
    inNotEmpty <- 0
    outNotEmpty <- 0

    #' Count for each feature how often it is present in each seed -------------
    for (seed in seeds) {
        #' count the features in the in-group
        if (!is.null(dataIn)) {
            inG <- subset(inGroupDomainDf,
                           inGroupDomainDf$seedID == seed)
            inG <- inG[str_detect(inG$seedID, inG$orthoID),]

            if (!empty(inG)) {
                inNotEmpty <- inNotEmpty + 1
                inGroupFeatures <-  plyr::count(inG, "feature")
                for (i in 1:nrow(inGroupFeatures)) {
                    for (j in 1:nrow(dataIn)) {
                        if (dataIn[j, 1] == inGroupFeatures[i, 1]) {
                            dataIn[j, 2] <-
                                dataIn[j, 2] + inGroupFeatures[i, 2]
                            dataIn[j, 3] <- dataIn[j,3] + 1
                        }
                    }
                }
            }
        }

        #' count the featueres in the out-group
        if (!is.null(dataOut)) {
            OutG <- subset(outGroupDomainDf,
                            outGroupDomainDf$seedID == seed)
            OutG <- OutG[str_detect(OutG$seedID, OutG$orthoID),]

            if (!empty(OutG)) {
                outNotEmpty <- outNotEmpty + 1
                outGroupFeatures <-  plyr::count(OutG, "feature")
                for (i in 1:nrow(outGroupFeatures)) {
                    for (j in 1:nrow(dataOut)) {
                        if (dataOut[j, 1] == outGroupFeatures[i, 1]) {
                            dataOut[j, 2] <-
                                dataOut[j, 2] + outGroupFeatures[i, 2]
                            dataOut[j, 3] <- dataOut[j,3] + 1
                        }
                    }
                }
            }
        }
    }

    if (domainsThreshold > 0) {
        thresholdIn  <- (domainsThreshold / 100) * inNotEmpty
        thresholdOut <- (domainsThreshold / 100) * outNotEmpty

        dataIn <- dataIn[sort(dataIn$feature),]
        dataOut <- dataOut[sort(dataOut$feature),]

        featuresIn <- dataIn$feature[!(dataIn$feature %in% dataOut$feature)]
        featuresOut <-
            dataOut$feature[!(dataOut$feature %in% dataIn$feature)]

        for (i in featuresOut) {
            dataIn <- rbind(dataIn, c(i, 0, 0))
        }

        for (i in featuresIn) {
            dataOut <- rbind(dataOut, c(i, 0, 0))
        }

        dataIn$proteins <- as.numeric(dataIn$proteins)
        dataOut$proteins <- as.numeric(dataOut$proteins)

        dataInKeep <- dataIn[dataIn$proteins >= thresholdIn,]
        dataOutKeep <- dataOut[dataOut$proteins >= thresholdOut, ]

        keep <- rbind(dataInKeep, dataOutKeep)
        keep <- keep[!duplicated(keep$feature),]

        dataIn <- dataIn[dataIn$feature %in% keep$feature,]
        dataOut <- dataOut[dataOut$feature %in% keep$feature,]
    }

    #' Calculate the average of appearances for each feature -------------------
    if (!is.null(dataIn)) {
        dataIn$amount <- as.numeric(dataIn$amount)
        dataIn$amount <- dataIn$amount / inNotEmpty
        dataIn$type <- "In-Group"
    }

    if (!is.null(dataOut)) {
        dataOut$amount <- as.numeric(dataOut$amount)
        dataOut$amount <- dataOut$amount / outNotEmpty
        dataOut$type <- "Out-Group"
    }

    #' Get the data for teh barplot --------------------------------------------
    if (is.null(dataIn) & !is.null(dataOut)) {
        dataBarplot <- dataOut
    }  else if (is.null(dataOut) & !is.null(dataIn)) {
        dataBarplot <- dataIn
    } else if (!is.null(dataIn) & !is.null(dataOut)) {
        dataBarplot <- rbind(dataIn, dataOut)
    } else {
        dataBarplot <- NULL
    }

    #' only show features that interest the user
    if (!("all" %in% interestingFeatures)) {
        featuresList <- NULL
        for (feature in interestingFeatures) {
            dataBarplot$feature <- as.character(dataBarplot$feature)
            subsetFeatures <- subset(dataBarplot$feature,
                                      startsWith(dataBarplot$feature, feature))
            featuresList <- append(featuresList, subsetFeatures)
        }

        if (is.null(featuresList)) return()
        #' only keep rows in which the feature begins with a element out of the
        #' interesing Features
        featuresList <- featuresList[!duplicated(featuresList)]
        dataBarplot <- subset(dataBarplot, dataBarplot$feature
                               %in% featuresList)
    }

    dataBarplot <- dataBarplot[sort(dataBarplot$feature),]

    #' generate the barplot ----------------------------------------------------
    if (!is.null(dataBarplot)) {
        barplotGC <- ggplot(dataBarplot,
                             aes(x = feature, y = amount, fill = type ),
                             main  = " ") +
            geom_bar(
                stat = "identity", position = position_dodge(), width = 0.5
            ) +
            scale_fill_grey() +
            labs(x = " ", y = "Average instances per protein", fill = "Group") +
            theme_minimal()

        barplotGC <- barplotGC +
            theme(axis.text.x = element_text(size = parameters$xSizeGC,
                                             angle = parameters$angleGC,
                                             hjust = 1),
                  axis.text.y = element_text(size = parameters$ySizeGC),
                  axis.title.y = element_text(size = parameters$ySizeGC),
                  legend.position = parameters$legendGC,
                  legend.text = element_text(size = parameters$legendSizeGC ),
                  legend.title = element_text(size = parameters$legendSizeGC))

        #' return the barplot --------------------------------------------------
        return(barplotGC)
    } else (return(NULL))
}

#' get the plots to download --------------------------------------------------
#' @export
#' @param  gene gene for which the plot is generated
#' @param interestingFeatures list of databases for which the features should
#' be included
#' @param parameters contains "showPValue","highlightSignificant",
#' "significance", "var1ID", "var2ID", "xSizeGC", "ySizeGC",
#' "interestingFeatures", "angleGC", "legendGC", "legendSizeGC"
#' @return arrange grop containing the plots for this gene
#' @author Carla Mölbert (carla.moelbert@gmx.de)
getPlotsToDownload <- function(
    gene, parameters, interestingFeatures, significantGenes, domainThreshold,
    selectedPoint
){
    infoGene <- subset(significantGenes,
                        significantGenes$geneID == gene)
    return(getMultiplot(infoGene, parameters, interestingFeatures,
                         domainsThreshold,
                         selectedPoint))
}

#' Get the pValues to print under the plot -----------------------------------
#' @export
#' @param pvalues list contianing the p-values for a specific gene
#' @return string containing the information about the p-values
#' @author Carla Mölbert (carla.moelbert@gmx.de)
getInfoPValues <- function(pvalues) {

    if (is.na(pvalues[1])) infoPValues <- "not enough information"
    else if (length(pvalues) == 1) {
        infoPValues <- paste("Kolmogorov-Smirnov-Test:",
                               as.character(pvalues[1]), sep = " " )
    } else{
        infoPValues1 <- paste("Kolmogorov-Smirnov-Test:",
                                as.character(pvalues[1]), sep = " " )
        infoPValues2 <- paste("Wilcoxon-Mann-Whitney-Test: ",
                                as.character(round(pvalues[2], 6)), sep = " ")
        infoPValues <- paste(infoPValues1, infoPValues2, sep = "\n")
    }

    infoPValues <- paste("pValues:", infoPValues, sep = "\n")
    return(infoPValues)
}
