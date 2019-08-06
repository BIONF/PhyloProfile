#' Profile plot
#' @param data data for heatmap plot (from reactive fn "dataHeat")
#' @param clusteredDataHeat clustered data (from reactive fn "clusteredDataHeat"
#' @param applyCluster choose clustered data or not (from input$applyCluster)
#' @param parameters plot parameters (colors, size, variable names, ...)
#' @param inSeq subset sequences for customized profile (input$inSeq)
#' @param inTaxa subset taxa for customized profile (input$inTaxa)
#' @param rankSelect selected taxonomy rank (input$rankSelect)
#' @param inSelect selected taxon name (input$inSelect)
#' @param taxonHighlight highlighted taxon (input$taxonHighlight)
#' @param geneHighlight highlighted gene (input$geneHighlight)
#' @param typeProfile either "mainProfile" or "customizedProfile"
#' @return info for selected point on the profile
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

source("R/functions.R")

createProfilePlotUI <- function(id) {
    ns <- NS(id)
    tagList(
        uiOutput(ns("plot.ui")),
        br(),
        downloadButton(ns("profileDownload"),"Download profile",
                       class = "butDL"),
        tags$head(
            tags$style(HTML(
                ".butDL{background-color:#476ba3;} .butDL{color: white;}"))
        )
    )
}

createProfilePlot <- function(input, output, session,
                                data, clusteredDataHeat,
                                applyCluster,
                                parameters,
                                inSeq, inTaxa,
                                rankSelect, inSelect,
                                taxonHighlight, geneHighlight,
                                typeProfile,
                                colorByGroup) {
    # data for heatmap ---------------------------------------------------------
    dataHeat <- reactive({
        if (is.null(data())) stop("Profile data is NULL!")

        if (typeProfile() == "customizedProfile") {
            if (is.null(inTaxa()) | is.null(inSeq())) return()

            dataHeat <- dataCustomizedPlot(data(), inTaxa(), inSeq())
            if (applyCluster() == TRUE) {
                dataHeat <- dataCustomizedPlot(
                    clusteredDataHeat(), inTaxa(), inSeq()
                )
            }
        } else {
            dataHeat <- dataMainPlot(data())
            if (applyCluster() == TRUE) {
                dataHeat <- dataMainPlot(clusteredDataHeat())
            }
        }
        return(dataHeat)
    })

    # render heatmap profile ---------------------------------------------------
    output$plot <- renderPlot({
        if (is.null(data())) stop("Profile data is NULL!")
        if (typeProfile() == "customizedProfile") {
            if (inSeq()[1] == "all" & inTaxa()[1] == "all") return()
        }

        highlightProfilePlot(
            dataHeat(),
            parameters(),
            taxonHighlight(),
            rankSelect(),
            geneHighlight()
        )
    })

    output$plot.ui <- renderUI({
        ns <- session$ns

        if (typeProfile() == "customizedProfile") {
            if (is.null(inSeq()[1]) | is.null(inTaxa()[1])) return()
            else if (inSeq()[1] == "all" & inTaxa()[1] == "all")  return()
        }

        shinycssloaders::withSpinner(
            plotOutput(
                ns("plot"),
                width = parameters()$width,
                height = parameters()$height,
                click = ns("plotClick")
            )
        )
    })

    output$profileDownload <- downloadHandler(
        filename = function() {
            c("profile.pdf")
        },
        content = function(file) {
            ggsave(
                file,
                plot = highlightProfilePlot(
                    dataHeat(), parameters(), "none", rankSelect(), "none"
                ),

                width = parameters()$width * 0.056458333,
                height = parameters()$height * 0.056458333,
                units = "cm", dpi = 300, device = "pdf", limitsize = FALSE
            )
        }
    )
    # get info of clicked point on heatmap plot --------------------------------
    selectedpointInfo <- reactive({

        # get selected supertaxon name
        taxaList <- getNameList()
        rankName <- rankSelect()
        inSelect <- taxaList$ncbiID[taxaList$fullName == inSelect()]

        dataHeat <- dataHeat()
        if (is.null(dataHeat)) stop("Data for heatmap is NULL!")

        if (typeProfile() == "customizedProfile") {
            # get sub-dataframe of selected taxa and sequences
            dataHeat$supertaxonMod <- substr(
                dataHeat$supertaxon,
                6,
                nchar(as.character(dataHeat$supertaxon))
            )

            if (is.null(inSeq()[1]) | is.null(inTaxa()[1]))  
                stop("Subset taxa or genes is NULL!")
            if (inTaxa()[1] == "all" & inSeq()[1] != "all") {
                # select data from dataHeat for selected sequences only
                dataHeat <- subset(dataHeat, geneID %in% inSeq())
            } else if (inSeq()[1] == "all" & inTaxa()[1] != "all") {
                # select data from dataHeat for selected taxa only
                dataHeat <- subset(dataHeat, supertaxonMod %in% inTaxa())
            } else {
                # select data from dataHeat for selected sequences and taxa
                dataHeat <- subset(dataHeat, geneID %in% inSeq()
                                   & supertaxonMod %in% inTaxa())
            }

            # drop all other supertaxon that are not in sub-dataframe
            dataHeat$supertaxon <- factor(dataHeat$supertaxon)
            dataHeat$geneID <- factor(dataHeat$geneID)
        }

        # get values
        if (is.null(input$plotClick$x)) return()
        else {
            # get cooridiate point
            if (parameters()$xAxis == "genes") {
                corX <- round(input$plotClick$y);
                corY <- round(input$plotClick$x)
            } else {
                corX <- round(input$plotClick$x);
                corY <- round(input$plotClick$y)
            }

            # get geneID
            genes <- levels(dataHeat$geneID)
            geneID <- toString(genes[corY])
            # get supertaxon (spec)
            supertaxa <- levels(dataHeat$supertaxon)
            spec <- toString(supertaxa[corX])
            # get var1, percentage of present species and var2 score
            var1 <- NA
            if (!is.na(dataHeat$var1[dataHeat$geneID == geneID
                                     & dataHeat$supertaxon == spec][1])) {
                var1 <- max(
                    na.omit(dataHeat$var1[dataHeat$geneID == geneID
                                          & dataHeat$supertaxon == spec])
                )
            }
            Percent <- NA
            if (!is.na(dataHeat$presSpec[dataHeat$geneID == geneID
                                         & dataHeat$supertaxon == spec][1])) {
                Percent <- {
                    max(
                        na.omit(
                            dataHeat$presSpec[dataHeat$geneID == geneID
                                              & dataHeat$supertaxon == spec]
                        )
                    )
                }
            }
            var2 <- NA
            if (!is.na(dataHeat$var2[dataHeat$geneID == geneID
                                     & dataHeat$supertaxon == spec][1])) {
                var2 <- {
                    max(na.omit(dataHeat$var2[dataHeat$geneID == geneID
                                              & dataHeat$supertaxon == spec]))
                }
            }

            # get ortholog ID
            orthoID <- dataHeat$orthoID[dataHeat$geneID == geneID
                                        & dataHeat$supertaxon == spec]
            if (length(orthoID) > 1) {
                orthoID <- paste0(orthoID[1], ",...")
            }

            if (is.na(Percent)) return()
            else {
                info <- c(geneID,
                          as.character(orthoID),
                          spec,
                          round(var1, 2),
                          round(Percent, 2),
                          round(var2, 2))
                return(info)
            }
        }
    })

    return(selectedpointInfo)
}
