#' Group Comparison
#' @param filteredDf phylogenetic profile data after filtering
#' @param inGroup list of taxon IDs for the in-group taxa
#' @param variable selected variable for calculating the differences betwee
#' in- and out-group ("var1" or "var2")
#' @param varName names of the additional variables
#' @param compareType type of comparison (statistical tests or mean comparison)
#' @param significanceLevel confidence level for stat tests or delta mean cutoff
#' @param plotParameters plot parameters including xSize, ySize, angle (x-axis
#' angle), legendPosition and legendSize, titleSize, flipPlot, mValue (mean or
#' median point), widthVar and heightVar (width and height for variable
#' distribution plot), widthFeature and heightFeature (width and height for
#' feature distribution plot), inGroupName and outGroupName.
#' @param domainDf dataframe of domain architectures
#' @param doCompare check if Compare button clicked (input$doCompare)
#' @param doUpdate check if Update button clicked (input$updateGC)
#' @return list of candidate genes and their p-values (or delta means)
#' @authors Carla Mölbert {carla.moelbert@gmx.de}, Vinh Tran
#' {tran@bio.uni-frankfurt.de}

groupComparisonUI <- function(id){
    ns <- NS(id)
    sidebarLayout(
        sidebarPanel(
            width = 3,
            uiOutput(ns("candidateGenes.ui")),
            shinyBS::bsPopover(
                "candidateGenes.ui",
                "",
                "Select gene to show the plots",
                "right"
            ),

            downloadButton(ns("downloadPlot"),"Download plot",
                           class = "butDL"),
            tags$head(
                tags$style(HTML(
                    ".butDL{background-color:#476ba3;} .butDL{color: white;}"))
            ),
            hr(),
            uiOutput(ns("featureTypeSelect.ui")),
            selectInput(
                ns("filterGainLoss"),
                "Filter gain/loss features",
                choices = c("gain", "loss", "none"),
                selected = "none"
            ),
            numericInput(
                ns("featureThreshold"),
                "Number of intances per protein (IPP) cutoff",
                min = 0, max = 9999, step = 0.1, value = 0.25
            ),
            sliderInput(
                ns("dIPPThreshold"),
                "Delta IPP / Sum IPP cutoff",
                min = 0 , max = 1, value = 0,
                step = 0.01, round = FALSE
            )
        ),
        mainPanel(
            width = 9,
            tags$style(
                HTML("#plotsUI { height:650px; overflow-y:scroll}")
            ),
            uiOutput(ns("downloadALl.ui")),
            uiOutput(ns("varPlots.ui")),
            uiOutput(ns("featurePlots.ui"))
        )
    )
}

groupComparison <- function (
    input, output, session,
    filteredDf,
    inGroup,
    variable,
    varName,
    compareType,
    significanceLevel,
    plotParameters,
    domainDf,
    doCompare,
    doUpdate
) {
    ### get candidate genes and their p-values
    candidateGenes <- reactive({
        if (is.null(inGroup())) return()
        if (is.null(variable()) | variable()[1] == "none") 
            stop("No variable available or selected!")

        if (compareType() == "Statistical tests") {
            pvalues <- compareTaxonGroups(
                filteredDf(),
                inGroup(),
                FALSE,
                variable(),
                significanceLevel()
            )
            return(pvalues[pvalues <= significanceLevel()])
        } else {
            deltaMean <- compareMedianTaxonGroups(
                filteredDf(),
                inGroup(),
                FALSE,
                variable()
            )
            return(deltaMean[deltaMean <= significanceLevel()])
        }
    })

    ### render list of candidate genes
    output$candidateGenes.ui <- renderUI({
        ns <- session$ns
        doCompare()
        isolate({
            candidateGenes <- names(candidateGenes())

            if (length(candidateGenes()) > 0) {
                candidateGenesList <- c("all", candidateGenes)
                selectInput(
                    ns("candidateGenes"), "Candidate gene(s):",
                    choices = candidateGenesList,
                    selected = candidateGenesList[2],
                    multiple = FALSE
                )
            } else {
                selectInput(ns("candidateGenes"), "Candidate gene(s):", "none")
            }
        })
    })

    ### generate data for variable distribution plotting
    plotDf <- reactive({
        if (is.null(candidateGenes())) return()
        if (length(input$candidateGenes) == 0 | input$candidateGenes == "none")
            return()

        if (input$candidateGenes[1] == "all") return()
        return(
            getVarDistributionData(
                filteredDf(), inGroup(), input$candidateGenes, varName()
            )
        )
    })

    ### render variable distribution plot(s)
    output$distributionPlot <- renderPlot({
        if (is.null(plotDf())) return()
        g <- plotVarDistribution(
            input$candidateGenes, candidateGenes(), compareType(), plotDf(),
            plotParameters()
        )
        grid.newpage()
        grid.draw(g)
    })

    output$varPlots.ui <- renderUI({
        if (is.null(plotDf())) return()
        ns <- session$ns
        plotOutput(
            ns("distributionPlot"),
            width = plotParameters()$widthVar,
            height = plotParameters()$heightVar
        )
    })

    ### get data for protein features plotting
    featureDf <- reactive({
        doUpdate()
        if (is.null(domainDf()) | is.null(plotDf())) return()
        if (length(input$candidateGenes) == 0 | input$candidateGenes == "none")
            return()

        featureDf <- dataFeatureTaxGroup(
            filteredDf(), domainDf(), inGroup(), input$candidateGenes
        )
        ippCutoff <- isolate(input$featureThreshold)
        dIPPCutoff <- isolate(input$dIPPThreshold)

        if (input$filterGainLoss == "loss") {
            loss <- featureDf$Feature[featureDf$Taxon_group == "Out-group"
                                      & featureDf$IPP >= ippCutoff]
            return(
                featureDf[featureDf$Feature %in% loss
                          & featureDf$dIPP < 0
                          & abs(featureDf$dIPP) >= dIPPCutoff,]
            )
        } else if (input$filterGainLoss == "gain") {
            gain <- featureDf$Feature[featureDf$Taxon_group == "In-group"
                                      & featureDf$IPP >= ippCutoff]
            return(
                featureDf[featureDf$Feature %in% gain
                          & featureDf$dIPP >= dIPPCutoff,]
            )
        } else {
            dfIn <- featureDf[featureDf$Taxon_group == "In-group",]
            dfOut <- featureDf[featureDf$Taxon_group == "Out-group",]
            mergedDf <- merge(dfIn, dfOut, by = "Feature", all = TRUE)
            mergedDf[is.na(mergedDf)] <- 0
            allFeature <- as.character(
                mergedDf$Feature[
                    mergedDf$IPP.x >= ippCutoff | mergedDf$IPP.y >= ippCutoff
                ]
            )
            return(featureDf[featureDf$Feature %in% allFeature
                             & abs(featureDf$dIPP) >= dIPPCutoff,])
        }
    })

    ### List with possible features for the selected gene
    output$featureTypeSelect.ui <- renderUI({
        ns <- session$ns
        if (is.null(featureDf())) {
            selectInput(ns("featureTypeSelect"),
                        "Feature type(s) of interest",
                        choices = "ALL",
                        selected = "ALL",
                        multiple = TRUE,
                        selectize = FALSE)
        } else {
            featureDf <- data.frame(
                do.call(
                    rbind, strsplit(as.character(featureDf()$Feature), "_")
                ), stringsAsFactors = FALSE
            )
            featureList <- unique(featureDf[,1])
            selectInput(ns("featureTypeSelect"),
                        "Feature type(s) of interest:",
                        choices = c("ALL", featureList),
                        selected = "ALL",
                        multiple = TRUE,
                        selectize = FALSE)
        }
    })

    ### filter feature data based on selected type of features
    featureDfSelected <- reactive({
        if (is.null(featureDf())) return()
        if ("ALL" %in% input$featureTypeSelect) return(featureDf())
        else {
            featureDf <- featureDf()
            selectedDf <- lapply(
                input$featureTypeSelect,
                function (x) featureDf[grep(x, featureDf$Feature),]
            )
            return(do.call(rbind, selectedDf))
        }
    })

    ### render feature distribution plot(s)
    output$featurePlot <- renderPlot({
        if (is.null(featureDfSelected())) return()
        featureDistTaxPlot(featureDfSelected(), plotParameters())
    })

    output$featurePlots.ui <- renderUI({
        if (is.null(featureDf())) return()
        ns <- session$ns
        plotOutput(
            ns("featurePlot"),
            width = plotParameters()$widthFeature,
            height = plotParameters()$heightFeature
        )
    })

    ### download plots for single candidate
    output$downloadPlot <- downloadHandler(
        filename = function() {
            paste0("GC_", input$candidateGenes, ".pdf")
        },
        content = function(file) {
            ggsave(
                file,
                plot = arrangeGrob(
                    plotVarDistribution(
                        input$candidateGenes,
                        candidateGenes(),
                        compareType(),
                        plotDf(),
                        plotParameters()
                    ),
                    featureDistTaxPlot(featureDfSelected(), plotParameters()),
                    ncol = 1
                ),
                dpi = 300, device = "pdf"
            )
        }
    )

    ### download plots for all candidates ###########
    output$downloadALl.ui <- renderUI({
        ns <- session$ns
        if (is.null(input$candidateGenes)) return()
        if (input$candidateGenes[1] == "all") {
            tagList(
                h3(em("For analyzing all candidate genes, please download them
                      to your computer! Current plot settings and all thresholds
                      will be applied!")),
                downloadButton(
                    ns("downloadAll"),"Download all plots", class = "butDL"
                )
            )
        }
    })

    output$downloadAll <- downloadHandler(
        filename = function() {
            paste0("GC_allGenes.zip")
        },
        content = function(file) {
            fs <- NULL
            for(gene in names(candidateGenes())) {
                fs <- c(fs, paste0("GC_", gene, ".pdf"))
                pdf(paste0("GC_", gene, ".pdf"))
                plot <- arrangeGrob(
                    plotVarDistribution(
                        gene,
                        candidateGenes(),
                        compareType(),
                        getVarDistributionData(
                            filteredDf(), inGroup(),
                            gene, varName()
                        ),
                        plotParameters()
                    ),
                    featureDistTaxPlot(
                        getSelectedFeatureData(
                            filteredDf(), domainDf(), inGroup(), gene,
                            input$featureThreshold, input$dIPPThreshold,
                            input$filterGainLoss, input$featureTypeSelect
                        ),
                        plotParameters()
                    ),
                    ncol = 1
                )
                grid.draw(plot)
                dev.off()
            }
            zip(zipfile = file, files = fs)
        },
        contentType = "application/zip"
    )

    ### disable/enable download functions
    observe({
        if (is.null(filteredDf())) shinyjs::disable("downloadPlot")
        else if (doCompare() == FALSE) shinyjs::disable("downloadPlot")
        else if (input$candidateGenes[1] == "all")
            shinyjs::disable("downloadPlot")
        else shinyjs::enable("downloadPlot")
    })
    
    ### return list of candidate genes
    outGenes <- reactive({
        return(names(candidateGenes()))
    })
    if (!is.null(outGenes)) return(outGenes)
}

################################# FUNCTIONS ####################################

#' Get data for variable distribution plot
#' @param data main (filtered) phylogenetic profile dataframe
#' @param inGroup list of in-group taxa (ncbi ID)
#' @param candidateGene selected candidate gene
#' @param varName names of up to 2 variables
#' @return dataframe contains values of the variables that are used for plotting
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
getVarDistributionData <- function(
    data, inGroup, candidateGene, varName
) {
    names(data)[names(data)=="var1"] <- varName[1]
    if (length(varName) == 2)
        names(data)[names(data)=="var2"] <- varName[2]
    return(
        dataVarDistTaxGroup(data, inGroup, candidateGene, varName)
    )
}

#' Plot variable distribution
#' @param candidateGene name of candidate gene
#' @param candidateGenesList named candidate gene list (used to get p-value or
#' delta median for selected candidate gene)
#' @param compareType "Statistical tests" or "Median values"
#' @param plotDf Data for plotting
#' @param plotParameters plot parameters
#' @return grid object. Use grid.draw() to plot.
#' @authors Vinh Tran {tran@bio.uni-frankfurt.de}
plotVarDistribution <- function(
    candidateGene, candidateGenesList, compareType, plotDf, plotParameters
) {
    pValue <- candidateGenesList[
        match(candidateGene, names(candidateGenesList))]
    if (compareType == "Statistical tests") {
        title = paste0(candidateGene, "   p-value = ", round(pValue, 3))
    } else {
        title <- paste0(candidateGene, "   delta-mean = ", round(pValue, 3))
    }

    g <- varDistTaxPlot(plotDf, c(plotParameters, title = title))
    return(g)
}

#' Get data for feature comparison plot
#' @param mainData main (filtered) phylogenetic profile dataframe
#' @param domainDf dataframe of domain architectures
#' @param inGroup list of in-group taxa (ncbi ID)
#' @param candidateGene name of candidate gene
#' @param ippCutoff cutoff for number of intances per protein
#' @param dIPPCutoff cutoff for the differences in the IPP of the same feature
#' between in- and out-group
#' @param filterGainLoss get only possible gain/loss features
#' @param featureTypeSelect type of features (pfam/smart/cast/etc.)
#' @return dataframe for feature comparison plot
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
getSelectedFeatureData <- function(
    mainData, domainDf, inGroup, candidateGene, ippCutoff, dIPPCutoff,
    filterGainLoss, featureTypeSelect
) {
    featureDf <- dataFeatureTaxGroup(
        mainData, domainDf, inGroup, candidateGene
    )

    if (filterGainLoss == "loss") {
        loss <- featureDf$Feature[featureDf$Taxon_group == "Out-group"
                                  & featureDf$IPP >= ippCutoff]
        featureDfOut <- featureDf[
            featureDf$Feature %in% loss & featureDf$dIPP < 0
            & abs(featureDf$dIPP) >= dIPPCutoff,
        ]
    } else if (filterGainLoss == "gain") {
        gain <- featureDf$Feature[featureDf$Taxon_group == "In-group"
                                  & featureDf$IPP >= ippCutoff]
        featureDfOut <- featureDf[
            featureDf$Feature %in% gain & featureDf$dIPP >= dIPPCutoff,
        ]
    } else {
        dfIn <- featureDf[featureDf$Taxon_group == "In-group",]
        dfOut <- featureDf[featureDf$Taxon_group == "Out-group",]
        mergedDf <- merge(dfIn, dfOut, by = "Feature", all = TRUE)
        mergedDf[is.na(mergedDf)] <- 0
        allFeature <- as.character(
            mergedDf$Feature[
                mergedDf$IPP.x >= ippCutoff | mergedDf$IPP.y >= ippCutoff
                ]
        )
        featureDfOut <- featureDf[
            featureDf$Feature %in% allFeature
            & abs(featureDf$dIPP) >= dIPPCutoff,
        ]
    }

    ### filter data based on selected feature types
    if ("ALL" %in% featureTypeSelect) return(featureDfOut)
    selectedDf <- lapply(
        featureTypeSelect,
        function (x) featureDfOut[grep(x, featureDfOut$Feature),]
    )
    return(do.call(rbind, selectedDf))
}
