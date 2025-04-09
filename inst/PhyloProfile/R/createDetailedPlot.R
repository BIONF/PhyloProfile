#' Detailed plot
#' @param data (from reactive fn "detailPlotDt")
#' @param var1ID name of variable 1 (from input$var1ID)
#' @param var2ID name of variable 2 (from input$var2ID)
#' @param detailedText text size (from input$detailedText)
#' @param detailedHeight plot height (from input$detailedHeight)
#' @return information of selected protein on detailed plot
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

createDetailedPlotUI <- function(id) {
    ns <- NS(id)
    tagList(
        uiOutput(ns("detailPlot.ui")),
        downloadButton(ns("downloadDetailed"), "Download plot",
                       class = "butDL"),
        tags$head(
            tags$style(HTML(
                ".butDL{background-color:#476ba3;} .butDL{color: white;}"))
        ),
        hr(),
        verbatimTextOutput(ns("detailClick"))
    )
}

createDetailedPlot <- function(
    input, output, session, data, var1ID, var2ID, detailedText, detailedHeight,
    font
){
    # simplify ortho IDs if they are in BIONF format ---------------------------
    plotData <- reactive({
        plotDf <- data()
        if (
            checkBionfFormat(
                plotDf$orthoID[1], plotDf$geneID[1],
                gsub('ncbi','',plotDf$abbrName[1])
            )
        ) {
            plotDf <- within(
                plotDf,
                orthoMod <- data.frame(
                    do.call(
                        'rbind', strsplit(as.character(orthoID),'|', fixed=TRUE)
                    )
                )
            )
            plotDf$orthoID <- plotDf$orthoMod$X3
        }
        return(plotDf)
    })

    # render detailed plot -----------------------------------------------------
    output$detailPlot <- renderPlot({
        detailPlot(plotData(), detailedText(), var1ID(), var2ID(), font())
    })

    output$detailPlot.ui <- renderUI({
        ns <- session$ns
        shinycssloaders::withSpinner(
            plotOutput(
                ns("detailPlot"),
                width = 800,
                height = detailedHeight(),
                click = ns("plotClickDetail")
            )
        )
    })

    output$downloadDetailed <- downloadHandler(
        filename = function() {
            c("detailedPlot.svg")
        },
        content = function(file) {
            g <- detailPlot(data(), detailedText(), var1ID(), var2ID(), font())
            ggsave(
                file,
                plot = g,
                width = 800 * 0.035,
                height = detailedHeight() * 0.035,
                units = "cm",
                dpi = 300,
                device = "svg",
                limitsize = FALSE
            )
        }
    )

    # get info of detailed plot ------------------------------------------------
    pointInfoDetail <- reactive({
        selDf <- data()
        selDf$orthoID <- as.character(selDf$orthoID)

        # if only one ortholog, get directly from data()
        if (nrow(selDf) == 1) {
            seedID <- as.character(selDf$geneID)
            orthoID <- as.character(selDf$orthoID)
            taxID <- as.character(selDf$abbrName)
        }
        # else, get from plotClickDetail
        else {
            if (is.null(input$plotClickDetail$x)) return(NULL)
            else {
                corX <- round(input$plotClickDetail$y)
                corY <- round(input$plotClickDetail$x)
            }
            # get pair of sequence IDs
            seedID <- as.character(selDf$geneID[!is.na(selDf$geneID)][1])
            orthoID <- as.character(selDf$orthoID[corX])
            taxID <- as.character(selDf$abbrName[corX])
        }

        # get var1, var2
        var1 <- as.list(selDf$var1[selDf$orthoID == orthoID])
        var1 <- as.character(var1[!is.na(var1)])
        var2 <- as.list(selDf$var2[selDf$orthoID == orthoID])
        var2 <- as.character(var2[!is.na(var2)])
        if (length(var2) == 0) var2 = "NA"

        ncbiID <- selDf[selDf$orthoID == orthoID, ]$abbrName
        ncbiID <- as.character(ncbiID[!is.na(ncbiID)][1])

        # return info
        if (is.na(orthoID))
            return(NULL)
        else {
            if (orthoID != "NA")
                return(c(seedID, orthoID, var1, var2, ncbiID))
        }
    })

    # * show info when clicking on detailed plot -------------------------------
    output$detailClick <- renderText({
        info <- pointInfoDetail() # info = seedID, orthoID, var1

        if (is.null(info)) paste("select ortholog")
        else {
            a <- paste0("seedID = ", info[1])
            b <- paste0("hitID = ", info[2])
            c <- ""
            if (var1ID() != "") {
                c <- paste0(var1ID(), " = ", info[3])
            }
            d <- ""
            if (var2ID() != "") {
                d <- paste0(var2ID(), " = ", info[4])
            }
            paste(a, b, c, d, sep = "\n")
        }
    })
    return(pointInfoDetail)
}

#' create detailed plot
#' @param selDf data for plotting  (from reactive fn "detailPlotDt")
#' @param detailedText text size (input$detailedText)
#' @param var1ID name of variable 1 (input$var1ID)
#' @param var2ID name of variable 2 (input$var2ID)
#' @param font font of text. Default = Arial"
#' @return detailed plot (ggplot object)
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

detailPlot <- function(selDf, detailedText, var1ID, var2ID, font = "Arial"){
    selDf$xLabel <- paste(selDf$orthoID, " (", selDf$fullName, ")", sep = "")

    # create joined DF for plotting var1 next to var2
    var1Df <- subset(selDf, select = c("xLabel", "var1"))
    var1Df$type <- var1ID
    colnames(var1Df) <- c("id", "score", "var")

    var2Df <- subset(selDf, select = c("xLabel", "var2"))
    var2Df$type <- var2ID
    colnames(var2Df) <- c("id", "score", "var")

    detailedDf <- rbind(var1Df, var2Df)
    if(any(is.na(detailedDf$score))) {
        detailedDf[
            detailedDf$id %in% detailedDf$id[is.na(detailedDf$score)],
        ]$score <- NA
    }
    
    # remove ONE missing variable
    if (nlevels(as.factor(detailedDf$var)) > 1) {
        detailedDf <- detailedDf[nchar(detailedDf$var) > 0, ]
    }

    # keep order of ID (xLabel)
    detailedDf$id <- factor(detailedDf$id, levels = unique(detailedDf$id))

    # create plot
    gp <- ggplot(detailedDf, aes(y = score, x = id, fill = var)) +
        geom_bar(stat = "identity", position = position_dodge(), na.rm = TRUE) +
        coord_flip() +
        labs(x = "") +
        labs(fill = "") +
        theme_minimal()
    gp <- gp + theme(axis.text.x = element_text(angle = 90, hjust = 1),
                     axis.text = element_text(size = detailedText),
                     axis.title = element_text(size = detailedText),
                     legend.text = element_text(size = detailedText),
                     text = element_text(family = font)
    )
    return(gp)
}
