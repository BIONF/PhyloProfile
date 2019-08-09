#' Protein domain architecture plot
#' @param pointInfo() info of clicked point
#' (from reactive fn "pointInfoDetail")
#' @param domainInfo() domain information
#' (from reactive fn "getDomainInformation")
#' @param labelArchiSize lable size (from input$labelArchiSize)
#' @param titleArchiSize title size (from input$titleArchiSize)
#' @param archiHeight plot height (from input$archiHeight)
#' @param archiWidth plot width (from input$archiWidth)
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

source("R/functions.R")

createArchitecturePlotUI <- function(id) {
    ns <- NS(id)
    tagList(
        uiOutput(ns("archiPlot.ui")),
        downloadButton(ns("archiDownload"), "Download plot", class = "butDL"),
        tags$head(
            tags$style(HTML(
                ".butDL{background-color:#476ba3;} .butDL{color: white;}"))
        ),
        textOutput(ns("selectedDomain"))
    )
}

createArchitecturePlot <- function(
    input, output, session,
    pointInfo, domainInfo,
    labelArchiSize, titleArchiSize, archiHeight, archiWidth
){
    output$archiPlot <- renderPlot({
        if (is.null(nrow(domainInfo()))) stop("Domain info is NULL!")
        g <- createArchiPlot(
            pointInfo(), domainInfo(), labelArchiSize(), titleArchiSize()
        )
        if (any(g == "No domain info available!")) {
            msgPlot()
        } else {
            grid.draw(g)
        }
    })

    output$archiPlot.ui <- renderUI({
        ns <- session$ns
        if (is.null(nrow(domainInfo()))) {
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
                    height = archiHeight(),
                    width = archiWidth(),
                    click = ns("archiClick")
                )
            )
        }
    })

    output$archiDownload <- downloadHandler(
        filename = function() {
            c("domains.pdf")
        },
        content = function(file) {
            g <- createArchiPlot(
                pointInfo(), domainInfo(), labelArchiSize(), titleArchiSize()
            )
            grid.draw(g)
            ggsave(
                file, plot = g,
                width = archiWidth() * 0.056458333,
                height = archiHeight() * 0.056458333,
                units = "cm", dpi = 300, device = "pdf", limitsize = FALSE
            )
        }
    )

    output$selectedDomain <- renderText({
        if (is.null(input$archiClick$y)) return("No domain selected!")
        convertY(unit(input$archiClick$y, "npc"), "native")
    })
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
