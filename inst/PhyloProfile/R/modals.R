
# * popup for plotting detailed plot -------------------------------------------
detailedPlotModal <- function() {
    modalDialog(
        title = "Detailed plot",
        size = "l",
        easyClose = TRUE,
        footer = modalButton("Close"),
        
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
        shinyBS::bsButton(
            "doDomainPlot", "Show domain architecture"#, disabled = TRUE
        ),
        uiOutput("checkDomainFiles"),
        br(),
        h4("Sequence:"),
        verbatimTextOutput("fasta"),
        br(),
        h4("Links:"),
        uiOutput("dbLink.ui")
    )
}

# * popup for plotting domain architecture plot --------------------------------
archiPlotModal <- function(id, title = "Domain architecture") {
    modalDialog(
        title = title,
        size = "l",
        easyClose = TRUE,
        footer = modalButton("Close"),
        createArchitecturePlotUI("archiPlot")
    )
}

# * popup for confirming parsing taxa from input file ------------------
parseConfirmModal <- function() {
    modalDialog(
        title = "Get taxonomy info",
        size = "s",
        easyClose = FALSE,
        footer = NULL,
        HTML(
            '<p>Fetching Missing Taxonomy Information and
      Post-processing.</p><p><em>This window will close
      automatically when everything is done. Please wait...</em></p>
      <p><strong><span style="color: #ff0000;">PLEASE RELOAD THIS
      TOOL WHEN FINISHED!!!</span></strong></p>'
        )
    )
}