#' @title Run PhyloProfile app
#' @export
#' @return shiny application of PhyloProfile
#' @import bioDist
#' @import Cairo
#' @importFrom colourpicker colourInput
#' @import dendextend
#' @import devtools
#' @import DT
#' @import energy
#' @import GenomeInfoDbData
#' @import gplots
#' @import GO.db
#' @import gtable
#' @import scales
#' @import shinyBS
#' @import shinycssloaders
#' @import tidyr
#' @import knitr
#' @rawNamespace import(rmarkdown, except = c(pdf_document, md_document, 
#'     html_document))
#' @import BiocStyle
#' @rawNamespace import(shinyjs, except = colourInput)
#' @examples
#' ?phyloprofile::runPhyloprofile
#' \dontrun{
#' runPhyloprofile()
#' }

runPhyloprofile <- function(){
    appDir <- system.file("phyloprofile", package = "phyloprofile")
    if (appDir == "") {
        stop(
            "Could not find apps director. Try re-installing `phyloprofile`.",
            call = FALSE
        )
    }

    shiny::runApp(
        appDir,
        launch.browser = TRUE,
        display.mode = "normal"
    )
}
