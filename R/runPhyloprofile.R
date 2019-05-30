#' Run PhyloProfile app
#' @export
#' @return A shiny application - GUI version of PhyloProfile
#' @import bioDist
#' @import BiocStyle
#' @importFrom colourpicker colourInput
#' @import dendextend
#' @import devtools
#' @import DT
#' @import energy
#' @import GenomeInfoDbData
#' @import gplots
#' @import GO.db
#' @import gtable
#' @import knitr
#' @rawNamespace import(rmarkdown, except = c(pdf_document, md_document,
#'     html_document))
#' @rawNamespace import(RCurl, except = c(reset, complete))
#' @import scales
#' @import shinyBS
#' @rawNamespace import(shinyjs, except = colourInput)
#' @import shinycssloaders
#' @import svMisc
#' @import tidyr
#' @examples
#' ?runPhyloProfile
#' \dontrun{
#' runPhyloProfile()
#' }

runPhyloProfile <- function(){
    appDir <- system.file("PhyloProfile", package = "PhyloProfile")
    if (appDir == "") {
        stop(
            "Could not find apps director. Try re-installing `PhyloProfile`.",
            call = FALSE
        )
    }

    shiny::runApp(
        appDir,
        launch.browser = TRUE,
        display.mode = "normal"
    )
}
