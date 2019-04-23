#' Run PhyloProfile app
#' @export
#' @return shiny application - GUI version of PhyloProfile
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
#' ?PhyloProfile::runPhyloProfile
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
