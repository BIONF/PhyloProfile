#' Run PhyloProfile app
#' @export
#' @return A shiny application - GUI version of PhyloProfile
#' @import BiocStyle
#' @import Cairo
#' @importFrom colourpicker colourInput
#' @import energy
#' @import ExperimentHub
#' @import shinyBS
#' @rawNamespace import(shinyjs, except = colourInput)
#' @import shinycssloaders
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
