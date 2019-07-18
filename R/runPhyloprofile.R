#' Run PhyloProfile app
#' @export
#' @return A shiny application - GUI version of PhyloProfile
#' @import BiocStyle
#' @import Cairo
#' @importFrom colourpicker colourInput
#' @import energy
#' @import ExperimentHub
#' @import scales
#' @import shinyBS
#' @rawNamespace import(shinyjs, except = colourInput)
#' @import shinycssloaders
#' @rawNamespace import(svMisc, except = package)
#' @import tidyr
#' @import tree
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
