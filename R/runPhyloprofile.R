#' Run PhyloProfile app
#' @export
#' @param configFile Configuration file for specifying path to input files,
#' taxonomy rank and reference taxon, and some other settings
#' @return A shiny application - GUI version of PhyloProfile
#' @import BiocStyle
#' @import DT
#' @importFrom colourpicker colourInput
#' @import energy
#' @import ExperimentHub
#' @import shinyBS
#' @import yaml
#' @rawNamespace import(RCurl, except = reset)
#' @rawNamespace import(shinyjs, except = colourInput)
#' @examples
#' ?runPhyloProfile
#' \dontrun{
#' runPhyloProfile()
#' }

runPhyloProfile <- function(configFile = NULL){
    appDir <- system.file("PhyloProfile", package = "PhyloProfile")
    if (appDir == "") {
        stop(
            "Could not find apps directory. Try re-installing `PhyloProfile`.",
            call = FALSE
        )
    }
    
    .GlobalEnv$configFile <- configFile
    on.exit(rm(configFile, envir=.GlobalEnv))

    shiny::runApp(
        appDir,
        launch.browser = TRUE,
        display.mode = "normal"
    )
}
