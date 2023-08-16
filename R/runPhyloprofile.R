#' Run PhyloProfile app
#' @export
#' @param configFile Configuration file for specifying path to input files,
#' taxonomy rank and reference taxon, and some other settings
#' @param host IP adress (e.g. host = "127.0.0.1")
#' @param port Port (e.g. port = 8888)
#' @return A shiny application - GUI version of PhyloProfile
#' @import BiocStyle
#' @import DT
#' @importFrom colourpicker colourInput
#' @import energy
#' @import ExperimentHub
#' @import shinyBS
#' @import shinycssloaders
#' @import shinyFiles
#' @import stringr
#' @import yaml
#' @rawNamespace import(RCurl, except = reset)
#' @rawNamespace import(shinyjs, except = colourInput)
#' @examples
#' ?runPhyloProfile
#' \dontrun{
#' runPhyloProfile()
#' }

runPhyloProfile <- function(configFile = NULL, host = NULL, port = NULL){
    appDir <- system.file("PhyloProfile", package = "PhyloProfile")
    if (appDir == "") {
        stop(
            "Could not find apps directory. Try re-installing `PhyloProfile`.",
            call = FALSE
        )
    }

    .GlobalEnv$configFile <- configFile
    on.exit(rm(configFile, envir=.GlobalEnv))
    
    i_host <- i_port <- NULL
    i_launchBrowser <- TRUE
    if (!is.null(configFile) && file.exists(configFile)) {
        configs <- yaml::read_yaml(configFile)
        i_host <- configs$host
        i_port <- configs$port
        i_launchBrowser <- configs$launchBrowser
    }
    if (!is.null(host)) i_host <- host
    if (!is.null(port)) i_port <- port
    
    if (!is.logical(i_launchBrowser)) i_launchBrowser <- TRUE
    if (!is.null(i_host) && !is.null(i_port)) {
        shiny::runApp(
            appDir,
            host = i_host, port = i_port, launch.browser = i_launchBrowser,
            display.mode = "normal"
        )
    } else {
        shiny::runApp(
            appDir,
            launch.browser = TRUE,
            display.mode = "normal"
        )
    }
}
