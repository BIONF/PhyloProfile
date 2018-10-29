#' @title Run PhyloProfile app
#' @export

run_phyloprofile <- function(){
	appDir <- system.file("shiny-apps", "phyloprofile",
	                       package = "phyloprofile")
	if (appDir == "") {
		stop("Could not find apps director. Try re-installing `phyloprofile`.",
				call = FALSE)
	}

	shiny::runApp(appDir,
	              launch.browser = TRUE,
	              display.mode = "normal")
}
