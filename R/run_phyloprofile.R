#' Startup script for running PhyloProfile
#' @export

run_phyloprofile <- function(){
	app_dir <- system.file("shiny-apps", "phyloprofile", 
	                       package = "phyloprofile")
	if (app_dir == "") {
		stop("Could not find apps director. Try re-installing `phyloprofile`.",
				call = FALSE)
	}

	shiny::runApp(appDir = getwd(), 
	              launch.browser = TRUE, 
	              display.mode = "normal")
}
