#' Function to keep user defined geneID order
#' @export
#' @param data data frame contains gene ID column
#' @param order TRUE or FALSE (from input$ordering)
#' @return data either sorted or non-sorted
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
unsort_id <- function(data, order){
    data$geneID <- as.factor(data$geneID)
    if (order == FALSE) {
        # keep user defined geneID order
        data$geneID <- factor(data$geneID, levels = unique(data$geneID))
    }
    return(data)
}

#' Check installed packages
#' and install missing packages automatically
#' @param packages list of packages need to be checked
#' @return none
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

install_packages <- function(packages){
    missing_packages <-
        packages[!(packages %in% installed.packages()[, "Package"])]
    if (length(missing_packages))
        install.packages(
            missing_packages,
            dependencies = TRUE,
            repos = "http://cran.us.r-project.org"
        )
}

install_packages_bioconductor <- function(packages){
    missing_packages <-
        packages[!(packages %in% installed.packages()[, "Package"])]
    if (length(missing_packages)) {
        source("https://bioconductor.org/biocLite.R")
        biocLite(missing_packages, ask = FALSE)
    }
}

#' Check internet connection
#' @return status of internet connection
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

has_internet <- function(){
    !is.null(curl::nslookup("r-project.org", error = FALSE))
}

#' #' Get last n characters from string x
#' substr_right <- function(x, n){
#'   substr(x, nchar(x) - n + 1, nchar(x))
#' }

# FUNCTIONS FOR RENDER UI ELEMENTS ============================================
create_slider_cutoff <- function(id, title, start, stop, var_id){
    if (is.null(var_id)) return()
    if (var_id == "") {
        sliderInput(id, title,
                    min = 1,
                    max = 1,
                    step = 0.025,
                    value = 1,
                    width = 200)
    } else {
        sliderInput(id, title,
                    min = 0,
                    max = 1,
                    step = 0.025,
                    value = c(start, stop),
                    width = 200)
    }
}

update_slider_cutoff <- function(session, id, title, new_var, var_id){
    if (is.null(var_id) || var_id == "") return()

    updateSliderInput(session, id, title,
                      value = new_var,
                      min = 0,
                      max = 1,
                      step = 0.025)
}

create_plot_size <- function(id, title, value) {
    numericInput(id,
                 title,
                 min = 100,
                 max = 3200,
                 step = 50,
                 value = value,
                 width = 100)
}

create_text_size <- function(id, title, value, width) {
    numericInput(id,
                 title,
                 min = 3,
                 max = 99,
                 step = 1,
                 value = value,
                 width = width)
}

create_select_gene <- function(id, list, selected) {
    selectInput(id,
                "",
                list,
                selected = selected,
                multiple = TRUE,
                selectize = FALSE)
}
