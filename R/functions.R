require(tictoc)

#' Calculate percentage of present species-------------------------------------
#' @export
#' @param taxa_md_data contains "geneID", "ncbiID", "orthoID",
#' "var1", "var2", "paralog", ...., and "supertaxon"
#' @return new data frame with % of present species
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

calc_pres_spec <- function(taxa_md_data, taxa_count){
  taxa_md_data <- taxa_md_data[taxa_md_data$orthoID != "NA", ]
  
  # get geneID and supertaxon
  gene_id_supertaxon <- subset(taxa_md_data,
                               select = c("geneID",
                                          "supertaxon",
                                          "paralog",
                                          "abbrName"))
  # remove duplicated rows
  gene_id_supertaxon <- gene_id_supertaxon[!duplicated(gene_id_supertaxon), ]
  
  # remove NA rows from taxa_md_data
  taxa_md_data_no_na <- taxa_md_data[taxa_md_data$orthoID != "NA", ]
  
  # count present frequency of supertaxon for each gene
  gene_supertaxon_count <- plyr::count(taxa_md_data_no_na,
                                       c("geneID", "supertaxon"))
  
  # merge with taxa_count to get total number of species of each supertaxon
  # and calculate presSpec
  pres_spec_dt <- merge(gene_supertaxon_count,
                        taxa_count,
                        by = "supertaxon",
                        all.x = TRUE)
  
  spec_count <- plyr::count(gene_id_supertaxon, c("geneID",
                                                  "supertaxon"))
  pres_spec_dt <- merge(pres_spec_dt,
                        spec_count, by = c("geneID",
                                           "supertaxon"))
  
  pres_spec_dt$presSpec <- pres_spec_dt$freq / pres_spec_dt$freq.y
  
  pres_spec_dt <- pres_spec_dt[pres_spec_dt$presSpec <= 1, ]
  pres_spec_dt <- pres_spec_dt[order(pres_spec_dt$geneID), ]
  pres_spec_dt <- pres_spec_dt[, c("geneID", "supertaxon", "presSpec")]
  
  # add absent supertaxon into pres_spec_dt
  gene_id_supertaxon <- subset(gene_id_supertaxon,
                               select = -c(paralog, abbrName))
  final_pres_spec_dt <- merge(pres_spec_dt,
                              gene_id_supertaxon,
                              by = c("geneID", "supertaxon"),
                              all.y = TRUE)
  final_pres_spec_dt$presSpec[is.na(final_pres_spec_dt$presSpec)] <- 0
  
  # remove duplicated rows
  final_pres_spec_dt <- final_pres_spec_dt[!duplicated(final_pres_spec_dt), ]
  
  # return final_pres_spec_dt
  return(final_pres_spec_dt)
}

#' Get list with all main taxanomy ranks --------------------------------------
#' @export
#' @param none
#' @return list of all main ranks
#' @author Carla Mölbert (carla.moelbert@gmx.de)

get_taxonomy_ranks <- function(){
  all_ranks <- list("Strain " = "05_strain",
                    "Species" = "06_species",
                    "Genus" = "10_genus",
                    "Family" = "14_family",
                    "Order" = "19_order",
                    "Class" = "23_class",
                    "Phylum" = "26_phylum",
                    "Kingdom" = "28_kingdom",
                    "Superkingdom" = "29_superkingdom",
                    "unselected" = "")
  return(all_ranks)
}

#' Get name list (from "data/taxonNamesReduced.txt") --------------------------
#' @export
#' @param none
#' @return list of all main ranks
#' @author Carla Mölbert (carla.moelbert@gmx.de)

get_name_list <- function(as_character, delete_duplicated) {
  name_list <- as.data.frame(read.table("data/taxonNamesReduced.txt",
                                        sep = "\t",
                                        header = T,
                                        fill = TRUE))
  if (as_character) {
    name_list$fullName <- as.character(name_list$fullName)
    name_list$rank <- as.character(name_list$rank)
  }
  if (delete_duplicated) {
    name_list <- name_list[!duplicated(name_list), ]
  }
  return(name_list)
}

#' Get taxa list (from "data/taxonomyMatrix.txt") -----------------------------
#' @export
#' @param none
#' @return list of all main ranks
#' @author Carla Mölbert (carla.moelbert@gmx.de)

get_taxa_list <- function(subset_taxa_check, subset_taxa){
  dt <- as.data.frame(read.table("data/taxonomyMatrix.txt",
                                 sep = "\t",
                                 header = T,
                                 stringsAsFactors = T))
  if (subset_taxa_check) {
    dt <- dt[dt$abbrName  %in% subset_taxa, ]
  }
  return(dt)
}


#' Function to keep user defined geneID order ---------------------------------
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
    biocLite(missing_packages)
  }
}

#' Check internet connection --------------------------------------------------
#' @return status of internet connection
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

has_internet <- function(){
  !is.null(curl::nslookup("r-project.org", error = FALSE))
}

#' #' Get last n characters from string x
#' substr_right <- function(x, n){
#'   substr(x, nchar(x) - n + 1, nchar(x))
#' }

###################### FUNCTIONS FOR RENDER UI ELEMENTS #######################

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


# PROFILE CLUSTERING =========================

#' needed for:  mutual information, pearson
install_packages_bioconductor("bioDist")
library(bioDist)

#' needed for: distance correlation 
install_packages("energy")
library(energy)


#' get the phylogenetic profiles ----------------------------------------------
#' @export
#' @param data
#' @param dist_method
#' @param var1_aggregate_by
#' @param var2_aggregate_by
#' @return dataframe containing phylogenetic profiles
#' @author Carla Mölbert (carla.moelbert@gmx.de)
get_data_clustering <- function(data,
                                profile_type,
                                var1_aggregate_by,
                                var2_aggregate_by){
  sub_data_heat <- subset(data, data$presSpec > 0)
  if (profile_type == "binary") {
    #' Binary Profiles --------------------------------------------------------
    sub_data_heat <- sub_data_heat[, c("geneID", "supertaxon", "presSpec")]
    print(head(sub_data_heat,20))
    sub_data_heat$presSpec[sub_data_heat$presSpec > 0] <- 1
    print(head(sub_data_heat,20))
    sub_data_heat <- sub_data_heat[!duplicated(sub_data_heat), ]
    wide_data <- spread(sub_data_heat, supertaxon, presSpec)
  
    
  }else {
    #' Profiles with FAS scores ----------------------------------------------
    var <- profile_type
    
    sub_data_heat <- sub_data_heat[, c("geneID", "supertaxon", var)]
    sub_data_heat <- sub_data_heat[!duplicated(sub_data_heat), ]
    
    #' aggreagte the values by the selected method
    if (var == "var1") aggregate_by <- var1_aggregate_by
    else aggregate_by <- var2_aggregate_by
    
    sub_data_heat <- aggregate(sub_data_heat[, var],
                               list(sub_data_heat$geneID,
                                    sub_data_heat$supertaxon),
                               FUN = aggregate_by)
    
    colnames(sub_data_heat) <- c("geneID", "supertaxon", var)
    
    wide_data <- spread(sub_data_heat, supertaxon, var)
  }
  
  dat <- wide_data[, 2:ncol(wide_data)]  # numerical columns
  rownames(dat) <- wide_data[, 1]
  dat[is.na(dat)] <- 0
  
  return(dat)
}

#' Get the distance matrix depending on the distance method--------------------
#' @export
#' @param profiles datafram containing phylogenetic profiles
#' @param dist_method distance method
#' @return distance matrix
#' @author Carla Mölbert (carla.moelbert@gmx.de)
get_distance_matrix <- function(profiles, method, cutoff){

  dist_methods <- c("euclidean", "maximum", "manhattan", "canberra", "binary")
  if (method %in% dist_methods) {
    distance_matrix <- dist(profiles, method = method)
  } else if (method == "distance_correlation") {
    matrix <- data.frame()
    for (i in 1:nrow(profiles)) { # rows
      for (j in 1:nrow(profiles)) { # columns
        if (i == j) {
          matrix[i,i] = 1 # if this cell is NA as.dist does not work probably 
          break
        }
        dist <- dcor(unlist(profiles[i,]), unlist(profiles[j,]))
        # Swich the value so that the profiles with a high correlation are 
        # clustered together
        matrix[i,j] <- 1 - dist 
      }
    }
    profile_names <- rownames(profiles)
    colnames(matrix) <- profile_names[1:length(profile_names) - 1]
    rownames(matrix) <- profile_names
    distance_matrix <- as.dist(matrix)
    
  } else if (method == "mutual_information") {
    empty_rows <- profiles[(rowSums(profiles) == 0),]
    profiles <- profiles[!(rowSums(profiles) == 0),] 
    distance_matrix <- mutualInfo(as.matrix(profiles))
    distance_matrix <- max(distance_matrix, na.rm=TRUE) - distance_matrix
  } else if (method == "pearson") {
    distance_matrix <-  cor.dist(as.matrix(profiles))
  }
  write.table(as.matrix(distance_matrix), file = paste0("../distance_matrix_",
                                                        method),
              col.names = TRUE, row.names = TRUE, quote = FALSE, sep = " \t")
  return(distance_matrix)
}