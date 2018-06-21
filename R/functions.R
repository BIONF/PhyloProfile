# unsorting function to keep user defined geneID order ------------------------
unsort_id <- function(data, order){
  data$geneID <- as.factor(data$geneID)
  if (order == FALSE) {
    # keep user defined geneID order
    data$geneID <- factor(data$geneID, levels = unique(data$geneID))
  }
  return(data)
}

# Essential functions =========================================================

# get list with all taxanomy ranks --------------------------------------------
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

# Get name list ("data/taxonNamesReduced.txt") --------------------------------
get_name_list <- function (as_character, delete_duplicated) {
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

# Get taxa list ("data/taxonomyMatrix.txt") -----------------------------------
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

# reverse string --------------------------------------------------------------
str_reverse <- function(x){
  sapply(lapply(strsplit(x, NULL), rev), paste, collapse = "")
}

# get last n characters from string x  ----------------------------------------
substr_right <- function(x, n){
  substr(x, nchar(x) - n + 1, nchar(x))
}

# check internet connection  --------------------------------------------------
has_internet <- function(){
  !is.null(curl::nslookup("r-project.org", error = FALSE))
}

# FUNCTION FOR CLUSTERING PROFILES --------------------------------------------
clustered_gene_list <- function(data, dist_method, cluster_method){

  # do clustering
  row.order <- hclust(dist(data, method = dist_method),
                      method = cluster_method)$order
  col.order <- hclust(dist(t(data), method = dist_method),
                      method = cluster_method)$order

  # re-order data accoring to clustering
  dat_new <- data[row.order, col.order]

  # return clustered gene ID list
  clustered_gene_ids <- as.factor(row.names(dat_new))
  clustered_gene_ids
}

# calculate percentage of present species -------------------------------------
calc_pres_spec <- function(taxa_md_data, taxa_count){
   # taxa_md_data = df("geneID",
   #                 "ncbiID",
   #                 "orthoID",
   #                 "var1",
   #                 "var2",
   #                 "paralog",
   #                 ....,
   #                 "supertaxon")
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

# Calculate the contigency table for the fisher exact test --------------------
get_table <- function(profile_1, profile_2){
  contigency_table <- data.frame(c(0,0), c(0,0))
  for(i in 1:length(profile_1)){
    if(profile_1[i] == 1){
      if(profile_2[i] == 1) {
        contigency_table[1,1] <- contigency_table[1,1] + 1
      } else {
        contigency_table[2,1] <- contigency_table[2,1] + 1
      }
    } else{
      if(profile_2[i] == 1) {
        contigency_table[1,2] <- contigency_table[1,2] + 1
      } else {
        contigency_table[2,2] <- contigency_table[2,2] + 1
      }
    }
  }
  contigency_table
}

# calculate percentage of present species -------------------------------------
calc_pres_spec <- function(taxa_md_data, taxa_count){
  # taxa_md_data = df("geneID",
  #                 "ncbiID",
  #                 "orthoID",
  #                 "var1",
  #                 "var2",
  #                 "paralog",
  #                 ....,
  #                 "supertaxon")
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


# GROUP COMPARISON ============================================================

# print list of available taxa ------------------------------------------------
taxa_select_gc <- function(rank_select_gc, subset_taxa){
  
  # if there is no rank set, there can not be any available taxa
  if (length(rank_select_gc) == 0) return()
  else{
    
    # load list of unsorted taxa
    dt <- get_taxa_list(TRUE, subset_taxa)
    
    # load list of taxon name
    name_list <- get_name_list(TRUE, FALSE)
    
    # get rank name from rank_select
    rank_name <- substr(rank_select_gc, 4, nchar(rank_select_gc))
    choice <- as.data.frame
    choice <- rbind(dt[rank_name])
    colnames(choice) <- "ncbiID"
    choice <- merge(choice, name_list, by = "ncbiID", all = FALSE)
    return(choice)
  }
}
