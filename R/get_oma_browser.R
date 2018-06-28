#' Get OMA information functions
#' using OmaDB package (https://github.com/klarakaleb/OmaDB)
library(OmaDB)

#' check OMA IDs or Uniprot IDs as valid Input --------------------------------
#' @export
#' @param id_list list of ids needs to be checked
#' @return list of invalid IDs (not readable for OMA)
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

check_oma_id <- function(ids) {
  invalid <- list()
  for (id in ids) {
    id <- as.character(id)
    data <- OmaDB::getData("protein", id)
    if (is.null(data$entry_nr)) invalid <- c(invalid, id)
  }
  return(invalid)
}

#' get the members for a OMA or Uniprot id ------------------------------------
#' @export
#' @param id oma ID
#' @param ortho_type type of OMA orthologs
#' @return list of ortholog members
#' @author Carla Mölbert {carla.moelbert@gmx.de}

get_members <- function(id, ortho_type) {
  # get the members of the Hierarchical Orthologous Group
  if (ortho_type == "HOG") {
    members <- OmaDB::getHOG(id = id,
                             level = "root",
                             members = TRUE)$members$omaid
  } 
  # get the members of the Ortholoug group
  else if (ortho_type == "OG") {
    members <- OmaDB::getData(type = "group",
                              id = id)$members$omaid
  } 
  # get the members of the Orthologous Pair
  else if (ortho_type == "PAIR") {
    members <- OmaDB::resolveURL(OmaDB::getData(type = "protein",
                                                id = id)$orthologs)$omaid
    # add query ID into output list
    seed <- OmaDB::getData("protein",id)$omaid
    members <- c(seed,members)
  }
  
  return(members)
}

#' get domain annotation for a domain URL of an OMA sequence ------------------
#' @export
#' @param domainURL URL address for domain annotation of ONE OMA id
#' @return data frame contains feature names, start and end positions
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

get_domain_from_url <- function(domainURL) {
  domains <- OmaDB::resolveURL(domainURL)$regions
  domains$feature <- NA
  domains$start <- NA
  domains$end <- NA
  for (i in 1:nrow(domains)) {
    pos <- unlist(strsplit(domains$location[i], ":"))
    domains[i,]$start <- pos[1]
    domains[i,]$end <- pos[2]
    
    if (nchar(domains$name[i]) > 0) {
      domains[i,]$feature <- paste0(domains$source[i],"_",domains$domainid[i]," (",domains$name[i],")")
    } else {
      domains[i,]$feature <- paste0(domains$source[i],"_",domains$domainid[i])
    }
    domains[i,]$feature <- gsub("#", "-", domains[i,]$feature)
  }
  return(domains[, c("feature","start","end")])
}

#' get taxonomy ID, sequence and annotation for an OMA sequence ---------------
#' @export
#' @param id oma ID
#' @return data frame contains above info for ONE OMA id
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

get_data_for_one_oma <- function(id) {
  omaDf <- data.frame(
    "ortho_id" = character(),
    "taxon_id" = character(),
    "seq" = character(),
    "length" = numeric(),
    "domains" = character(),
    stringsAsFactors = FALSE
  )
  
  # get ncbi taxonomy id
  spec_name <- substr(id, 1, 5)
  taxon_id <- paste0("ncbi",
                     OmaDB::getTaxonomy(members = spec_name,
                                        newick = FALSE)$id)
  # get raw data
  raw <- OmaDB::getData("protein",id)
  
  # get sequence
  seq <- as.character(raw$sequence)
  
  # get sequence length
  length <- raw$sequence_length
  
  # get annotation
  domainDf <- get_domain_from_url(raw$domains)
  domainDf_join <- c(domainDf, sep = "#")
  domains <- paste(unlist(do.call(paste, domainDf_join)), collapse = "\t")
  
  # return data frame contains all info
  omaDf[1,] <- c(id, taxon_id, seq, length, domains)
  return(omaDf)
}

#' get all oma info for list of IDs and save into a data frame (final_oma_df) -
#' @export
#' @param id_list list of input OMA/UniProt ids
#' @param ortho_type type of OMA ortholog
#' @return data frame contains all OMA information (as in get_data_for_one_oma)
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

# get_oma_browser <- function(id_list, ortho_type) {
#   final_omaDf <- data.frame()
#   
#   for (seed_id in id_list) {
#     # get members
#     members <- get_members(seed_id, ortho_type)
#     oma_seed_id <- OmaDB::getData("protein",seed_id)$omaid
#     # get all data
#     for (ortho in members) {
#       orthoDf <- get_data_for_one_oma(ortho)
#       orthoDf$seed <- seed_id
#       if (ortho == oma_seed_id) {
#         orthoDf$ortho_id <- seed_id
#       }
#       final_omaDf <- rbind(final_omaDf, orthoDf)
#     }
#   }
#   
#   return(final_omaDf)
# }

#' create phylogenetic profile from full OMA dataframe
#' @export
#' @param final_oma_df obtained OMA data for a list of input IDs
#' @return profile in long format
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

create_profile_from_oma <- function(final_oma_df) {
  profile_df <- final_oma_df[, c("seed", "taxon_id", "ortho_id")]
  colnames(profile_df) <- c("geneID", "ncbiID", "orthoID")
  return(profile_df[!duplicated(profile_df), ])
}

#' create domain annotation dataframe for ONE protein
#' @export
#' @param domain_id 
#' @param ortho_id
#' @param length
#' @param domain_list list of all domains and their positions for this protein
#' @return domain annotation in a dataframe
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

create_domain_df <- function(domain_id, ortho_id, length, domain_list) {
  domain_df = data.frame(
    seedID = character(),
    orthoID = character(),
    length = numeric(),
    feature = character(),
    start = numeric(),
    end = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:length(domain_list)) {
    anno_info <- strsplit(domain_list[i], "#")[[1]]
    domain_df[i,] <- c(domain_id,
                       ortho_id,
                       length,
                       anno_info)
  }
  
  return(domain_df)
}

#' create domain annotation dataframe for complete OMA data
#' @export
#' @param final_oma_df obtained OMA data for a list of input IDs
#' @return domain annotation in a dataframe to input into phyloprofile
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

get_all_domains_oma <- function(final_oma_df) {
  oma_domain_df <- data.frame()
  for (i in 1:nrow(final_oma_df)) {
    domainID <- paste0(final_oma_df[i,]$seed, "#", final_oma_df[i,]$ortho_id)
    
    seed_line <- final_oma_df[final_oma_df$ortho_id == final_oma_df[i,]$seed, ]
    seed_domains <- strsplit(as.character(seed_line$domains), "\t")[[1]]
    seed_domainDf <- create_domain_df(domainID,
                                      seed_line$ortho_id,
                                      seed_line$length,
                                      seed_domains)
    oma_domain_df <- rbind(oma_domain_df, seed_domainDf)
      
    ortho_domains <- strsplit(as.character(final_oma_df[i,]$domains), "\t")[[1]]
    ortho_domainDf <- create_domain_df(domainID,
                                       final_oma_df[i,]$ortho_id,
                                       final_oma_df[i,]$length,
                                       ortho_domains)
    oma_domain_df <- rbind(oma_domain_df, ortho_domainDf)
  }
  
  oma_domain_df$length <- as.numeric(oma_domain_df$length)
  oma_domain_df$start <- as.numeric(oma_domain_df$start)
  oma_domain_df$end <- as.numeric(oma_domain_df$end)
  return(oma_domain_df)
}

#' get all fasta sequences for complete OMA data
#' @export
#' @param final_oma_df obtained OMA data for a list of input IDs
#' @return dataframe contains all fasta sequences
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

get_all_fasta_oma <- function(final_oma_df) {
  oma_fasta_df <- data.frame()
  
  fasta_df <- final_oma_df[, c("ortho_id", "seq")]
  for (i in 1:nrow(fasta_df)) {
    seq_id <- as.character(fasta_df$ortho_id[i])
    seq <- as.character(fasta_df$seq[i])
    fasta_out <- paste(paste0(">", seq_id),
                       seq,
                       sep = "\n")
    oma_fasta_df <- rbind(oma_fasta_df, as.data.frame(fasta_out))
  }
  
  return(oma_fasta_df[!duplicated(oma_fasta_df), ])
}

#' get selected fasta sequences from the complete OMA data
#' @export
#' @param final_oma_df obtained OMA data for a list of input IDs
#' @param seq_id sequence need to be returned
#' @return required sequence in fasta format
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

get_selected_fasta_oma <- function(final_oma_df, seq_id) {
  selected_df <- subset(final_oma_df[, c("ortho_id", "seq")],
                        final_oma_df$ortho_id == seq_id)
  header <- paste0(">", selected_df$ortho_id[1])
  return(paste(header, selected_df$seq[1], sep = "\n"))
}