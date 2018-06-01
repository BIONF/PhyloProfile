library(devtools)
if (!require("roma")) install_github("trvinh/roma")
library ("roma")

# OMA IDs or Uniprot IDs as Input =============================================
# check if an object is a OMA or Uniport id -----------------------------------
check_oma_id <- function(id){
  id <- as.character(id)
  data_id <- getData("protein", id)
  if (!is.null(data_id$entry_nr)) return (TRUE)
  else {
    return (FALSE)
  }
}

# get the members for a OMA or Uniprot id -------------------------------------
get_members <- function(id, output_type){
  if (output_type == "HOG"){
    # get the members of the Hierarchical Orthologous Group
    members <- getAttribute(getHOG(id = id, level = "root", members = TRUE),
                            "members")
  } else if (output_type == "OG"){
    # get the members of the Ortholoug group
    members <- getAttribute(getData(type = "group", id = id),
                            "members")
  } else if (output_type == "PAIR"){
    # get the members of the Orthologous Pair
    data <- getData(type = "protein", id = id)
    members <- resolveURL(data$orthologs)
  }
  return(members)
}

# ncbiID =  "ncbi" + taxonomy ID ----------------------------------------------
get_ncbi_id <- function(oma_id){
  species_name <- substr(oma_id, 1, 5)
  taxonomy_id <- getTaxonomy(species_name)$id
  ncbi_id <- paste0("ncbi", taxonomy_id)
}

# Transform the input of ids to a dataframe in long format --------------------
oma_ids_to_long <- function(oma_ids, output_type){
  # "geneID ncbiID orthoID"
  long_dataframe <- data.frame()
  row_nr <- 1

  for(id in oma_ids){
    if(check_oma_id(id)){############## NEU (Makeing sure that it is a valid id, the function is the first in the script)
      start_id = Sys.time() ############## NEU (For the time stamp)
      members <- get_members(id, output_type)
      
      gene_id <- paste0("OG_", id)
      oma_id <- getAttribute(getData("protein", id), "omaid")############## NEU
      
      ncbi <- get_ncbi_id(oma_id)############## NEU (need an oma id not an uniprot id)
      long_dataframe[row_nr,1] <- gene_id
      long_dataframe[row_nr,2] <- ncbi
      long_dataframe[row_nr,3] <- oma_id ############## NEU (so it is consistent)
      row_nr <- row_nr +1

      # Get orthoID and ncbiID for each member of the hogs
      for (i in 1:nrow(members)){
        member <- members[i, ]
        ortho_id <- member$omaid # use the oma ID as ortho ID
        
        # Data for the current member (ortho ID)
        member_data <- getData("protein", ortho_id)
        
        ncbi_id <- get_ncbi_id(member_data$omaid)
        

        # New line for the long format
        long_dataframe[row_nr,1] <- gene_id
        long_dataframe[row_nr,2] <- ncbi_id
        long_dataframe[row_nr,3] <-ortho_id
        row_nr <- row_nr + 1
      }
      end_id <- Sys.time()############## NEU (for the time stamp)
      time <- end_id - start_id############## NEU (for the time stamp)
      #print(paste("runtime for", id, "with", nrow(members),"members:", time, sep = " "  ))############## NEU (time stamp for each input id)
    } else {############## NEU
      print(paste(id, "is not a valid oma or uniprot id"))############## NEU 
    } ############## NEU
  }
  colnames(long_dataframe) <- c("geneID", "ncbiID", "orthoID")
  return(long_dataframe)
}

# Get the fasta file for a list of ids ----------------------------------------
oma_ids_to_fasta <- function(oma_ids, output_type){
  lines_fasta <- c()
  for (id in oma_ids){
    members <- get_members(id, output_type)
    gene_id <- paste0("OG_", id)

    for (i in 1:nrow(members)){
      member <- members[i, ]
      ortho_id <- member$omaid # use the oma ID as ortho ID

      # Data for the current member (ortho ID)
      member_data <- getData("protein", ortho_id)
      ncbi_id <- get_ncbi_id(member_data$omaid)

      # New lines for the fasta format
      header_sequence <- paste0(">", gene_id, "|", ncbi_id, "|", ortho_id)
      sequence <- as.character(member_data$sequence)
      lines_fasta <- append(lines_fasta, header_sequence)
      lines_fasta <- append(lines_fasta, sequence)
    }

  }
  return(lines_fasta)
  # fasta_file <- file("output.txt")
  # writeLines(lines_fasta, fasta_file)
  # close(fasta_file)
}

long_to_fasta <- function(long){
  lines_fasta <- c()
  for(row in 1:nrow(long)){
    gene_id <- long$geneID[row]
    ncbi_id <- long$ncbiID[row]
    ortho_id <- long$orthoID[row]

    header_sequence <- paste0(">", gene_id, "|", ncbi_id, "|", ortho_id)
    sequence <- as.character(getData("protein", ortho_id)$sequence)
    lines_fasta <- append(lines_fasta, header_sequence)
    lines_fasta <- append(lines_fasta, sequence)
  }
  return(lines_fasta)
}

get_fasta_oma <- function(seq_id, group_id, long_df){
  selected_df <- subset(long_df,
                        long_df$geneID == group_id & long_df$orthoID == seq_id)
  header <- paste0(">",
                   selected_df$geneID,
                   "|",
                   selected_df$ncbiID,
                   "|",
                   selected_df$orthoID)
  seq <- as.character(getData("protein",  selected_df$orthoID)$sequence)
  return(paste(header, seq, sep = "\n"))
}

# Get the domain file for a list of ids ---------------------------------------
oma_ids_to_domain <- function(oma_ids, output_type){
  domain_data <- data.frame()
  row_nr <- 0

  for (id in oma_ids){
    members <- get_members(id, output_type)
    gene_id <- paste0("OG_", id)

    for (i in 1:nrow(members)){
      member <- members[i, ]
      ortho_id <- member$omaid # use the oma ID as ortho ID

      # seedID = geneID#orthoID
      seed_id <- paste0(gene_id, "#", ortho_id)

      # Data for the current member (ortho ID)
      member_data <- getData("protein", ortho_id)

      # length of the sequence
      length <- member_data$sequence_length

      # Informations about the domain
      domains <- resolveURL(member_data$domains)
      regions <- domains$regions
      regions$feature <- paste(regions$source, regions$name, sep = " ")

      for (i in 1:nrow(regions)){
        row_nr <- row_nr + 1
        domain <- regions[i, ]

        location <- unlist(strsplit(domain$location, ":"))

        domain_data[row_nr,1] <- seed_id
        domain_data[row_nr,2] <- ortho_id
        domain_data[row_nr,3] <- length
        domain_data[row_nr,4] <- domain$feature
        domain_data[row_nr,5] <- location[1]
        domain_data[row_nr,6] <- location[2]
      }
    }
  }
  colnames(domain_data) <- c("seedID",
                             "orthoID",
                             "length",
                             "feature",
                             "start",
                             "end")
  return(domain_data)
}

long_to_domain <- function(long){
  domain_data <- data.frame()
  row_nr <- 0
  for(row in 1:nrow(long)){
    gene_id <- long$geneID[row]
    ortho_id <- long$orthoID[row]

    # seedID = geneID#orthoID
    seed_id <- paste0(gene_id, "#", ortho_id)

    # Data for the current member (ortho ID)
    member_data <- getData("protein", ortho_id)

    # length of the sequence
    length <- member_data$sequence_length

    # Informations about the domain
    domains <- resolveURL(member_data$domains)
    regions <- domains$regions
    regions$feature <- paste(regions$source, regions$name, sep = " ")
    if(!is.null(nrow(regions))){############## NEU (Do catch when there are no regions found)
      for (i in 1:nrow(regions)){
        row_nr <- row_nr + 1
        domain <- regions[i, ]

        location <- unlist(strsplit(domain$location, ":")) 

        domain_data[row_nr,1] <- seed_id
        domain_data[row_nr,2] <- ortho_id
        domain_data[row_nr,3] <- length
        domain_data[row_nr,4] <- domain$feature
        domain_data[row_nr,5] <- location[1]
        domain_data[row_nr,6] <- location[2]
      }
    }############## NEU
  }
  colnames(domain_data) <- c("seedID",
                             "orthoID",
                             "length",
                             "feature",
                             "start",
                             "end")
  
  domain_data$start <- as.integer(domain_data$start) ############## NEU
  domain_data$end <- as.integer(domain_data$end) ############## NEU
  return(domain_data)
}
