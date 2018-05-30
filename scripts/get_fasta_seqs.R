

# Fasta output ----------------------------------------------------------------
# input$main_input, input$demo_data
fasta_out_data <- function(data_out, filein, demo_data,
                           input_type_fasta,
                           one_seq_fasta,
                           path,
                           dir_format,
                           file_ext,
                           id_format,
                           long_df){
  
  fasta_out_df <- data.frame()
  
  # check main input
  if (!is.null(filein)){
    input_type <- check_input_vadility(filein)
  } else{
    input_type <- "NA"
  }
  
  # get seqs for AMPK-TOR and microsporidia ONLINE demo data
  if (demo_data == "ampk-tor" | demo_data == "lca-micros"){
    fasta_url <- paste0("https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/demo/fasta_file/concatenatedSeq.fa")
    if (demo_data == "ampk-tor"){
      fasta_url <- paste0("https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/expTestData/ampk-tor/ampk-tor.extended.fa")
    }
    
    if (url.exists(fasta_url)){
      # load fasta file
      fa_file <- as.data.frame(read.table(fasta_url,
                                          sep = "\t",
                                          header = F,
                                          fill = T,
                                          stringsAsFactors = FALSE,
                                          quote = ""))
      fa_df <- data.frame("seqID" = fa_file$V1[grepl(">", fa_file$V1)],
                          "seq" = fa_file$V1[!grepl(">", fa_file$V1)],
                          stringsAsFactors = FALSE)
      
      # get sequences
      for (j in 1:nrow(data_out)){
        seq_id <- as.character(data_out$orthoID[j])
        group_id <- as.character(data_out$geneID[j])
        
        seq <- as.character(fa_df$seq[fa_df$seqID == paste0(">", seq_id)])
        fasta_out <- paste(paste0(">", group_id, "|", seq_id),
                           seq,
                           sep = "\n")
        fasta_out_df <- rbind(fasta_out_df, as.data.frame(fasta_out))
      }
    } else {
      fasta_out <- paste0(fasta_url, " not found!!!")
      fasta_out_df <- rbind(fasta_out_df, as.data.frame(fasta_out))
    }
  }
  
  # get seqs for fasta main input
  if (input_type == "fasta"){
    file <- filein$datapath
    fasta_file <- readAAStringSet(file)
    
    seq_name <- names(fasta_file)
    sequence <- paste(fasta_file)
    # data frame contains all sequences from input file
    fa <- data.frame(seq_name, sequence)
    
    for (j in 1:nrow(data_out)){
      seq_id <- paste0(as.character(data_out$geneID[j]),
                       "|ncbi",
                       as.character(data_out$ncbiID[j]),
                       "|",
                       as.character(data_out$orthoID[j]))
      
      seq <- fa$sequence[pmatch(seq_id, fa$seq_name)]
      
      if (length(seq[1]) < 1){
        fasta_out <- paste0(seq_id,
                            " not found in ",
                            file,
                            "! Please check again!")
      } else{
        fasta_out <- paste(paste0(">", seq_id), seq[1], sep = "\n")
      }
      fasta_out_df <- rbind(fasta_out_df, as.data.frame(fasta_out))
    }
  }
  
  if (input_type == "oma"){
    for(j in 1:nrow(data_out)){
      seq_id <- as.character(data_out$orthoID[j])
      group_id <- as.character(data_out$geneID[j])
      fasta_out <- get_fasta_oma(seq_id, group_id, long_df)
      fasta_out_df <- rbind(fasta_out_df, as.data.frame(fasta_out))
    }
    
    
  } else {
    # get seqs for extended.fa
    if (demo_data == "none" & input_type_fasta == "Concatenated fasta file"){
      if (!is.null(one_seq_fasta)){
        fas_in <- one_seq_fasta
        file <- toString(fas_in$datapath)
        
        # read fasta file and save sequences into dataframe
        fasta_file <- readAAStringSet(file)
        
        seq_name <- names(fasta_file)
        sequence <- paste(fasta_file)
        # data frame contains all sequences from input file
        fa <- data.frame(seq_name, sequence)
        
        # get selected sequences
        for (j in 1:nrow(data_out)){
          seq_id <- as.character(data_out$orthoID[j])
          group_id <- as.character(data_out$geneID[j])
          seq <- fa$sequence[pmatch(seq_id, fa$seq_name)]
          flag <- 1
          if (is.na(seq)){
            seq_id <- paste0(as.character(data_out$geneID[j]),
                             "|ncbi",
                             as.character(data_out$ncbiID[j]),
                             "|",
                             as.character(data_out$orthoID[j]))
            seq <- fa$sequence[pmatch(seq_id, fa$seq_name)]
            flag <- 0
          }
          if (length(seq[1]) < 1){
            fasta_out <- paste0(seq_id,
                                " not found in ",
                                file,
                                "! Please check the header format in FASTA file!")
          } else{
            if (!is.na(seq[1])){
              if (flag == 1){
                fasta_out <- paste(paste0(">", group_id, "|", seq_id),
                                   seq[1],
                                   sep = "\n")
              } else {
                fasta_out <- paste(paste0(">", seq_id), seq[1], sep = "\n")
              }
              
            } else {
              fasta_out <- {
                paste0(seq_id,
                       " not found in uploaded FASTA file!!! Please check again!!!")
              }
            }
          }
          fasta_out_df <- rbind(fasta_out_df, as.data.frame(fasta_out))
        }
      } else {
        if (input_type != "fasta"){
          fasta_out <- {
            paste0("Please provide FASTA file(s) in Input & settings page!")
          }
          fasta_out_df <- rbind(fasta_out_df, as.data.frame(fasta_out))
        }
      }
    }
    
    # get seqs for other cases (input offline fasta files in a folder)
    if (demo_data == "none"
        & input_type_fasta == "Fasta folder"
        & input_type != "fasta"){
      if (path != ""){
        
        # get list of species IDs
        if (id_format == 1){
          spec_df <- {
            as.data.frame(str_split_fixed(str_reverse(as.character(data_out$orthoID)),
                                          ":", 2))
          }
          spec_df$spec_id <- str_reverse(as.character(spec_df$V2))
        } else if (id_format == 2){
          spec_df <- {
            as.data.frame(str_split_fixed(str_reverse(as.character(data_out$orthoID)),
                                          "@", 2))
          }
          
          spec_df$spec_id <- str_reverse(as.character(spec_df$V2))
        } else if (id_format == 3){
          spec_df <- {
            as.data.frame(str_split_fixed(str_reverse(as.character(data_out$orthoID)),
                                          "|", 2))
          }
          spec_df$spec_id <- str_reverse(as.character(spec_df$V2))
        }
        
        # read all specices FASTA files at once
        fa <- data.frame()
        for (i in 1:length(levels(as.factor(spec_df$specID)))){
          spec_id <- as.character(levels(as.factor(spec_df$specID))[i])
          
          # full path fasta file
          file <- paste0(path, "/", spec_id, ".", file_ext)
          if (dir_format == 2){
            file <- paste0(path, "/", spec_id, "/", spec_id, ".", file_ext)
          }
          
          # read fasta file and save sequences into dataframe
          if (file.exists(file)){
            fasta_file <- readAAStringSet(file)
            
            seq_name <- names(fasta_file)
            sequence <- paste(fasta_file)
            # data frame contains all sequences from input file
            fa <- rbind(fa, data.frame(seq_name, sequence))
          }
        }
        # now get selected sequences
        if (nrow(fa) > 0){
          for (j in 1:nrow(data_out)){
            seq_id <- as.character(data_out$orthoID[j])
            group_id <- as.character(data_out$geneID[j])
            
            seq <- fa$sequence[pmatch(seq_id, fa$seq_name)]
            
            if (length(seq[1]) < 1){
              fasta_out <- {
                paste0(seq_id,
                       " not found in ",
                       file,
                       "! Please check id_format in FASTA config again!")
              }
            } else{
              fasta_out <- paste(paste0(">", group_id, "|", seq_id),
                                 seq[1],
                                 sep = "\n")
            }
            fasta_out_df <- rbind(fasta_out_df, as.data.frame(fasta_out))
          }
        } else {
          fasta_out <- paste0("No fasta file has been found in ",
                              path,
                              "!!! Please check the full path to FASTA folder and the id_format (header format) in FASTA config again!!!")
          fasta_out_df <- rbind(fasta_out_df,
                                as.data.frame(fasta_out))
        }
      } else {
        fasta_out <- {
          paste0("Please provide FASTA files in Input & settings page!")
        }
        fasta_out_df <- rbind(fasta_out_df,
                              as.data.frame(fasta_out))
      }
    }
  }
  # remove duplicated sequences
  fasta_out_df <- fasta_out_df[!duplicated(fasta_out_df), ]

  return(fasta_out_df)
}