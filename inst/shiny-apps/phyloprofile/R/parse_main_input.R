#' Parse main input functions
#'
#' These functions check the validity of the main input file
#' and convert different input format into long format.

#' check validity of main input file
#' @export
#' @param filein input file
#' @return input file format or type of error
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

source("R/get_oma_browser.R")

check_input_vadility <- function(filein) {
  input_dt <- as.data.frame(read.table(file = filein$datapath,
                                       sep = "\t",
                                       header = FALSE,
                                       check.names = FALSE,
                                       comment.char = "",
                                       fill = TRUE))

  if (is.na(input_dt[1, ncol(input_dt)])) {
    return("moreCol")
  } else {
    names(input_dt) <- as.character(unlist(input_dt[1, ]))

    # XML format (starts with <?xml)
    if (grepl("<?xml", colnames(input_dt)[1])) {
      return("xml")
    }
    # FASTA format (starts with ">" )
    else if (grepl(">", colnames(input_dt)[1]) == TRUE) {
      return("fasta")
    }
    # LONG or WIDE format (starts with "geneID")
    else {
      if (grepl("geneID", colnames(input_dt)[1])) {
        # LONG format
        if (is.na(pmatch("ncbi", colnames(input_dt)[3])) ||
            is.na(pmatch("ncbi", colnames(input_dt)[4])) ||
            is.na(pmatch("ncbi", colnames(input_dt)[5]))) {
          return("long")
        }
        # WIDE format
        else {
          tmp <- input_dt[input_dt == ""][1]
          if (!is.na(tmp) & tmp == "") {
            return("emptyCell")
          } else {
            return("wide")
          }
        }
      }
      # OMA ids
      else {
				invalid_oma <- check_oma_id(levels(input_dt[,1]))
        if (length(invalid_oma) == 0) {
          return("oma")
        } else {
          return(invalid_oma)
        }
      }
    }
  }
}

#' parse orthoXML input file
#' @export
#' @param input_file
#' @return a data frame that contains input data in long-format
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

xml_parser <- function(input_file){
  cmd <- paste("python ",
               getwd(),
               "/scripts/orthoxmlParser.py",
               " -i ",
               input_file,
               sep = "")

  df_in <- as.data.frame(read.table(text = system(cmd, intern = TRUE)))

  # the first row will be the header
  colnames(df_in) <- as.character(unlist(df_in[1, ]))

  df_in <- subset(df_in[df_in$geneID != "geneID", ])
  df_in <- droplevels(df_in)

  return(df_in)
}

#' parse input file in fasta format
#' @export
#' @param input_file
#' @return a data frame that contains input data in long-format
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

fasta_parser <- function(input_file){
  cmd <- paste("python ",
               getwd(),
               "/scripts/fastaParser.py",
               " -i ",
               input_file,
               sep = "")

  df_in <- as.data.frame(read.table(text = system(cmd, intern = TRUE)))

  # the first row will be the header
  colnames(df_in) <- as.character(unlist(df_in[1, ]))
  df_in <- subset(df_in[df_in$geneID != "geneID", ])
  df_in <- droplevels(df_in)

  # remove var1 and var2 columns if they are all NAs
  if (all(is.na(df_in$var2))) {
    df_in <- subset(df_in, select = -c(var2) )
  }
  if (all(is.na(df_in$var1))) {
    df_in <- subset(df_in, select = -c(var1) )
  }

  return(df_in)
}

#' transform wide format input into long format
#' @export
#' @param input_file
#' @return a data frame that contains input data in long-format
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

wide_to_long <- function(input_file){
  wide_dataframe <- as.data.frame(read.table(file = input_file,
                                             sep = "\t",
                                             header = TRUE,
                                             check.names = FALSE,
                                             comment.char = ""))
  long_dataframe <- data.frame()
  row_nr_long <- 0
  ncbi_ids <- colnames(wide_dataframe)

  for (row_nr in 1:nrow(wide_dataframe)) {
    geneID <- wide_dataframe[row_nr, 1]
    for (column_nr in 2:ncol(wide_dataframe)) {
      current_cell <- as.character(wide_dataframe[row_nr, column_nr])
      new_row_info <- unlist(strsplit(current_cell, "#"))
      row_nr_long <- row_nr_long + 1
      long_dataframe[row_nr_long, 1] <- geneID
      long_dataframe[row_nr_long, 2] <- ncbi_ids[column_nr]
      long_dataframe[row_nr_long, 3] <- new_row_info[1]
      long_dataframe[row_nr_long, 4] <- new_row_info[2]
      long_dataframe[row_nr_long, 5] <- new_row_info[3]
    }
  }

  colnames(long_dataframe) <- c("geneID", "ncbiID", "orthoID", "var1", "var2")
  return(long_dataframe)
}

#' create final dataframe for all kinds of input file
#' @export
#' @param input_file
#' @return a data frame that contains input data in long-format
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

create_long_matrix <- function(input_file){
  if (input_file[1] == "lca-micros") {
    long_dataframe <-
      read.table(
        "https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/demo/test.main.long",
        sep = "\t",
        header = TRUE,
        fill = TRUE,
        stringsAsFactors = FALSE
      )
  } else if (input_file[1] == "ampk-tor") {
    long_dataframe <-
      read.table(
        "https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/expTestData/ampk-tor/ampk-tor.phyloprofile",
        sep = "\t",
        header = TRUE,
        fill = TRUE,
        stringsAsFactors = FALSE
      )
  } else {
    filein <- input_file
    if (is.null(filein)) return()
    input_type <- check_input_vadility(filein)

    # XML
    if (input_type == "xml") {
      long_dataframe <- xml_parser(filein$datapath)
    }
    # FASTA
    else if (input_type == "fasta") {
      long_dataframe <- fasta_parser(filein$datapath)
    }
    # LONG
    else if (input_type == "long") {
      long_dataframe <- as.data.frame(read.table(file = filein$datapath,
                                                 sep = "\t",
                                                 header = TRUE,
                                                 check.names = FALSE,
                                                 comment.char = ""))
    }
    # WIDE
    else if (input_type == "wide") {
      long_dataframe <- wide_to_long(filein$datapath)
    }
    # # OMA
    # else if (input_type == "oma") {
    #   # do not load the data before the button is loaded
    #   if (is.null(input$get_data_oma)) return()
    #
    #   isolate({
    #     # do not generate data befor the button was clicked
    #     if (input$get_data_oma[1] == 0) return()
    #     oma_type <- input$selected_oma_type
    #     if (is.null(oma_type)) return()
    #     oma_ids <- as.data.frame(read.table(file = filein$datapath,
    #                                         sep = "\t",
    #                                         header = FALSE,
    #                                         check.names = FALSE,
    #                                         comment.char = ""))
    #
    #     oma_ids[, 1] <- as.character(oma_ids[, 1])
    #     long_dataframe <- oma_ids_to_long(oma_ids[, 1], oma_type)
    #   })
    # }
    else {
      return(NULL)
    }
  }

  # make sure all columns have the same type (factor)
  for (i in 1:ncol(long_dataframe)) {
    long_dataframe[, i] <- as.factor(long_dataframe[, i])
  }

  # long_dataframe$orthoID <- gsub("\\|",":",long_dataframe$orthoID) 
  return(long_dataframe)
}
