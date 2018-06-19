source("R/parse_main_input.R")
source("R/get_oma_browser.R")

#' Parse domain input file
#' and return output as a data frame
#'
#' @main_input
#' @input_type
#' @demo_data
#' @anno_location
#' @file_domain_input
#' @domainPath
#' @session
#' @datapath

parse_domain_input <- function(main_input,
                               input_type,
                               demo_data,
                               anno_location,
                               file_domain_input,
                               domainPath,
                               session,
                               datapath){
  domains <- data.frame()
  files <- c()

  if(input_type == "oma"){
    domains <- long_to_domain(main_input)
  }

  # get domain file(s) from folder
  if(demo_data == "lca-micros" | anno_location == "from folder" ){
    genes <- unlist(main_input$geneID)
    genes <- unique(genes)

    for (gene in genes){
      file <- get_domain_file(gene,
                              demo_data,
                              anno_location,
                              file_domain_input,
                              domainPath,
                              session,
                              datapath)
      files <- append(files, file)
    }
  }
  # or from singel file
  else{
    file <- get_domain_file(NULL,
                            demo_data,
                            anno_location,
                            file_domain_input,
                            domainPath,
                            session,
                            datapath)
    files <- c(file)
  }

  # parse domain file
  if (demo_data == "lca-micros" | demo_data == "ampk-tor"){
    for(file in files){
      domain_df <- as.data.frame(read.csv(file,
                                          sep = "\t",
                                          header = F,
                                          comment.char = "",
                                          stringsAsFactors = FALSE,
                                          quote = ""))
      domains <- rbind(domains, domain_df)
    }

    if (ncol(domains) == 5){
      colnames(domains) <- c("seedID",
                             "orthoID",
                             "feature",
                             "start",
                             "end")
    } else if (ncol(domains) == 6){
      colnames(domains) <- c("seedID",
                             "orthoID",
                             "feature",
                             "start",
                             "end",
                             "weight")
    } else if (ncol(domains) == 7){
      colnames(domains) <- c("seedID",
                             "orthoID",
                             "feature",
                             "start",
                             "end",
                             "weight",
                             "path")
    }
    domains$length <- max(domains$end)

  } else {
    for(file in files){
      if (file != FALSE){
        exeptions <- c("noFileInput", "noSelectHit",
                       "noSelectHit", "noFileInFolder")
        if(!(file %in% exeptions)){
          domain_df <- as.data.frame(read.table(file,
                                                sep = "\t",
                                                header = FALSE,
                                                comment.char = ""))
          domains <- rbind(domains, domain_df)
        }
      }
    }

    if (ncol(domains) == 5){
      colnames(domains) <- c("seedID",
                             "orthoID",
                             "feature",
                             "start",
                             "end")
    } else if (ncol(domains) == 6){
      colnames(domains) <- c("seedID",
                             "orthoID",
                             "length",
                             "feature",
                             "start",
                             "end")
    } else if (ncol(domains) == 7){
      colnames(domains) <- c("seedID",
                             "orthoID",
                             "length",
                             "feature",
                             "start",
                             "end",
                             "weight")
    } else if (ncol(domains) == 8){
      colnames(domains) <- c("seedID",
                             "orthoID",
                             "length",
                             "feature",
                             "start",
                             "end",
                             "weight",
                             "path")
    }
  }

  print(head(domains))

  domains$seedID <- gsub("\\|",":",domains$seedID)
  domains$orthoID <- gsub("\\|",":",domains$orthoID)

  return(domains)
}

#' Get domain file
#'
#' @seed
#' @demo_data
#' @anno_location
#' @file_domain
#' @domain_path
#' @session
#' @datapath

get_domain_file <- function(seed,
                            demo_data,
                            anno_location,
                            file_domain,
                            domain_path,
                            session,
                            datapath){

  # for demo data
  if (demo_data == "lca-micros" | demo_data == "ampk-tor"){
    updateButton(session, "do_domain_plot", disabled = FALSE)
    if (demo_data == "lca-micros"){
      file_domain <- {
        suppressWarnings(paste0("https://github.com/BIONF/phyloprofile-data/blob/master/demo/domain_files/",
                                seed,
                                ".domains?raw=true"))
      }
    } else {
      file_domain <- {
        suppressWarnings(paste0("https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/expTestData/ampk-tor/ampk-tor.domains_F"))
      }
    }
  }
  # for user input file
  else {
    if (anno_location == "from file"){
      # file_domain <- fileDomain_input
      if (is.null(file_domain)){
        file_domain <- "noFileInput"
      } else {
        if (is.null(info)){
          file_domain <- "noSelectHit"
          updateButton(session, "do_domain_plot", disabled = TRUE)
        } else {
          updateButton(session, "do_domain_plot", disabled = FALSE)
          file_domain <- file_domain$datapath
        }
      }
    } else {
      if (is.null(info)){
        file_domain <- "noSelectHit"
        updateButton(session, "do_domain_plot", disabled = TRUE)
      } else {
        ### check file extension
        all_extension <- c("txt", "csv", "list", "domains", "architecture")
        flag <- 0
        for (i in 1:length(all_extension)){

          file_domain <- paste0(domain_path,
                                "/",
                                seed,
                                ".",
                                all_extension[i])
          if (file.exists(file_domain) == TRUE){
            updateButton(session,
                         "do_domain_plot",
                         disabled = FALSE)
            flag <- 1
            break ()
          }
        }

        if (flag == 0){
          file_domain <- "noFileInFolder"
          updateButton(session,
                       "do_domain_plot",
                       disabled = TRUE)
        }
      }
    }
  }
  return (file_domain)
}
