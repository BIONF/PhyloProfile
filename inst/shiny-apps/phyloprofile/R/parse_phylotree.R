#' Functions for working with phylogenetic trees

#' Check the validity of input newick tree
#' @export
#' @param filein input newick tree
#' @param main_input main phylogenetic profile input
#' @param subset_taxa list of all input taxa
#' @return checking results (1 = missing parenthesis; 2 = missing comma;
#' 3 = tree has singleton; or list of missing taxa in profile)
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

check_newick <- function(filein, main_input, subset_taxa){
  tree <- read.table(file = filein$datapath,
                     header = FALSE,
                     check.names = FALSE,
                     comment.char = "",
                     fill = FALSE)

  # get tree structure
  tree_struc <- gsub(regex("\\w"), "", as.character(tree$V1))

  open <- str_count(tree_struc, "\\(")
  close <- str_count(tree_struc, "\\)")
  comma <- str_count(tree_struc, "\\,")
  singleton <- str_count(tree_struc, "\\(\\)")

  # return check
  if (is.null(main_input)) {
    return(0) # don't check if main input is absent
  } else {
    if (singleton > 0) {
      return(3) # tree contains singleton
    }
    if (open != close) {
      return(1) # missing parenthesis
    } else {
      if ((comma - open) > 1 | (comma - open) < 0) {
        return(2) # missing comma
      } else {
        # get list of tips
        node_string <- gsub(regex("\\W+"), "#", as.character(tree$V1))
        node_list <- unlist(strsplit(node_string, "#"))
        # list of input taxa
        input_taxa <- subset_taxa

        missing_taxa <- list()
        j <- 1
        for (i in 1:length(node_list)) {
          if (nchar(node_list[i]) > 0 & !(node_list[i] %in% input_taxa)) {
            missing_taxa[[j]] <- node_list[i]
            j <- j + 1
          }
        }
        if (length(missing_taxa) > 0) {
          # contains taxa that not exist in main input
          return(paste(missing_taxa, collapse = "; "))
        } else {
          return(0)
        }
      }
    }
  }
  return(0)
}

#' Create rooted tree from a phylogenetic profile matrix
#' @export
#' @param df phylogenetic profile used for creating tree
#' @param root_taxon taxon used for rooting tree
#' @return rooted tree based on root_taxon
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

create_rooted_tree <- function(df, root_taxon){
  # calculate distance matrix
  taxdis <- tryCatch(taxa2dist(df), error = function(e) e)
  # create tree
  tree <- as.phylo(hclust(taxdis))
  # root tree
  tree <- root(tree, outgroup = root_taxon, resolve.root = TRUE)
  # return
  return(tree)
}

#' Sort supertaxa list based on chosen reference taxon
#' @export
#' @param tree taxonomy tree
#' @return list of taxa sorted according to the taxonomy tree
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

sort_taxa_from_tree <- function(tree){
  is_tip <- tree$edge[, 2] <= length(tree$tip.label)
  ordered_tips <- tree$edge[is_tip, 2]
  taxon_list <- rev(tree$tip.label[ordered_tips])
  return(taxon_list)
}

#' Funciton used to replace read.tree function of APE package
#' in case input tree has singletons
#' Read a Newick string with node labels & (possible) singles
#' @export
#' @param
#' @return
#' @author Liam J. Revell 2013

read.newick <- function(file = "", text){
  # check to see if reading from file
  if (file != "") text <- scan(file, sep = "\n", what = "character")
  if (length(text) > 1) {
    tree <- lapply(text, newick)
    class(tree) <- "multiPhylo"
  } else tree <- newick(text)
  return(tree)
}

#' Main Newick string function
#' @export
#' @param
#' @return
#' @author Liam J. Revell 2013

newick <- function(text){
  text <- unlist(strsplit(text, NULL))
  tip.label <- vector(mode = "character")
  node.label <- vector(mode = "character")
  edge <- matrix(c(1, NA), 1, 2)
  edge.length <- vector()
  currnode <- 1
  Nnode <- currnode
  i <- j <- k <- 1
  while (text[i] != ";") {
    if (text[i] == "(") {
      if (j > nrow(edge)) edge <- rbind(edge, c(NA, NA))
      edge[j, 1] <- currnode
      i <- i + 1
      # is the next element a label?
      if (is.na(match(text[i],
                      c("(", ")", ",", ":", ";")))) {
        temp <- get_label(text, i)
        tip.label[k] <- temp$label
        i <- temp$end
        edge[j, 2] <- k
        k <- k + 1
        # is there a branch length?
        if (text[i] == ":") {
          temp <- get_edge_length(text, i)
          edge.length[j] <- temp$edge.length
          i <- temp$end
        }
      } else if (text[i] == "(") {
        Nnode <- Nnode + 1 # creating a new internal node
        currnode <- Nnode
        edge[j, 2] <- currnode # move to new internal node
      }
      j <- j + 1
    } else if (text[i] == ")") {
      i <- i + 1
      # is the next element a label?
      if (is.na(match(text[i], c("(", ")", ",", ":", ";")))) {
        temp <- get_label(text, i)
        node.label[currnode] <- temp$label
        i <- temp$end
      }
      # is there a branch length?
      if (text[i] == ":") {
        temp <- get_edge_length(text, i)
        if (currnode > 1) {
          ii <- match(currnode, edge[, 2])
          edge.length[ii] <- temp$edge.length
        } else root.edge <- temp$edge.length
        i <- temp$end
      }
      # move down the tree
      if (currnode > 1) currnode <- edge[match(currnode, edge[, 2]), 1]
    } else if (text[i] == ",") {
      if (j > nrow(edge)) edge <- rbind(edge, c(NA, NA))
      edge[j, 1] <- currnode
      i <- i + 1
      # is the next element a label?
      if (is.na(match(text[i],
                      c("(", ")", ",", ":", ";")))) {
        temp <- get_label(text, i)
        tip.label[k] <- temp$label
        i <- temp$end
        edge[j, 2] <- k
        k <- k + 1
        # is there a branch length?
        if (text[i] == ":") {
          temp <- get_edge_length(text, i)
          edge.length[j] <- temp$edge.length
          i <- temp$end
        }
      } else if (text[i] == "(") {
        Nnode <- Nnode + 1 # creating a new internal node
        currnode <- Nnode
        edge[j, 2] <- currnode # move to internal node
      }
      j <- j + 1
    }
  }
  Ntip <- k - 1
  edge[edge > 0] <- edge[edge > 0] + Ntip
  edge[edge < 0] <- edge[edge < 0]
  edge.length[is.na(edge.length)] <- 0
  if (length(edge.length) == 0) edge.length <- NULL
  node.label[is.na(node.label)] <- ""
  if (length(node.label) == 0) node.label <- NULL

  # assemble into "phylo" object
  tree <- list(edge = edge,
               Nnode = as.integer(Nnode),
               tip.label = tip.label,
               edge.length = edge.length,
               node.label = node.label)
  class(tree) <- "phylo"
  return(tree)
}

#' Get node label
#' @export
#' @param
#' @return list of node labels
#' @author Liam J. Revell 2013

get_label <- function(text, start, stop.char = c(",", ":", ")")){
  i <- 0
  label <- vector()
  while (is.na(match(text[i + start], stop.char))) {
    label[i + 1] <- text[i + start]
    i <- i + 1
  }
  return(list(label = paste(label, collapse = ""), end = i + start))
}

#' Gets branch length
#' @export
#' @param
#' @return list of branch lengths
#' @author Liam J. Revell 2013

get_edge_length <- function(text, start){
  i <- start + 1; m <- 1
  temp <- vector()
  while (is.na(match(text[i], c(",", ")", ";")))) {
    temp[m] <- text[i]
    i <- i + 1
    m <- m + 1
  }
  return(list(edge.length = as.numeric(paste(temp, collapse = "")),
              end = i))
}
