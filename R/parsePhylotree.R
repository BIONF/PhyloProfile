# Functions for working with phylogenetic trees

#' Check the validity of input newick tree
#' @export
#' @param tree input newick tree
#' @param inputTaxonID list of all input taxon IDs for the phylogenetic profiles
#' @return Possible formatting error of input tree. 0 = suitable tree for using
#' with PhyloProfile, 1 = missing parenthesis; 2 = missing comma; 
#' 3 = tree has singleton; or a list of taxa that do not exist in the input 
#' phylogenetic profile.
#' @importFrom stringr regex
#' @importFrom stringr str_count
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @seealso \code{\link{getInputTaxaID}} for getting input taxon IDs,
#' \code{\link{ppTree}} for an example of input tree
#' @examples
#' data("ppTree", package="PhyloProfile")
#' checkNewick(ppTree, c("ncbi3702", "ncbi3711", "ncbi7029"))

checkNewick <- function(tree, inputTaxonID = NULL){
    if (missing(tree)) return("No tree given!")
    # get tree structure
    treeStruc <- gsub(stringr::regex("\\w"), "", as.character(tree$V1))

    open <- stringr::str_count(treeStruc, "\\(")
    close <- stringr::str_count(treeStruc, "\\)")
    comma <- stringr::str_count(treeStruc, "\\,")
    singleton <- stringr::str_count(treeStruc, "\\(\\)")

    if (singleton > 0) return(3) # tree contains singleton
    if (open != close) return(1) # missing parenthesis
    else {
        if ((comma - open) > 1 | (comma - open) < 0) {
            # return(2) # missing comma
        } else {
            # get list of tips
            nodeString <- gsub(regex("\\W+"), "#", as.character(tree$V1))
            nodeList <- unlist(strsplit(nodeString, "#"))
            # list of input taxa
            inputTaxa <- inputTaxonID

            # get missing taxa
            setdiff(inputTaxa, nodeList)
            missingTaxa <- setdiff(nodeList, inputTaxa)
            missingTaxa <- missingTaxa[lapply(missingTaxa, nchar) > 0]
            
            if (length(missingTaxa) > 0) {
                # contains taxa that not exist in main input
                return(paste(missingTaxa, collapse = "; "))
            } else return(0)
        }
    }
    return(0)
}

#' Create rooted tree from a taxonomy matrix
#' @export
#' @param df data frame contains taxonomy matrix used for generating tree 
#' (see distDf in example)
#' @param rootTaxon taxon used for rooting the taxonomy tree
#' @importFrom stats hclust
#' @importFrom ape as.phylo
#' @importFrom ape root
#' @return A rooted taxonomy tree as an object of class "phylo".
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @seealso \code{\link{taxa2dist}} for distance matrix generation from a
#' taxonomy matrix, \code{\link{getTaxonomyMatrix}} for getting taxonomy
#' matrix, \code{\link{ppTaxonomyMatrix}} for a demo taxonomy matrix data
#' @examples
#' data("ppTaxonomyMatrix", package = "PhyloProfile")
#' # prepare matrix for calculating distances
#' distDf <- subset(ppTaxonomyMatrix, select = -c(ncbiID, fullName))
#' row.names(distDf) <- distDf$abbrName
#' distDf <- distDf[, -1]
#' # create taxonomy tree rooted by ncbi10090
#' createRootedTree(distDf, "ncbi10090")

createRootedTree <- function(df, rootTaxon){
    if (missing(df)) return("No taxonomy matrix given!")
    # calculate distance matrix
    taxdis <- tryCatch(taxa2dist(df), error = function(e) e)
    # create tree
    tree <- ape::as.phylo(stats::hclust(taxdis))
    # root tree
    if (missing(rootTaxon)) rootTaxon = tree$tip.label[1]
    tree <- ape::root(tree, outgroup = rootTaxon, resolve.root = TRUE)
    # return
    return(tree)
}

#' Get sorted supertaxon list based on a rooted taxonomy tree
#' @export
#' @param tree an "phylo" object for a rooted taxonomy tree
#' @return A list of sorted taxa obtained the input taxonomy tree.
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @seealso \code{\link{ppTaxonomyMatrix}} for a demo taxonomy matrix data
#' @examples
#' data("ppTaxonomyMatrix", package = "PhyloProfile")
#' # prepare matrix for calculating distances
#' distDf <- subset(ppTaxonomyMatrix, select = -c(ncbiID, fullName))
#' row.names(distDf) <- distDf$abbrName
#' distDf <- distDf[, -1]
#' # create taxonomy tree rooted by ncbi10090
#' rootedTree <- createRootedTree(distDf, "ncbi10090")
#' # get taxon list sorted from tree
#' sortTaxaFromTree(rootedTree)

sortTaxaFromTree <- function(tree){
    if (missing(tree)) return("No tree given!")
    isTip <- tree$edge[, 2] <= length(tree$tip.label)
    orderedTips <- tree$edge[isTip, 2]
    taxonList <- rev(tree$tip.label[orderedTips])
    return(taxonList)
}
