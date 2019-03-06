#' Get the dataframe with the significant genes
#' @export
#' @usage getSignificantGenes(inGroup, selectedGenesList, rank, var,
#'     useCommonAncestor, referenceTaxon, parameters, ncbiIDList,
#'     dataFull, rightFormatFeatures, domains)
#' @param inGroup list of taxa
#' @param selectedGenesList list of genes
#' @param rank selected taxonamy rank
#' @param var variable for which to  calculate the significance
#' @param useCommonAncestor boolean if the common anchestor should be used
#' @param referenceTaxon taxon which is used as reference
#' @param parameters contains "showPValue","highlightSignificant",
#' "significance", "var1ID", "var2ID", "xSizeGC", "ySizeGC",
#' "interestingFeatures", "angleGC", "legendGC", "legendSizeGC"
#' @param ncbiIDList list of ncbi ids
#' @param dataFull full processed main data
#' @param rightFormatFeatures boolean if the features have the right format
#' @param domains dataframe holding the domain input
#' @return dataframe with the significant genes
#' @importFrom stats p.adjust
#' @author Carla Mölbert (carla.moelbert@gmx.de)
#' @examples
#' data("mainLongRaw", package="phyloprofile")
#' data("fullProcessedProfile", package="phyloprofile")
#' inGroup <- c("Aeropyrum pernix", "Agrobacterium fabrum")
#' selectedGenesList <- "OG_1017"
#' rank <- "species"
#' var <- ""
#' useCommonAncestor <- FALSE
#' referenceTaxon <- "Homo sapiens"
#' ncbiIDList <- getInputTaxaID(mainLongRaw)
#' dataFull <- fullProcessedProfile
#' rightFormatFeatures <- TRUE
#' domainFile <- system.file(
#'     "extdata", "domainFiles/OG_1017.domains",
#'     package = "phyloprofile", mustWork = TRUE
#' )
#' domains <- parseDomainInput("OG_1017", domainFile, "file")
#' parameters <- list(
#'     "showPValue" = TRUE,
#'     "highlightSignificant" = FALSE,
#'     "significance" = 0.05,
#'     "var1ID" = "1st variable",
#'     "var2ID" = "2nd variable",
#'     "xSizeGC" = 9,
#'     "ySizeGC" = 9,
#'     "interestingFeatures" = NULL,
#'     "angleGC" = 90,
#'     "legendGC" = "none",
#'     "legendSizeGC" = 9,
#'     "pValuesSize" = 9
#' )
#' getSignificantGenes(
#'     inGroup,
#'     selectedGenesList,
#'     rank,
#'     var,
#'     useCommonAncestor,
#'     referenceTaxon,
#'     parameters,
#'     ncbiIDList,
#'     dataFull,
#'     rightFormatFeatures,
#'     domains
#' )

getSignificantGenes <- function(
    inGroup,
    selectedGenesList,
    rank,
    var,
    useCommonAncestor,
    referenceTaxon,
    parameters,
    ncbiIDList,
    dataFull,
    rightFormatFeatures,
    domains
){
    if (is.null(inGroup) | length(selectedGenesList) == 0) return()
    significanceLevel <- parameters$significance

    nameList <- getNameList() # load name List
    taxaList <- getTaxonomyMatrix(FALSE, ncbiIDList) # load unsorted taxa

    # Get the rank and the in-group
    # if there is more than one element in the inGroup -> use common anchestor
    if (useCommonAncestor == TRUE) {
        ancestor <- getCommonAncestor(inGroup, rank,
                                        nameList, taxaList, ncbiIDList)
        if (is.null(ancestor)) return("No common ancestor found")
        inGroup <- ancestor[1]
        rank <- ancestor[2]
    } else {
        # rank <- substring(rank, 4)
        rank <- rank
    }
    if (is.na(rank)) return("No common ancestor found")

    # provide the empty data frame
    if (var == "Both") {
        significantGenesDf <- data.frame(
            geneID = character(),
            inGroup = I(list()),
            outGroup = I(list()),
            pvalues1 = I(list()),
            pvalues2 = I(list()),
            features = I(list()),
            databases = I(list()))
    } else {
        significantGenesDf <- data.frame(
            geneID = character(),
            inGroup = I(list()),
            outGroup = I(list()),
            pvalues = I(list()),
            features = I(list()),
            databases = I(list()))
    }

    # Get the list of genes to look at
    if (is.element("all", selectedGenesList)) {
        genes <- dataFull$geneID
        genes <- genes[!duplicated(genes)]
    } else {
        genes <- selectedGenesList
    }
    genes <- sort(genes)

    # Subset depending on the rank and the inGroup
    selectedSubset <- getSelectedSubset(rank, inGroup, nameList, taxaList)
    selectedSubset <- subset(
        selectedSubset, !selectedSubset$fullName == referenceTaxon
    )

    # Check for each gene if it is significant
    for (gene in genes) {
        print(paste("Analyzing the distribution of", gene, "..."))
        # Processing the dataframes for in- and out-group
        selectedGeneDf <- subset(dataFull, dataFull$geneID == gene)

        inGroupDf <- {
            subset(
                selectedGeneDf,
                selectedGeneDf$abbrName %in% selectedSubset$abbrName
            )
        }
        outGroupDf <- {
            subset(
                selectedGeneDf,
                !(selectedGeneDf$abbrName %in% selectedSubset$abbrName)
            )
        }
        outGroupDf <- {
            subset(outGroupDf, !outGroupDf$fullName == referenceTaxon)
        }

        # Generate and check the pValues for the gene
        pvalue <- getPValues(inGroupDf, outGroupDf, var, gene, parameters)
        newRow <- data.frame(
            geneID = gene,
            inGroup = NA,
            outGroup = NA,
            pvalues = NA,
            features = NA
        )
        newRow$inGroup <- list(inGroupDf)
        newRow$outGroup <- list(outGroupDf)

        if (var == "Both") {
            newRow$pvalues1 <- pvalue[1]
            newRow$pvalues2 <- pvalue[2]
        } else {
            newRow$pvalues <- pvalue
        }

        features  <- getFeatures(gene, domains)
        newRow$features <- list(features)
        if (rightFormatFeatures) {
            newRow$databases <- list(getPrefixFeatures(features))
        }
        significantGenesDf <- rbind(significantGenesDf, newRow)
    }

    if (var == "Both") {
        significantGenesDf$pvalues1 <- {
            p.adjust(
                significantGenesDf$pvalues1, method = "holm",
                n = length(significantGenesDf$pvalues1)
            )
        }

        significantGenesDf$pvalues2 <- {
            p.adjust(
                significantGenesDf$pvalues2, method = "holm",
                n = length(significantGenesDf$pvalues2)
            )
        }

        significantGenesDf <-
            significantGenesDf[
                significantGenesDf$pvalues1 <= significanceLevel
                | significantGenesDf$pvalues2 <= significanceLevel ,
            ]

        significantGenesDf <-
            significantGenesDf[
                !is.na(significantGenesDf$pvalues1) |
                !is.na(significantGenesDf$pvalues1),
            ]
    } else {
        significantGenesDf$pvalues <- {
            p.adjust(
                significantGenesDf$pvalues, method = "holm",
                n = length(significantGenesDf$pvalues)
            )
        }

        significantGenesDf <- {
            significantGenesDf[
                significantGenesDf$pvalues <= significanceLevel,
            ]
        }

        significantGenesDf <- {
            significantGenesDf[!is.na(significantGenesDf$pvalues) ,]
        }
    }

    # return the significant genes
    if (nrow(significantGenesDf) != 0) {
        significantGenesDf$var <- var
        significantGenesDf$rank <- rank
        return(significantGenesDf)
    } else {
        return("No candidate genes found")
    }
}

#' Generate the inGroup
#' @export
#' @usage getCommonAncestor(inGroup, rank, nameList, selectedInGroup,
#'     ncbiIDList)
#' @param inGroup list of taxa
#' @param rank selected rank
#' @param nameList contains "ncbiID", "fullName", "rank", "parentID"
#' @param selectedInGroup contains "abbrName, "ncbiID", fullName", "strain",
#' "genus"..
#' @param ncbiIDList list of input taxon IDs
#' @return common anchestor
#' @author Carla Mölbert (carla.moelbert@gmx.de)
#' @examples
#' inGroup <- c("Aeropyrum pernix", "Agrobacterium fabrum")
#' rank <- "species"
#' nameListFile <- system.file(
#'     "extdata", "data/taxonNamesFull.txt",
#'     package = "phyloprofile", mustWork = TRUE
#' )
#' nameList <- as.data.frame(data.table::fread(nameListFile))
#' selectedInGroup <- getTaxonomyMatrix(FALSE, NULL)
#' ncbiIDList <- c(56636, 1176649)
#' getCommonAncestor(
#'     inGroup,
#'     rank,
#'     nameList,
#'     selectedInGroup,
#'     ncbiIDList
#' )

getCommonAncestor <- function(inGroup,
                                rank,
                                nameList,
                                selectedInGroup,
                                ncbiIDList){

    allRanks <- getTaxonomyRanks()

    selectedInGroup <- {
        selectedInGroup[!duplicated(selectedInGroup), ]
    }

    # ranks were all elements of the inGroup might be in the same taxon
    # possibleRanks <- allRanks[allRanks >= rank]
    possibleRanks <- allRanks[match(rank, allRanks):(length(allRanks)-1)]
    position <-  1
    if (length(inGroup) == 1) rank <- rank #substring(rank, 4)

    # find the common ancestor of all taxa in the inGroup
    while (length(inGroup) > 1 & position < length(possibleRanks)) {

        currentRank <- as.character(possibleRanks[position][1])
        nextRank <- as.character(possibleRanks[position + 1][1])

        # dataframe with all elements with fitting rank
        dfInGroup <- subset(nameList, nameList$rank == currentRank)

        # subset of dfInGroup  with elements that belong to the inGroup
        dfInGroup <- subset(dfInGroup, dfInGroup$fullName %in% inGroup)

        # get all elements which could belong to the in-group
        possibleInGroup <- subset(selectedInGroup,
                                    select = c(currentRank, nextRank))
        possibleInGroup <- {
            possibleInGroup[
                possibleInGroup[,currentRank] %in% dfInGroup$ncbiID,
            ]
        }
        possibleInGroup <- possibleInGroup[!duplicated(possibleInGroup), ]

        # only consider elements that have the next higher rank
        subsetNextRank <- taxaSelectGC(nextRank, ncbiIDList)
        subsetNextRank <- subsetNextRank[!duplicated(subsetNextRank), ]
        subsetNextRank <- {
            subset(
                subsetNextRank,
                subsetNextRank$ncbiID %in% possibleInGroup[, nextRank]
            )
        }
        inGroup <- subsetNextRank$fullName
        position <- position + 1
        rank <- nextRank
    }

    # Return the in-group and the rank
    if (position > length(possibleRanks)) return()
    return(c(inGroup, rank))
}


#' Get the subset depending on the choosen rank
#' @param rank selected rank
#' @param inGroup list of taxa
#' @param nameList contains "ncbiID", "fullName", "rank", "parentID"
#' @param taxaList contains "abbrName, "ncbiID", fullName", "strain", "genus"
#' @return list of prefixes for the features
#' @author Carla Mölbert (carla.moelbert@gmx.de)
getSelectedSubset <- function(rank, inGroup, nameList, taxaList){
    # Look if the fullName is in the inGroup
    nameList$fullName <- as.character(nameList$fullName)
    nameListRank <- subset(nameList, nameList$rank == rank)
    inGroupSubset <- subset(nameList, nameList$fullName %in% inGroup)

    # Look if it has the right rank
    selectedSubset <- taxaList[taxaList[, rank] %in% inGroupSubset$ncbiID,]

    return(selectedSubset)
}

#' Decide if the gene is significant
#' @param inGroup contains "supertaxon", "geneID", "ncbiID", "orthoID",
#' "var1", "var2", "paralog", "abbrName", "taxonID", "fullname",
#' "supertaxonID", "rank", "category", "presSpec", "mVar1", "mVar2"
#' @param outGroup  as in-group but with information containing the out-group
#' @param variable variable(s) to claculate the plots for
#' @param gene gene to calculate the p-values for
#' @param parameters contains "showPValue","highlightSignificant",
#' "significance", "var1ID", "var2ID", "xSizeGC", "ySizeGC",
#' "interestingFeatures", "angleGC", "legendGC", "legendSizeGC"
#' @return return the pvalues
#' @author Carla Mölbert (carla.moelbert@gmx.de)

getPValues <- function(inGroup, outGroup, variable, gene, parameters){
    significanceLevel <- parameters$significance

    # get the p-values for both variables
    if (variable == "Both") {
        var1 <- parameters$var1
        var2 <- parameters$var2

        pvalues1 <- calculatePValue(
            inGroup$var1,
            outGroup$var1,
            significanceLevel
        )

        pvalues2 <- calculatePValue(
            inGroup$var2,
            outGroup$var2,
            significanceLevel
        )
        pvalues <- list(pvalues1, pvalues2)
        return(pvalues)
    }
    # get the p-values for one variable
    else {
        # Check which variable is selected and get the pValues
        if (variable == parameters$var1ID) {
            pvalues <- calculatePValue(
                inGroup$var1,
                outGroup$var1,
                significanceLevel
            )
        } else {
            pvalues <- calculatePValue(
                inGroup$var2,
                outGroup$var2,
                significanceLevel
            )
        }

        return(pvalues)
    }
}

#' calculate the pValues
#' @param varIn list of values for the variable concerning the in-group
#' @param varOut list of values for the variable concerning the out-group
#' @param significanceLevel significant cutoff for statistical test
#' @return return the pvalues
#' @importFrom stats ks.test
#' @importFrom stats wilcox.test
#' @author Carla Mölbert (carla.moelbert@gmx.de)

calculatePValue <- function(varIn, varOut, significanceLevel){
    # delete all entrys that are NA
    varIn <- varIn[!is.na(varIn)]
    varOut <- varOut[!is.na(varOut)]

    # if there is no data in one of the groups the p-value is NULL
    if (length(varIn) == 0) return(NA)
    else if (length(varOut) == 0) return(NA)
    else {
        # * Kolmogorov-Smirnov Test
        # H0 : The two samples have the same distribution
        ks <- suppressWarnings(
            ks.test(unique(varIn), unique(varOut), exact = FALSE)
        )

        pValue <- ks$p.value # probabilitiy to recet H0, if it is correct

        if (pValue < significanceLevel) pvalue <- c(pValue)

        else {
            # * Wilcoxon-Mann-Whitney Test
            # H0: the samples have the same location parameters

            wilcox <- suppressWarnings(
                wilcox.test(varIn,
                            varOut,
                            alternative = "two.sided",
                            #exact = FALSE,
                            paired = FALSE)
            )
            pValueWilcox <- wilcox$p.value
            # pvalue <- c(pValue, pValueWilcox)
            pvalue <- c(pValueWilcox)
        }

        # perm <- jmuOutlier:: perm.test(
        #     unique(varIn), unique(varOut),
        #     alternative = c("two.sided"),
        #     mu = 0, # Hypothesis: Samples do not differ
        #     paired = FALSE,
        #     all.perms = TRUE, # Tries to get a exact p value
        #     num.sim = 1000000,
        #     plot = FALSE, # does not plot
        #     stat = mean
        # ) # compaires the means of the distributions
        # pvalue <- perm$p.value

        # return the calculated pvalues ----------------------------------------
        return(pvalue)
    }
}

#' get the list with all the features in the gene
#' @param selectedGene gene to get the feartures for
#' @param domains contains "seedID", "orthoID", "feature", "start",   "end"
#' @return dataframe for the specific gene containing "seedID",  "orthoID",
#' "feature", "start",   "end"
#' @author Carla Mölbert (carla.moelbert@gmx.de)
getFeatures <- function(selectedGene, domains){
    subsetDomains <- {
        subset(
            domains,
            substr(
                domains$seedID,
                1,
                nchar(as.character(selectedGene))
            ) == selectedGene
        )
    }
    subsetDomains <- subsetDomains[!duplicated(subsetDomains), ]
    return(subsetDomains)
}

#' Get the database for each feature in a specific gene
#' @param data contains "seedID", "orthoID", "feature", "start", "end"
#' @return list of prefixes for the features
#' @author Carla Mölbert (carla.moelbert@gmx.de)
getPrefixFeatures <- function(data){
    features <- data$feature
    choices <- gsub("_.*", "", features)
    choices <- choices[!duplicated(choices)]
    return(choices)
}

#' print list of available taxa
#' @param rankSelectGC rank selected for group compariosn
#' @param inputTaxonID contains "seedID",  "orthoID", "feature", "start", "end"
#' @return avilable taxa containing "ncbiID", "fullName", "rank", "parentID"
#' @author Carla Mölbert (carla.moelbert@gmx.de)
taxaSelectGC <- function(rankSelectGC, inputTaxonID){
    # if there is no rank set, there can not be any available taxa
    if (length(rankSelectGC) == 0) return()
    else{
        # load list of unsorted taxa
        if (is.null(inputTaxonID)) dt <- getTaxonomyMatrix(
            FALSE, inputTaxonID
        )
        else dt <- getTaxonomyMatrix(TRUE, inputTaxonID)

        # load list of taxon name
        nameList <- phyloprofile::getNameList()

        # get rank name from rankSelect
        if (substr(rankSelectGC,3,3) == "_") {
            # rankName <- substr(rankSelectGC, 4, nchar(rankSelectGC))
            rankName <- rankSelectGC
        }
        else rankName <- rankSelectGC

        choice <- as.data.frame
        choice <- rbind(dt[rankName])
        colnames(choice) <- "ncbiID"
        choice <- merge(choice, nameList, by = "ncbiID", all = FALSE)
        return(choice)
    }
}

#' Create a list with all main taxanomy ranks
#' @export
#' @return A list of all main ranks (from strain to superkingdom)
#' @author Carla Mölbert (carla.moelbert@gmx.de)
#' @examples
#' getTaxonomyRanks()

getTaxonomyRanks <- function(){
    allRanks <- list(
        "Strain " = "strain",
        "Species" = "species",
        "Genus" = "genus",
        "Family" = "family",
        "Order" = "order",
        "Class" = "class",
        "Phylum" = "phylum",
        "Kingdom" = "kingdom",
        "Superkingdom" = "superkingdom",
        "unselected" = ""
    )

    return(allRanks)
}
