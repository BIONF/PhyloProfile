#' Calculate the phylogenetic gene age from the phylogenetic profiles
#' @export
#' @usage estimateGeneAge(processedProfileData, rankName, refTaxon,
#'     var1CO, var2CO, percentCO)
#' @param processedProfileData dataframe contains the full processed
#' phylogenetic profiles (see ?fullProcessedProfile or ?parseInfoProfile)
#' @param rankName working taxonomy rank (e.g. "species", "genus", "family")
#' @param refTaxon reference taxon name (e.g. "Homo sapiens", "Homo" or
#' "Hominidae")
#' @param var1CO cutoff for var1
#' @param var2CO cutoff for var2
#' @param percentCO cutoff for percentage of species present in each
#' supertaxon
#' @return A dataframe contains estimated gene ages for the seed proteins.
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @seealso \code{\link{parseInfoProfile}} for creating a full processed
#' profile dataframe; \code{\link{getNameList}} and
#' \code{\link{getTaxonomyMatrix}} for getting taxonomy info,
#' \code{\link{fullProcessedProfile}} for a demo input dataframe
#' @examples
#' data("fullProcessedProfile", package="PhyloProfile")
#' rankName <- "class"
#' refTaxon <- "Mammalia"
#' processedProfileData <- fullProcessedProfile
#' var1Cutoff <- c(0, 1)
#' var2Cutoff <- c(0, 1)
#' percentCutoff <- c(0, 1)
#' estimateGeneAge(
#'     processedProfileData,
#'     rankName,
#'     refTaxon,
#'     var1Cutoff, var2Cutoff, percentCutoff
#' )

estimateGeneAge <- function(
    processedProfileData, rankName, refTaxon, var1CO, var2CO, percentCO
){
    rankList <- c("family", "class", "phylum", "kingdom", "superkingdom","root")
    # get selected (super)taxon ID
    taxList <- getNameList()
    superID <- taxList[
        taxList$fullName == refTaxon & taxList$rank == rankName,]$ncbiID
    # full non-duplicated taxonomy data
    Dt <- getTaxonomyMatrix(FALSE, NULL)
    # subset of taxonomy data, containing only ranks from rankList
    subDt <- Dt[, c("abbrName", rankList)]
    # get (super)taxa IDs for one of representative species
    firstLine <- Dt[Dt[, rankName] == superID, ][1, ]
    supFirstLine <- firstLine[, c("abbrName", rankList)]
    # compare each taxon IDs with selected taxon & create a "category" DF
    catList <- lapply(
        seq(nrow(subDt)), function (x) {
            cat <- subDt[x, ] %in% supFirstLine
            cat <- paste0(cat, collapse = "")})
    catDf <- data.frame(ncbiID = as.character(subDt$abbrName),
                        cat = do.call(rbind, catList), stringsAsFactors = FALSE)
    catDf$cat <- gsub("TRUE", "1", catDf$cat)
    catDf$cat <- gsub("FALSE", "0", catDf$cat)
    # get main input data
    mdData <- droplevels(processedProfileData)
    mdData <- mdData[, c("geneID","ncbiID","orthoID","var1","var2","presSpec")]
    # add "category" into mdData
    mdDtExt <- merge(mdData, catDf, by = "ncbiID", all.x = TRUE)
    mdDtExt$var1[mdDtExt$var1 == "NA" | is.na(mdDtExt$var1)] <- 0
    mdDtExt$var2[mdDtExt$var2 == "NA" | is.na(mdDtExt$var2)] <- 0
    # remove cat for "NA" orthologs and also for orthologs that dont fit cutoffs
    if (nrow(mdDtExt[mdDtExt$orthoID == "NA" | is.na(mdDtExt$orthoID), ]) > 0)
        mdDtExt[mdDtExt$orthoID == "NA" | is.na(mdDtExt$orthoID),]$cat <- NA
    mdDtExt <- mdDtExt[stats::complete.cases(mdDtExt), ]
    # filter by %specpres, var1, var2 ..
    mdDtExt <- subset(
        mdDtExt, mdDtExt$var1 >= var1CO[1] & mdDtExt$var1 <= var1CO[2]
        & mdDtExt$var2 >= var2CO[1] & mdDtExt$var2 <= var2CO[2]
        & mdDtExt$presSpec >=percentCO[1] & mdDtExt$presSpec<= percentCO[2])
    # get the furthest common taxon with selected taxon for each gene
    geneAgeDf <- as.data.frame(tapply(mdDtExt$cat, mdDtExt$geneID, min))
    data.table::setDT(geneAgeDf, keep.rownames = TRUE)[]
    data.table::setnames(geneAgeDf, seq_len(2), c("geneID", "cat"))  #col names
    row.names(geneAgeDf) <- NULL   # remove row names
    ### convert cat into geneAge
    geneAgeDf$age[geneAgeDf$cat == "0000001"] <- "07_LUCA"
    geneAgeDf$age[geneAgeDf$cat == "0000011" | geneAgeDf$cat == "0000010"] <-
        paste0(
            "06_", taxList$fullName[taxList$ncbiID == supFirstLine$superkingdom
                                    & taxList$rank == "superkingdom"])
    geneAgeDf$age[geneAgeDf$cat == "0000111"] <- paste0(
        "05_", taxList$fullName[
            taxList$ncbiID == supFirstLine$kingdom & taxList$rank == "kingdom"])
    geneAgeDf$age[geneAgeDf$cat == "0001111"] <- paste0(
        "04_", taxList$fullName[
            taxList$ncbiID == supFirstLine$phylum & taxList$rank == "phylum"])
    geneAgeDf$age[geneAgeDf$cat == "0011111"] <- paste0(
        "03_", taxList$fullName[
            taxList$ncbiID == supFirstLine$class & taxList$rank == "class"])
    geneAgeDf$age[geneAgeDf$cat == "0111111"] <- paste0(
        "02_", taxList$fullName[
            taxList$ncbiID == supFirstLine$family & taxList$rank == "family"])
    geneAgeDf$age[geneAgeDf$cat == "1111111"] <- paste0(
        "01_", taxList$fullName[
            taxList$fullName == refTaxon & taxList$rank == rankName])
    # return geneAge data frame
    geneAgeDf <- geneAgeDf[, c("geneID", "cat", "age")]
    geneAgeDf$age[is.na(geneAgeDf$age)] <- "Undef"
    return(geneAgeDf)
}

#' Create data for plotting gene ages
#' @param geneAgeDf data frame containing estimated gene ages for seed proteins
#' @return A dataframe for plotting gene age plot containing the absolute number
#' and percentage of genes for each calculated evolutionary ages and the
#' corresponding position for writting those number on the plot.
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @export
#' @seealso \code{\link{estimateGeneAge}}
#' @examples
#' geneAgeDf <- data.frame(
#' geneID = c("100136at6656", "100265at6656", "101621at6656", "103479at6656"),
#' cat = c("0000001", "0000011", "0000001", "0000011"),
#' age = c("07_LUCA", "06_Eukaryota", "07_LUCA", "06_Eukaryota")
#' )
#' geneAgePlotDf(geneAgeDf)

geneAgePlotDf <- function(geneAgeDf){
    plotDf <- plyr::count(geneAgeDf, c("age"))
    plotDf$percentage <- round(plotDf$freq / sum(plotDf$freq) * 100)
    plotDf$pos <- cumsum(plotDf$percentage) - (0.5 * plotDf$percentage)
    return(plotDf)
}

#' Create gene age plot
#' @param geneAgePlotDf data frame required for plotting gene age (see
#' ?geneAgePlotDf)
#' @param geneAgeText text size
#' @return A gene age distribution plot as a ggplot2 object
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @import ggplot2
#' @export
#' @seealso \code{\link{estimateGeneAge}} and \code{\link{geneAgePlotDf}}
#' @examples
#' geneAgePlotDf <- data.frame(
#'     age = "07_LUCA",
#'     frea = 2,
#'     percentage = 100,
#'     pos = 50
#' )
#' geneAgeText <- 1
#' createGeneAgePlot(geneAgePlotDf, geneAgeText)

createGeneAgePlot <- function(geneAgePlotDf, geneAgeText){
    age <- NULL
    percentage <- NULL
    pos <- NULL
    freq <- NULL

    p <- ggplot(geneAgePlotDf, aes(fill = age, y = percentage, x = 1)) +
        geom_bar(stat = "identity") +
        scale_y_reverse() +
        coord_flip() +
        theme_minimal()
    p <- p + geom_text(
        data = geneAgePlotDf,
        aes(x = 1, y = 100 - pos, label = paste0(freq, "\n", percentage, "%")),
        size = 4 * geneAgeText
    )
    p <- p + theme(
            legend.position = "bottom",
            legend.title = element_blank(),
            legend.text = element_text(size = 12 * geneAgeText),
            axis.title = element_blank(), axis.text = element_blank()
        ) +
        scale_fill_brewer(palette = "Spectral") +
        guides(
            fill = guide_legend(
                nrow = max(round(nrow(geneAgePlotDf) / 3, 0), 1), byrow = TRUE
            )
        )
    return(p)
}
