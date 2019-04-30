#' Calculate the phylogenetic gene age from the phylogenetic profiles
#' @export
#' @usage estimateGeneAge(processedProfileData, rankName, refTaxon,
#'     var1Cutoff, var2Cutoff, percentCutoff)
#' @param processedProfileData dataframe contains the full processed
#' phylogenetic profiles (see ?fullProcessedProfile or ?parseInfoProfile)
#' @param rankName working taxonomy rank (e.g. "species", "genus", "family")
#' @param refTaxon reference taxon name (e.g. "Homo sapiens", "Homo" or
#' "Hominidae")
#' @param var1Cutoff cutoff for var1
#' @param var2Cutoff cutoff for var2
#' @param percentCutoff cutoff for percentage of species present in each
#' supertaxon
#' @return A dataframe contains estimated gene ages for the seed proteins.
#' @importFrom data.table setDT
#' @importFrom data.table setnames
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
    processedProfileData,
    rankName, refTaxon,
    var1Cutoff, var2Cutoff, percentCutoff
){
    rankList <- c(
        "family", "class", "phylum", "kingdom", "superkingdom", "root"
    )

    # get selected (super)taxon ID
    taxaList <- getNameList()
    superID <- 
        taxaList[
            taxaList$fullName == refTaxon & taxaList$rank == rankName,
        ]$ncbiID

    # full non-duplicated taxonomy data
    Dt <- getTaxonomyMatrix(FALSE, NULL)

    # subset of taxonomy data, containing only ranks from rankList
    subDt <- Dt[, c("abbrName", rankList)]

    # get (super)taxa IDs for one of representative species
    # get all taxon info for 1 representative
    firstLine <- Dt[Dt[, rankName] == superID, ][1, ]
    subFirstLine <- firstLine[, c("abbrName", rankList)]

    # compare each taxon ncbi IDs with selected taxon
    # and create a "category" data frame
    catList <- lapply(
        seq(nrow(subDt)),
        function (x) {
            cat <- subDt[x, ] %in% subFirstLine
            cat <- paste0(cat, collapse = "")
            cat <- gsub("TRUE", "1", cat)
            gsub("FALSE", "0", cat)
        }
    )
    
    catDf <- data.frame(
        ncbiID = as.character(subDt$abbrName),
        cat = do.call(rbind, catList),
        stringsAsFactors = FALSE
    )

    # get main input data
    mdData <- droplevels(processedProfileData)
    mdData <- mdData[, c(
        "geneID", "ncbiID", "orthoID", "var1", "var2", "presSpec"
    )]

    ### add "category" into mdData
    mdDataExtended <- merge(mdData, catDf, by = "ncbiID", all.x = TRUE)

    mdDataExtended$var1[mdDataExtended$var1 == "NA"
                        | is.na(mdDataExtended$var1)] <- 0
    mdDataExtended$var2[mdDataExtended$var2 == "NA"
                        | is.na(mdDataExtended$var2)] <- 0

    # remove cat for "NA" orthologs
    # and also for orthologs that do not fit cutoffs
    if (nrow(mdDataExtended[mdDataExtended$orthoID == "NA"
        | is.na(mdDataExtended$orthoID), ]) > 0) {
        mdDataExtended[mdDataExtended$orthoID == "NA"
        | is.na(mdDataExtended$orthoID), ]$cat <- NA
    }

    mdDataExtended <- mdDataExtended[complete.cases(mdDataExtended), ]

    # filter by %specpres, var1, var2 ..
    mdDataExtended$cat[mdDataExtended$var1 < var1Cutoff[1]] <- NA
    mdDataExtended$cat[mdDataExtended$var1 > var1Cutoff[2]] <- NA
    mdDataExtended$cat[mdDataExtended$var2 < var2Cutoff[1]] <- NA
    mdDataExtended$cat[mdDataExtended$var2 > var2Cutoff[2]] <- NA
    mdDataExtended$cat[mdDataExtended$presSpec < percentCutoff[1]] <- NA
    mdDataExtended$cat[mdDataExtended$presSpec > percentCutoff[2]] <- NA

    mdDataExtended <- mdDataExtended[complete.cases(mdDataExtended), ]

    ### get the furthest common taxon with selected taxon for each gene
    geneAgeDf <- as.data.frame(
        tapply(mdDataExtended$cat, mdDataExtended$geneID, min)
    )

    setDT(geneAgeDf, keep.rownames = TRUE)[]
    setnames(geneAgeDf, seq_len(2), c("geneID", "cat"))  # rename columns
    row.names(geneAgeDf) <- NULL   # remove row names

    ### convert cat into geneAge
    geneAgeDf$age[geneAgeDf$cat == "0000001"] <- "07_LUCA"
    geneAgeDf$age[geneAgeDf$cat == "0000011" | geneAgeDf$cat == "0000010"] <-
        paste0(
            "06_",
            as.character(
                taxaList$fullName[
                    taxaList$ncbiID == subFirstLine$superkingdom
                    & taxaList$rank == "superkingdom"
                ]
            )
        )

    geneAgeDf$age[geneAgeDf$cat == "0000111"] <-
        paste0(
            "05_",
            as.character(
                taxaList$fullName[
                    taxaList$ncbiID == subFirstLine$kingdom
                    & taxaList$rank == "kingdom"
                ]
            )
        )

    geneAgeDf$age[geneAgeDf$cat == "0001111"] <-
        paste0(
            "04_",
            as.character(
                taxaList$fullName[
                    taxaList$ncbiID == subFirstLine$phylum
                    & taxaList$rank == "phylum"
                ]
            )
        )

    geneAgeDf$age[geneAgeDf$cat == "0011111"] <-
        paste0(
            "03_",
            as.character(
                taxaList$fullName[
                    taxaList$ncbiID == subFirstLine$class
                    & taxaList$rank == "class"
                ]
            )
        )

    geneAgeDf$age[geneAgeDf$cat == "0111111"] <-
        paste0(
            "02_",
            as.character(
                taxaList$fullName[
                    taxaList$ncbiID == subFirstLine$family
                    & taxaList$rank == "family"
                ]
            )
        )

    geneAgeDf$age[geneAgeDf$cat == "1111111"] <-
        paste0(
            "01_",
            as.character(
                taxaList$fullName[
                    taxaList$fullName == refTaxon
                    & taxaList$rank == rankName
                ]
            )
        )

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
#' geneID = c("OG_1017", "OG_1019"),
#' cat = c("0000001", "0000001"),
#' age = c("07_LUCA", "07_LUCA")
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
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 scale_y_reverse
#' @importFrom ggplot2 coord_flip
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 scale_fill_brewer
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 guide_legend
#' @export
#' @seealso \code{\link{estimateGeneAge}} and \code{\link{geneAgePlotDf}}
#' @examples
#' geneAgePlotDf <- data.frame(
#'     age = "OG_1017",
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
