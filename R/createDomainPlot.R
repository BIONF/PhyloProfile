#' Create protein's domain architecure plot
#' @description Create architecture plot for both seed and orthologous protein.
#' If domains of ortholog are missing, only architecture of seed protein will
#' be plotted. NOTE: seed protein ID is the one being shown in the profile plot,
#' which normally is also the orthologous group ID.
#' @export
#' @usage createArchiPlot(info = NULL, domainDf = NULL, labelArchiSize = 12,
#'     titleArchiSize = 12)
#' @param info a list contains seed and ortholog's IDs
#' @param domainDf dataframe contains domain info for the seed and ortholog.
#' This including the seed ID, orthologs IDs, sequence lengths, feature names,
#' start and end positions, feature weights (optional) and the status to
#' determine if that feature is important for comparison the architecture
#' between 2 proteins* (e.g. seed protein vs ortholog) (optional).
#' @param labelArchiSize lable size (in px). Default = 12.
#' @param titleArchiSize title size (in px). Default = 12.
#' @importFrom gridExtra arrangeGrob
#' @import ggplot2
#' @return A domain plot as arrangeGrob object. Use grid::grid.draw(plot) to
#' render.
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @seealso \code{\link{domainPlotting}}, \code{\link{sortDomains}},
#' \code{\link{parseDomainInput}}, \code{\link{getQualColForVector}}
#' @examples
#' seedID <- "101621at6656"
#' orthoID <- "101621at6656|AGRPL@224129@0|224129_0:001955|1"
#' info <- c(seedID, orthoID)
#' domainFile <- system.file(
#'     "extdata", "domainFiles/101621at6656.domains",
#'     package = "PhyloProfile", mustWork = TRUE
#' )
#' domainDf <- parseDomainInput(seedID, domainFile, "file")
#' plot <- createArchiPlot(info, domainDf, 9, 9)
#' grid::grid.draw(plot)

createArchiPlot <- function(
    info = NULL,
    domainDf = NULL,
    labelArchiSize = 12, titleArchiSize = 12
){
    if (is.null(info) | is.null(domainDf)) return(ggplot() + theme_void())
    orthoID <- NULL

    # info
    group <- as.character(info[1])
    ortho <- as.character(info[2])

    # get sub dataframe based on selected groupID and orthoID
    ortho <- gsub("\\|", ":", ortho)
    grepID <- paste(group, "#", ortho, sep = "")

    subdomainDf <- domainDf[grep(grepID, domainDf$seedID), ]
    subdomainDf$feature <- as.character(subdomainDf$feature)

    if (nrow(subdomainDf) < 1) return(paste0("ERR-0"))
    else {
        # ortho domains df
		orthoDf <- subdomainDf[subdomainDf$orthoID == ortho,]
        # seed domains df
        seedDf <- subdomainDf[subdomainDf$orthoID != ortho,]
        if (nrow(seedDf) == 0) seedDf <- orthoDf

        seed <- as.character(seedDf$orthoID[1])

        # return ERR-0 if seedDf and orthoDf are empty
        if (nrow(seedDf) == 0) return(paste0("ERR-0"))

        # change order of one dataframe's features
        # based on order of other df's features
        if (length(orthoDf$feature) < length(seedDf$feature)) {
            orderedOrthoDf <- orthoDf[order(orthoDf$feature), ]
            orderedSeedDf <- sortDomains(orderedOrthoDf, seedDf)
        } else {
            orderedSeedDf <- seedDf[order(seedDf$feature), ]
            orderedOrthoDf <- sortDomains(orderedSeedDf, orthoDf)
        }

        # join weight values and feature names
        if ("weight" %in% colnames(orderedOrthoDf)) {
            orderedOrthoDf$yLabel <- paste0(
                orderedOrthoDf$feature,
                " (",
                round(orderedOrthoDf$weight, 2),
                ")"
            )
        } else {
            orderedOrthoDf$yLabel <- orderedOrthoDf$feature
        }
        if ("weight" %in% colnames(orderedSeedDf)) {
            orderedSeedDf$yLabel <- paste0(
                orderedSeedDf$feature,
                " (",
                round(orderedSeedDf$weight, 2),
                ")"
            )
        } else {
            orderedSeedDf$yLabel <- orderedSeedDf$feature
        }

        # create color scheme for all features
        # the same features in seed & ortholog will have the same colors
        featureSeed <- levels(as.factor(orderedSeedDf$feature))
        featureOrtho <- levels(as.factor(orderedOrthoDf$feature))
        allFeatures <- c(featureSeed, featureOrtho)
        allColors <- getQualColForVector(allFeatures)
        colorScheme <- structure(allColors, .Names = allFeatures)

        # plotting
        sep <- "|"

        if ("length" %in% colnames(subdomainDf)) {
            plotOrtho <- domainPlotting(
                orderedOrthoDf,
                ortho,
                sep,
                labelArchiSize,
                titleArchiSize,
                min(subdomainDf$start),
                max(c(subdomainDf$end, subdomainDf$length)),
                colorScheme
            )
            plotSeed <- domainPlotting(
                orderedSeedDf,
                seed,
                sep,
                labelArchiSize,
                titleArchiSize,
                min(subdomainDf$start),
                max(c(subdomainDf$end, subdomainDf$length)),
                colorScheme
            )

        } else{
            plotOrtho <- domainPlotting(
                orderedOrthoDf,
                ortho,
                sep,
                labelArchiSize,
                titleArchiSize,
                min(subdomainDf$start),
                max(subdomainDf$end),
                colorScheme
            )
            plotSeed <- domainPlotting(
                orderedSeedDf,
                seed,
                sep,
                labelArchiSize,
                titleArchiSize,
                min(subdomainDf$start),
                max(subdomainDf$end),
                colorScheme
            )
        }

        # grid.arrange(plotSeed,plotOrtho,ncol=1)
        if (ortho == seed) {
            g <- gridExtra::arrangeGrob(plotSeed, ncol = 1)
        } else {
            seedHeight <- length(levels(as.factor(orderedSeedDf$feature)))
            orthoHeight <- length(levels(as.factor(orderedOrthoDf$feature)))

            g <- gridExtra::arrangeGrob(plotSeed, plotOrtho, ncol = 1,
                                        heights = c(seedHeight, orthoHeight))
        }
        return(g)
    }
}

#' Create architecure plot for a single protein
#' @usage domainPlotting(df, geneID = "GeneID", sep = "|", labelSize = 12,
#'     titleSize = 12, minStart = NULL, maxEnd = NULL, colorScheme)
#' @param df domain dataframe for ploting containing the seed ID, ortholog ID,
#' ortholog sequence length, feature names, start and end positions,
#' feature weights (optional) and the status to determine if that feature is
#' important for comparison the architecture between 2 proteins* (e.g. seed
#' protein vs ortholog) (optional).
#' @param geneID ID of seed or orthologous protein
#' @param sep separate indicator for title. Default = "|".
#' @param labelSize lable size. Default = 12.
#' @param titleSize title size. Default = 12.
#' @param minStart the smallest start position of all domains
#' @param maxEnd the highest stop position of all domains
#' @param colorScheme color scheme for all domain types
#' @return Domain plot of a single protein as a ggplot object.
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @seealso \code{\link{getQualColForVector}},
#' \code{\link{parseDomainInput}}
#' @import ggplot2
#' @examples
#' \dontrun{
#' # get domain data
#' domainFile <- system.file(
#'     "extdata", "domainFiles/101621at6656.domains",
#'     package = "PhyloProfile", mustWork = TRUE
#' )
#' domainDf <- parseDomainInput(seedID, domainFile, "file")
#' df <- domainDf[
#'     domainDf$orthoID == "101621at6656|AGRPL@224129@0|224129_0:001955|1",]
#' # create color scheme for all domain types
#' allFeatures <- levels(as.factor(df$feature))
#' allColors <- getQualColForVector(allFeatures)
#' colorScheme <- structure(
#'     allColors,
#'     .Names = allFeatures
#' )
#' # other parameters
#' geneID <- "AGRPL@224129@0|224129_0:001955|1"
#' sep <- "|"
#' labelSize <- 9
#' titleSize <- 9
#' minStart <- min(df$start)
#' maxEnd <- max(df$end)
#' # do plotting
#' domainPlotting(
#'     df,
#'     geneID,
#'     sep,
#'     labelSize, titleSize,
#'     minStart, maxEnd,
#'     colorScheme
#' )
#' }

domainPlotting <- function(
    df = NULL,
    geneID = "GeneID",
    sep = "|",
    labelSize = 12, titleSize = 12,
    minStart = NULL, maxEnd = NULL,
    colorScheme = NULL
){
    feature <- NULL
    end <- NULL
    start <- NULL
    # parse parameters
    if (is.null(df)) return(ggplot() + theme_void())
    if (is.null(minStart)) minStart <- min(df$start)
    if (is.null(maxEnd)) maxEnd <- max(df$end)
    if (is.null(colorScheme)) {
        allFeatures <- levels(as.factor(df$feature))
        allColors <- getQualColForVector(allFeatures)
        colorScheme <- structure(
            allColors,
            .Names = allFeatures
        )
    }

    gg <- ggplot(df, aes(y = feature, x = end, color = as.factor(feature))) +
        geom_segment(
            data = df,
            aes(y = feature, yend = feature, x = minStart, xend = maxEnd),
            color = "white",
            size = 0
        ) +
        scale_color_manual(values = colorScheme)

    # draw lines for representing sequence length
    if ("length" %in% colnames(df)) {
        gg <- gg + geom_segment(
            data = df,
            aes(x = 0, xend = length, y = feature, yend = feature),
            size = 1,
            color = "#b2b2b2"
        )
    }

    # draw line and points
    gg <- gg + geom_segment(
        data = df, aes(x = start, xend = end, y = feature, yend = feature),
        size = 1.5
    )
    gg <- gg + geom_point(
        data = df, aes(y = feature, x = start),
        color = "#b2b2b2", size = 3, shape = 3
    )
    gg <- gg + geom_point(
        data = df, aes(y = feature, x = end),
        color = "#edae52", size = 3, shape = 5
    )

    # draw dashed line for domain path
    gg <- gg + geom_segment(
        data = df[df$path == "Y", ],
        aes(x = start, xend = end, y = feature, yend = feature),
        size = 3, linetype = "dashed"
    )

    # theme format
    titleMod <- gsub(":", sep, geneID)
    gg <- gg + scale_y_discrete(
        expand = c(0.075, 0), breaks = df$feature, labels = df$yLabel
    )
    gg <- gg + labs(title = paste0(titleMod), y = "Feature")
    gg <- gg + theme_minimal()
    gg <- gg + theme(panel.border = element_blank())
    gg <- gg + theme(axis.ticks = element_blank())
    gg <- gg + theme(
        plot.title = element_text(face = "bold", size = titleSize)
    )
    gg <- gg + theme(plot.title = element_text(hjust = 0.5))
    gg <- gg + theme(
        legend.position = "none", axis.title.x = element_blank(),
        axis.text.y = element_text(size = labelSize),
        axis.title.y = element_text(size = labelSize),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank()
    )
    # return plot
    return(gg)
}

#' Sort one domain dataframe based on the other domain dataframe
#' @description Sort domain dataframe of one protein (either seed or ortholog)
#' based on the dataframe of the its paired protein, in order to bring the
#' common domain feature in the same order which make it easy for comparing.
#' @param seedDf data of seed protein
#' @param orthoDf data of ortholog protein
#' @return Dataframe contains sorted domain list.
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @examples
#' \dontrun{
#' # get domain data
#' domainFile <- system.file(
#'     "extdata", "domainFiles/101621at6656.domains",
#'     package = "PhyloProfile", mustWork = TRUE
#' )
#' domainDf <- parseDomainInput(seedID, domainFile, "file")
#' # get seedDf and orthoDf
#' subDf <- domainDf[
#'     domainDf$seedID ==
#'     "101621at6656#101621at6656:AGRPL@224129@0:224129_0:001955:1",]
#' orthoDf <- subDf[subDf$orthoID == "101621at6656:DROME@7227@1:Q9VG04",]
#' seedDf <- subDf[subDf$orthoID != "101621at6656:DROME@7227@1:Q9VG04",]
#' # sort
#' sortDomains(seedDf, orthoDf)
#' }

sortDomains <- function(seedDf, orthoDf){
    if (is.null(seedDf) | is.null(orthoDf)) return()
    orderNo <- NULL
    # get list of features in seedDf
    featureList <- as.data.frame(levels(as.factor(seedDf$feature)))
    colnames(featureList) <- c("feature")
    # and add order number to each feature
    featureList$orderNo <- seq(length(featureList$feature))

    # merge those info to orthoDf
    orderedOrthoDf <- merge(orthoDf, featureList, all.x = TRUE)

    # sort orthoDf
    index <- with(orderedOrthoDf, order(orderNo))
    orderedOrthoDf <- orderedOrthoDf[index, ]

    #turn feature column into a character vector
    orderedOrthoDf$feature <- as.character(orderedOrthoDf$feature)
    #then turn it back into an ordered factor (to keep this order when plotting)
    orderedOrthoDf$feature <- factor(
        orderedOrthoDf$feature, levels = unique(orderedOrthoDf$feature)
    )
    #return sorted df
    return(orderedOrthoDf)
}
