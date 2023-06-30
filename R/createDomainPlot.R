#' Create protein's domain architecure plot
#' @description Create architecture plot for both seed and orthologous protein.
#' If domains of ortholog are missing, only architecture of seed protein will
#' be plotted. NOTE: seed protein ID is the one being shown in the profile plot,
#' which normally is also the orthologous group ID.
#' @export
#' @usage createArchiPlot(info = NULL, domainDf = NULL, labelArchiSize = 12,
#'     titleArchiSize = 12, showFeature = "all", seqIdFormat = "unknown")
#' @param info a list contains seed and ortholog's IDs
#' @param domainDf dataframe contains domain info for the seed and ortholog.
#' This including the seed ID, orthologs IDs, sequence lengths, feature names,
#' start and end positions, feature weights (optional) and the status to
#' determine if that feature is important for comparison the architecture
#' between 2 proteins* (e.g. seed protein vs ortholog) (optional).
#' @param labelArchiSize lable size (in px). Default = 12.
#' @param titleArchiSize title size (in px). Default = 12.
#' @param showFeature choose to show all, common or unique features. 
#' Default = "all"
#' @param seqIdFormat sequence ID format (either bionf or unknown). 
#' Default = "unknown"
#' @importFrom gridExtra arrangeGrob
#' @import ggplot2
#' @return A domain plot as arrangeGrob object. Use grid::grid.draw(plot) to
#' render.
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @seealso \code{\link{singleDomainPlotting}}, \code{\link{sortDomains}},
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
    info = NULL, domainDf = NULL, labelArchiSize = 12, titleArchiSize = 12,
    showFeature = "all", seqIdFormat = "unknown"
){
    if (is.null(info) | is.null(domainDf)) return(ggplot() + theme_void())
    group <- as.character(info[1])
    ortho <- as.character(info[2])
    # get sub dataframe based on selected groupID and orthoID
    group <- gsub("\\|", ":", group)
    ortho <- gsub("\\|", ":", ortho)
    grepID <- paste(group, "#", ortho, sep = "")
    subdomainDf <- domainDf[grep(grepID, domainDf$seedID), ]
    subdomainDf$feature <- as.character(subdomainDf$feature)
    orthoID <- NULL
    if (nrow(subdomainDf) < 1) return(paste0("No domain info available!"))
    else {
        # ortho & seed domains df
        orthoDf <- subdomainDf[subdomainDf$orthoID == ortho,]
        seedDf <- subdomainDf[subdomainDf$orthoID != ortho,]
        # filter common features
        if (!(showFeature == "all")) {
            allFeats <- c(
                levels(as.factor(orthoDf$feature)), 
                levels(as.factor(seedDf$feature))
            )
            countFeats <- as.data.frame(table(allFeats))
            commonFeats <- countFeats$allFeats[countFeats$Freq > 1]
            if (showFeature == "common") {
                orthoDf <- orthoDf[orthoDf$feature %in% commonFeats,]
                seedDf <- seedDf[seedDf$feature %in% commonFeats,]
            } else {
                orthoDf <- orthoDf[!(orthoDf$feature %in% commonFeats),]
                seedDf <- seedDf[!(seedDf$feature %in% commonFeats),]
            }
        }
        # final check
        if (nrow(seedDf) == 0) seedDf <- orthoDf
        seed <- as.character(seedDf$orthoID[1])
        if (nrow(seedDf) == 0) return(paste0("No domain info available!"))
        if (nrow(orthoDf) == 0) {
            ortho <- seed
            orthoDf <- seedDf
        }
        # change order of one df's features based on order of other df's
        if (length(orthoDf$feature) < length(seedDf$feature)) {
            if (nrow(orthoDf) > 0) {
                orderedOrthoDf <- orthoDf[order(orthoDf$feature), ]
                orderedSeedDf <- sortDomains(orderedOrthoDf, seedDf)
            }
        } else {
            orderedSeedDf <- seedDf[order(seedDf$feature), ]
            orderedOrthoDf <- sortDomains(orderedSeedDf, orthoDf)
        }
        # join weight values and feature names
        if ("weight" %in% colnames(orderedOrthoDf)) {
            orderedOrthoDf$yLabel <- paste0(
                orderedOrthoDf$feature," (",round(orderedOrthoDf$weight, 2),")")
        } else orderedOrthoDf$yLabel <- orderedOrthoDf$feature
        if ("weight" %in% colnames(orderedSeedDf)) {
            orderedSeedDf$yLabel <- paste0(
                orderedSeedDf$feature," (",round(orderedSeedDf$weight, 2),")")
        } else orderedSeedDf$yLabel <- orderedSeedDf$feature
        # simplify seq IDs if they are in bionf format
        if (seqIdFormat == "bionf") {
            seedTmp <- strsplit(as.character(seed),':', fixed = TRUE)[[1]]
            seed <- paste0(seedTmp[2], " - ", seedTmp[3])
            orthoTmp <- strsplit(as.character(ortho),':', fixed = TRUE)[[1]]
            ortho <- paste0(orthoTmp[2], " - ", orthoTmp[3])
        }
        # plotting
        minStart <- min(subdomainDf$start)
        maxEnd <- max(subdomainDf$end)
        if ("length" %in% colnames(subdomainDf))
            maxEnd <- max(c(subdomainDf$end, subdomainDf$length))
        g <- pairDomainPlotting(
            seed, ortho, orderedSeedDf, orderedOrthoDf, minStart, maxEnd,
            labelArchiSize, titleArchiSize)
        return(g)
    }
}

#' Create architecure plot for a single protein
#' @usage singleDomainPlotting(df, geneID = "GeneID", sep = "|", labelSize = 12,
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
#' PhyloProfile:::singleDomainPlotting(
#'     df,
#'     geneID,
#'     sep,
#'     labelSize, titleSize,
#'     minStart, maxEnd,
#'     colorScheme
#' )
#' }

singleDomainPlotting <- function(
    df = NULL, geneID = "GeneID", sep = "|", labelSize = 12, titleSize = 12,
    minStart = NULL, maxEnd = NULL, colorScheme = NULL
){
    feature <- end <- start <- NULL
    # parse parameters
    if (is.null(df)) return(ggplot() + theme_void())
    if (is.null(minStart)) minStart <- min(df$start)
    if (is.null(maxEnd)) maxEnd <- max(df$end)
    if (is.null(colorScheme)) {
        colorScheme <- structure(
            getQualColForVector(levels(as.factor(df$feature))),
            .Names = levels(as.factor(df$feature)))}
    gg <- ggplot(df, aes(y = feature, x = end, color = as.factor(feature))) +
        geom_segment(
            data = df, color = "white", linewidth = 0,
            aes(y = feature, yend = feature, x = minStart, xend = maxEnd)) +
        scale_color_manual(values = colorScheme)
    # draw lines for representing sequence length
    if ("length" %in% colnames(df))
        gg <- gg + geom_segment(
            data = df, linewidth = 1, color = "#b2b2b2",
            aes(x = 0, xend = length, y = feature, yend = feature))
    # draw line and points
    gg <- gg + geom_segment(
        data = df, aes(x = start, xend = end, y = feature, yend = feature),
        linewidth = 1.5) +
        geom_point(data = df, aes(y = feature, x = start),
                    color = "#b2b2b2", size = 3, shape = 3) +
        geom_point(data = df, aes(y = feature, x = end),
                    color = "#edae52", size = 3, shape = 5)
    # draw dashed line for domain path
    gg <- gg + geom_segment(
        data = df[df$path == "Y", ], linewidth = 3, linetype = "dashed",
        aes(x = start, xend = end, y = feature, yend = feature))
    # theme format
    gg <- gg + scale_y_discrete(
        expand = c(0.075, 0), breaks = df$feature, labels = df$yLabel)
    gg <- gg + labs(title = paste0(gsub(":", sep, geneID)), y = "Feature")
    gg <- gg + theme_minimal() + theme(panel.border = element_blank())
    gg <- gg + theme(axis.ticks = element_blank())
    gg <- gg + theme(plot.title = element_text(face = "bold", size = titleSize))
    gg <- gg + theme(plot.title = element_text(hjust = 0.5))
    gg <- gg + theme(
        legend.position = "none", axis.title.x = element_blank(),
        axis.text.y = element_text(size = labelSize),
        axis.title.y = element_text(size = labelSize),
        panel.grid.minor.x=element_blank(), panel.grid.major.x=element_blank())
    return(gg)
}

#' Create architecure plot for a pair of seed and ortholog protein
#' @usage pairDomainPlotting(seed, ortho, seedDf, orthoDf, minStart, maxEnd,
#'     labelSize, titleSize)
#' @param seed Seed ID
#' @param ortho Ortho ID
#' @param seedDf domain dataframe for seed domains containing the seed ID,
#' ortholog ID, sequence length, feature names, start and end positions,
#' feature weights (optional) and the status to determine if that feature is
#' important for comparison the architecture between 2 proteins* (e.g. seed
#' protein vs ortholog) (optional).
#' @param orthoDf domain dataframe for ortholog domains (same format as seedDf).
#' @param minStart the smallest start position of all domains
#' @param maxEnd the highest stop position of all domains
#' @param labelSize lable size. Default = 12.
#' @param titleSize title size. Default = 12.
#' @return Domain plot of a pair proteins as a arrangeGrob object.
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @examples
#' \dontrun{
#' seed <- "101621at6656"
#' ortho <- "101621at6656|AGRPL@224129@0|224129_0:001955|1"
#' ortho <- gsub("\\|", ":", ortho)
#' grepID <- paste(seed, "#", ortho, sep = "")
#' domainFile <- system.file(
#'     "extdata", "domainFiles/101621at6656.domains",
#'     package = "PhyloProfile", mustWork = TRUE
#' )
#' domainDf <- parseDomainInput(seed, domainFile, "file")
#' subdomainDf <- domainDf[grep(grepID, domainDf$seedID), ]
#' subdomainDf$feature <- as.character(subdomainDf$feature)
#' orthoDf <- subdomainDf[subdomainDf$orthoID == ortho,]
#' seedDf <- subdomainDf[subdomainDf$orthoID != ortho,]
#' minStart <- min(subdomainDf$start)
#' maxEnd <- max(c(subdomainDf$end, subdomainDf$length))
#' g <- pairDomainPlotting(seed,ortho,seedDf,orthoDf,minStart,maxEnd,9,9)
#' grid::grid.draw(g)
#' }

pairDomainPlotting <- function(
    seed = NULL, ortho = NULL, seedDf = NULL, orthoDf = NULL,
    minStart = 0, maxEnd = 999, labelSize = 12, titleSize = 12
) {
    if(is.null(seed) | is.null(ortho) | is.null(seedDf) | is.null(orthoDf))
        stop("Seed/Ortho ID or domain dataframe is NULL!")
    # create color scheme, so that the same features in seed & ortholog will
    # have the same colors
    featureSeed <- levels(as.factor(seedDf$feature))
    featureOrtho <- levels(as.factor(orthoDf$feature))
    allFeatures <- c(featureSeed, featureOrtho)
    allColors <- getQualColForVector(allFeatures)
    colorScheme <- structure(allColors, .Names = allFeatures)
    # plot
    sep <- "|"
    plotOrtho <- singleDomainPlotting(
        orthoDf, ortho, sep, labelSize, titleSize, minStart, maxEnd,colorScheme)
    plotSeed <- singleDomainPlotting(
        seedDf, seed, sep, labelSize, titleSize, minStart, maxEnd, colorScheme)
    if (ortho == seed) {
        g <- gridExtra::arrangeGrob(plotSeed, ncol = 1)
    } else {
        seedHeight <- length(levels(as.factor(seedDf$feature)))
        orthoHeight <- length(levels(as.factor(orthoDf$feature)))
        g <- gridExtra::arrangeGrob(
            plotSeed, plotOrtho, ncol = 1, heights = c(seedHeight, orthoHeight)
        )
    }
    return(g)
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
#' PhyloProfile:::sortDomains(seedDf, orthoDf)
#' }

sortDomains <- function(seedDf, orthoDf){
    if (is.null(seedDf) | is.null(orthoDf))
        stop("Domain data for seed & ortholog cannot be NULL!")
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

