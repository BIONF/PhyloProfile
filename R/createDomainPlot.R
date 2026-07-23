#' Create protein's domain architecure plot
#' @description Create architecture plot for both seed and orthologous protein.
#' If domains of ortholog are missing, only architecture of seed protein will
#' be plotted. NOTE: seed protein ID is the one being shown in the profile plot,
#' which normally is also the orthologous group ID.
#' @export
#' @usage createArchiPlot(info, domainDf, labelArchiSize, titleArchiSize,
#'     legendArchiSize, showScore, showWeight, namePosition, firstDist, 
#'     nameType, nameSize, segmentSize, nameColor, labelPos, colorType, 
#'     ignoreInstanceNo, currentNCBIinfo, featureClassSort, featureClassOrder, 
#'     colorPalette, resolveOverlap, font)
#' @param info A list contains seed and ortholog's IDs
#' @param domainDf Dataframe contains domain info for the seed and ortholog.
#' This including the seed ID, orthologs IDs, sequence lengths, feature names,
#' start and end positions, feature weights (optional) and the status to
#' determine if that feature is important for comparison the architecture
#' between 2 proteins* (e.g. seed protein vs ortholog) (optional).
#' @param labelArchiSize Lable size (in px). Default = 12.
#' @param titleArchiSize Title size (in px). Default = 12.
#' @param legendArchiSize Title size (in px). Default = 12.
#' @param showScore Show/hide E-values and Bit-scores. Default = NULL (hide)
#' @param showWeight Show/hide feature weights. Default = NULL (hide)
#' @param namePosition list of positions for domain names, choose from "plot",
#' "legend" or "axis". Default: "plot"
#' @param firstDist Distance of the first domain to plot title. Default = 0.5
#' @param nameType Type of domain names, either "Texts" or "Labels" (default)
#' @param nameSize Size of domain names. Default = 3
#' @param segmentSize Height of domain segment. Default = 5
#' @param nameColor Color of domain names (for Texts only). Default = "black"
#' @param labelPos Position of domain names (for Labels only). Choose from
#' @param colorType Choose to color "all", "shared", "unique" features or color
#' by "Feature type". Default = "all"
#' @param ignoreInstanceNo Ignore number of feature instances while identifying
#' shared or unique features. Default = FALSE
#' @param currentNCBIinfo Dataframe of the pre-processed NCBI taxonomy
#' data. Default = NULL (will be automatically retrieved from PhyloProfile app)
#' @param featureClassSort Choose to sort features. Default = "Yes"
#' @param featureClassOrder vector of ordered feature classes
#' @param colorPalette Choose between "Paired", "Set1", "Set2", "Set3",
#' "Accent", "Dark2" for the color pallete
#' @param resolveOverlap Choose to merge non-overlapped features of a feature
#' type into one line. Default = "Yes"
#' @param font font of text. Default = Arial"
#' @return A domain plot as arrangeGrob object. Use grid::grid.draw(plot) to
#' render.
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @importFrom stringr str_count
#' @seealso \code{\link{singleDomainPlotting}},\code{\link{pairDomainPlotting}},
#' \code{\link{sortDomains}}, \code{\link{parseDomainInput}}
#' @examples
#' seedID <- "101621at6656"
#' orthoID <- "101621at6656|AGRPL@224129@0|224129_0:001955|1"
#' info <- c(seedID, orthoID)
#' domainFile <- system.file(
#'     "extdata", "domainFiles/101621at6656.domains",
#'     package = "PhyloProfile", mustWork = TRUE
#' )
#' domainDf <- parseDomainInput(seedID, domainFile, "file")
#' domainDf$feature_id_mod <- domainDf$feature_id
#' domainDf$feature_id_mod <- gsub("SINGLE", "LCR", domainDf$feature_id_mod)
#' domainDf$feature_id_mod[domainDf$feature_type == "coils"] <- "Coils"
#' domainDf$feature_id_mod[domainDf$feature_type == "seg"] <- "LCR"
#' domainDf$feature_id_mod[domainDf$feature_type == "tmhmm"] <- "TM"
#' plot <- createArchiPlot(info, domainDf, font = "sans")
#' grid::grid.draw(plot)

createArchiPlot <- function(
    info = NULL, domainDf = NULL, labelArchiSize = 12, titleArchiSize = 12,
    legendArchiSize = 12, showScore = NULL, showWeight = NULL, 
    namePosition = "plot", firstDist = 0.5, nameType = "Labels", 
    nameSize = 3, segmentSize = 5,nameColor = "#000000", labelPos = "Above", 
    colorType = "Unique", ignoreInstanceNo = FALSE, currentNCBIinfo = NULL,
    featureClassSort = "Yes", featureClassOrder = NULL,
    colorPalette = "Paired", resolveOverlap = "Yes", font = "Arial"
){
    if (is.null(info) | is.null(domainDf)) return(ggplot() + theme_void())
    group <- as.character(info[1])
    ortho <- as.character(info[2])
    # get sub dataframe based on selected groupID and orthoID
    group <- gsub("\\|", ":", group)
    ortho <- gsub("\\|", ":", ortho)
    grepID <- paste(group, "#", ortho, sep = "")
    subdomainDf <- domainDf
    subdomainDf$feature <- as.character(subdomainDf$feature)
    subdomainDf <- subdomainDf[!duplicated(subdomainDf), ]
    orthoID <- NULL

    if (nrow(subdomainDf) < 1) return(paste0("No domain info available!"))
    else {
        # get minStart and maxEnd
        minStart <- min(subdomainDf$start)
        maxEnd <- max(subdomainDf$end)
        if ("length" %in% colnames(subdomainDf))
            maxEnd <- max(c(subdomainDf$end, subdomainDf$length))
        # ortho & seed domains df
        orthoDf <- subdomainDf[subdomainDf$orthoID == ortho,]
        seedDf <- subdomainDf[subdomainDf$orthoID != ortho,]
        if (nrow(seedDf) == 0) seedDf <- orthoDf
        seed <- as.character(seedDf$orthoID[1])
        if (nrow(seedDf) == 0) return(paste0("No domain info available!"))

        # simplify seed/ortho seq IDs if they are in bionf format
        if (!is.null(currentNCBIinfo)) {
            if (
                stringr::str_count(seed, ":") >= 2 &
                stringr::str_count(seed, "@") >= 2
            ) {
                seedTmp <- strsplit(as.character(seed),':', fixed = TRUE)[[1]]
                seedSpec <-
                    strsplit(as.character(seedTmp[2]),'@', fixed = TRUE)[[1]][2]
                seed <- paste0(
                    id2name(seedSpec, currentNCBIinfo)[,2], " - ", seedTmp[3]
                )
                if (ortho != seed) {
                    orthoTmp <- strsplit(
                        as.character(ortho),':', fixed = TRUE)[[1]]
                    orthoSpec <-
                        strsplit(
                            as.character(orthoTmp[2]),'@',fixed = TRUE)[[1]][2]
                    ortho <- paste0(
                        id2name(orthoSpec, currentNCBIinfo)[,2],
                        " - ", orthoTmp[3]
                    )
                }
            }
        }

        # add feature colors
        featureColorDf <- addFeatureColors(
            seedDf, orthoDf, colorType, colorPalette, ignoreInstanceNo
        )
        seedDf <- featureColorDf[[1]]
        orthoDf <- featureColorDf[[2]]

        # resolve (non)overlapped features
        if (resolveOverlap == "Yes") {
            seedDf <- resolveOverlapFeatures(seedDf)
            orthoDf <- resolveOverlapFeatures(orthoDf)
        } else {
            seedDf$featureOri <- seedDf$feature
            seedDf$featureOri <- as.character(seedDf$featureOri)
            orthoDf$featureOri <- orthoDf$feature
            orthoDf$featureOri <- as.character(orthoDf$featureOri)
        }

        if (nrow(orthoDf) > 0) {
            if (all.equal(seedDf, orthoDf)[1] == TRUE) featureClassSort <- "No"
            # sort features
            if (featureClassSort == "Yes") {
                # change order of one df's features based on order of other df's
                if (length(orthoDf$feature) < length(seedDf$feature)) {
                    orderedOrthoDf <- orthoDf[order(orthoDf$feature), ]
                    orderedSeedDf <- sortDomains(orderedOrthoDf, seedDf)
                    orderedOrthoDf <- sortDomains(orderedSeedDf, orderedOrthoDf)
                } else {
                    orderedSeedDf <- seedDf[order(seedDf$feature), ]
                    orderedOrthoDf <- sortDomains(orderedSeedDf, orthoDf)
                    orderedSeedDf <- sortDomains(orderedOrthoDf, orderedSeedDf)
                }
            } else {
                # change order based on list of feature types
                orderedSeedDf <- sortDomainsByList(seedDf, featureClassOrder)
                orderedOrthoDf <- sortDomainsByList(orthoDf, featureClassOrder)
            }
            # plotting
            g <- pairDomainPlotting(
                seed, ortho, orderedSeedDf, orderedOrthoDf, minStart, maxEnd,
                labelArchiSize, titleArchiSize, legendArchiSize, showScore, 
                showWeight, namePosition, firstDist, nameType, nameSize, 
                segmentSize, nameColor, labelPos, colorPalette, font)
        } else {
            orderedSeedDf <- sortDomainsByList(seedDf, featureClassOrder)
            # plotting
            g <- pairDomainPlotting(
                seed, seed, orderedSeedDf, orderedSeedDf, minStart, maxEnd,
                labelArchiSize, titleArchiSize, legendArchiSize, showScore, 
                showWeight, namePosition, firstDist, nameType, nameSize, 
                segmentSize, nameColor, labelPos, colorPalette, font)
        }
        return(g)
    }
}


#' Create architecure plot for a single protein
#' @usage singleDomainPlotting(df, geneID, sep, labelSize, titleSize, minStart,
#'     maxEnd, colorPalette, showScore, showWeight, namePosition, firstDist,
#'     nameType, nameSize, segmentSize, nameColor, labelPos, font)
#' @param df Domain dataframe for ploting containing the seed ID, ortholog ID,
#' ortholog sequence length, feature names, start and end positions,
#' feature weights (optional) and the status to determine if that feature is
#' important for comparison the architecture between 2 proteins* (e.g. seed
#' protein vs ortholog) (optional)
#' @param geneID ID of seed or orthologous protein
#' @param sep Separate indicator for title. Default = "|"
#' @param labelSize Lable size. Default = 12
#' @param titleSize Title size. Default = 12
#' @param minStart The smallest start position of all domains
#' @param maxEnd The highest stop position of all domains
#' @param colorPalette Color pallete. Default = Paired"
#' @param showScore Show/hide E-values and Bit-scores. Default = NULL (hide)
#' @param showWeight Show/hide feature weights. Default = NULL (hide)
#' @param namePosition List of positions for domain names, choose from "plot",
#' "legend" or "axis". Default: "plot"
#' @param firstDist Distance of the first domain to plot title. Default = 0.5
#' @param nameType Type of domain names, either "Texts" or "Labels" (default)
#' @param nameSize Size of domain names. Default = 3
#' @param segmentSize Height of domain segment. Default = 5
#' @param nameColor Color of domain names (for Texts only). Default = "black"
#' @param labelPos Position of domain names (for Labels only). Choose from
#' "Above" (default), "Below" or "Inside" the domain bar
#' @param font font of text. Default = Arial"
#' @return Domain plot of a single protein as a ggplot object.
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stringr str_wrap
#' @seealso \code{\link{parseDomainInput}}
#' @import ggplot2
#' @examples
#' seed <- "101621at6656"
#' ortho <- "101621at6656|AGRPL@224129@0|224129_0:001955|1"
#' ortho <- gsub("\\|", ":", ortho)
#' grepID <- paste(seed, "#", ortho, sep = "")
#' domainFile <- system.file(
#'     "extdata", "domainFiles/101621at6656.domains",
#'     package = "PhyloProfile", mustWork = TRUE
#' )
#' domainDf <- parseDomainInput(seed, domainFile, "file")
#' domainDf$feature_id_mod <- domainDf$feature_id
#' subdomainDf <- domainDf[grep(grepID, domainDf$seedID), ]
#' subdomainDf$feature <- as.character(subdomainDf$feature)
#' orthoDf <- subdomainDf[subdomainDf$orthoID == ortho,]
#' seedDf <- subdomainDf[subdomainDf$orthoID != ortho,]
#' minStart <- min(subdomainDf$start)
#' maxEnd <- max(c(subdomainDf$end, subdomainDf$length))
#' # resolve overlapping domains
#' seedDf <- PhyloProfile:::resolveOverlapFeatures(seedDf)
#' orthoDf <- PhyloProfile:::resolveOverlapFeatures(orthoDf)
#' # add feature colors
#' featureColorDf <- PhyloProfile:::addFeatureColors(seedDf, orthoDf)
#' seedDf <- featureColorDf[[1]]
#' orthoDf <- featureColorDf[[2]]
#' # do plot
#' g <- PhyloProfile:::singleDomainPlotting(
#'     seedDf, seed, minStart = minStart, maxEnd = maxEnd, font = "sans"
#' )
#' grid::grid.draw(g)

singleDomainPlotting <- function(
        df = NULL, geneID = "GeneID", sep = "|", labelSize = 12, titleSize = 12,
        minStart = NULL, maxEnd = NULL, colorPalette = "Set2",
        showScore = NULL, showWeight = NULL, namePosition = "plot",
        firstDist = 0.5, nameType = "Labels", nameSize = 3, segmentSize = 5,
        nameColor = "#000000", labelPos = "Above", font = "Arial"
){
    feature <- feature_id_mod <- end <- start <- featureOri <- NULL
    evalue <- bitscore <- weight <- NULL

    # parse parameters
    if (is.null(df)) return(ggplot() + theme_void())
    if (is.null(minStart)) minStart <- min(df$start)
    if (is.null(maxEnd)) maxEnd <- max(df$end)
    if ("color" %in% colnames(df)) {
        colorScheme <- structure(
            df$color, .Names = df$featureOri
        )
    } else {
        colorScheme <- structure(
            head(
                suppressWarnings(
                    RColorBrewer::brewer.pal(
                        nlevels(as.factor(df$feature)), colorPalette
                    )
                ),
                levels(as.factor(df$feature))
            ),
            .Names = levels(as.factor(df$feature)))
    }
    # initiate ggplot object
    gg <- ggplot(df, aes(y = feature, x = end))

    # draw lines for representing sequence length
    if ("length" %in% colnames(df)) {
        gg <- gg + geom_segment(
            data = df, linewidth = 0.6, color = "#b2b2b2", alpha = 0.5,
            aes(x = 0, xend = length, y = feature, yend = feature))
    }

    # draw features
    gg <- gg + geom_segment(
        data = df, aes(
            x = start, xend = end, y = feature, yend = feature,
            color = as.factor(featureOri)
        ), linewidth = segmentSize, alpha = 1.0) +
        scale_color_manual(values = colorScheme)

    # add feature names
    if ("plot" %in% namePosition) {
        if (nameType == "Labels") {
            if (labelPos == "Above") {
                gg <- gg + geom_label(
                    aes(
                        label = stringr::str_wrap(feature_id_mod),
                        x = (start+end)/2
                    ),
                    color = "black", vjust = -0.5, size = nameSize, family=font
                )
            } else if (labelPos == "Below") {
                gg <- gg + geom_label(
                    aes(
                        label = stringr::str_wrap(feature_id_mod),
                        x = (start+end)/2
                    ),
                    color = "black", vjust = 1.5, size = nameSize, family = font
                )
            } else {
                gg <- gg + geom_label(
                    aes(
                        label = stringr::str_wrap(feature_id_mod),
                        x = (start+end)/2
                    ),
                    color = "black", size = nameSize, family = font
                )
            }
        } else {
            gg <- gg + geom_text(
                aes(
                    label = stringr::str_wrap(feature_id_mod),
                    x = (start+end)/2
                ),
                color = nameColor, check_overlap = TRUE, size = nameSize,
                family = font
            )
        }
    }
    # add scores if selected
    maxVjust <- 2
    if ("Bit-score" %in% showScore & "E-value" %in% showScore) {
        gg <- gg + geom_text(
            aes(
                label = ifelse(
                    evalue == "NA" & bitscore == "NA", "", stringr::str_wrap(
                        paste0("E: ", evalue, "; B: ", bitscore)
                    )
                ), x = (start+end)/2
            ), color = "black", vjust = maxVjust, check_overlap = TRUE
        )
    } else if ("Bit-score" %in% showScore) {
        gg <- gg + geom_text(
            aes(
                label = ifelse(bitscore == "NA", "", stringr::str_wrap(
                    paste0("B: ", bitscore)
                )), x = (start+end)/2
            ), color = "black", vjust = maxVjust, check_overlap = TRUE
        )
    } else if ("E-value" %in% showScore) {
        gg <- gg + geom_text(
            aes(
                label = ifelse(evalue == "NA", "", stringr::str_wrap(
                    paste0("E: ", evalue)
                )), x = (start+end)/2
            ), color = "black", vjust = maxVjust, check_overlap = TRUE
        )
    }
    # add weight if selected
    if ("Weight" %in% showWeight) {
        if ("E-value" %in% showScore |  "Bit-score" %in% showScore) {
            maxVjust <- 3.5
        }
        gg <- gg + geom_text(
            aes(
                label = ifelse(weight == "NA", "", stringr::str_wrap(
                    paste0("W: ", weight)
                )), x = (start+end)/2
            ), color = "black", vjust = maxVjust, check_overlap = TRUE
        )
    }
    # theme format
    gg <- gg + labs(
        title = paste0(gsub(":", sep, geneID)), color = "Feature"
    )
    gg <- gg + theme(
        axis.line = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.minor.x=element_blank(), panel.grid.major.x=element_blank(),
        panel.background = element_blank()
    )
    gg <- gg + theme(
        plot.title = element_text(face = "bold", size = titleSize, hjust = 0.5)
    )
    gg <- gg + theme(
        axis.title = element_blank(),
        axis.text.x = element_text(size = labelSize))
    gg <- gg + theme(text = element_text(family = font))
    # add feature names on the axis (if required)
    if ("axis" %in% namePosition) {
        if ("plot" %in% namePosition | "legend" %in% namePosition) {
            gg <- gg + scale_y_discrete(
                expand = c(0.075, 0), breaks = df$feature,
                labels = df$feature_type)
        } else {
            gg <- gg + scale_y_discrete(
                expand = c(0.075, 0), breaks = df$feature, labels = df$feature)
        }
        gg <- gg + theme(axis.text.y = element_text(size = labelSize))
    } else {
        gg <- gg + theme(axis.text.y = element_blank())
    }
    # add legend (if required)
    if ("legend" %in% namePosition) {
        gg <- gg + theme(legend.position = "bottom")
    } else {
        gg <- gg + theme(legend.position = "none")
    }
    # add space on the top of the plot (for feature name)
    gg <- gg + coord_cartesian(
        clip = "off", ylim = c(1, nlevels(as.factor(df$feature)) + firstDist)
    )
    # scale x-axis to the length of the longest protein
    gg <- gg + xlim(0, maxEnd)
    return(gg)
}


#' Create architecure plot for a pair of seed and ortholog protein
#' @usage pairDomainPlotting(seed, ortho, seedDf, orthoDf, minStart, maxEnd,
#'     labelSize, titleSize, legendSize, showScore, showWeight, namePosition, 
#'     firstDist, nameType, nameSize, segmentSize, nameColor, labelPos, 
#'     colorPalette, font)
#' @param seed Seed ID
#' @param ortho Ortho ID
#' @param seedDf domain dataframe for seed domains containing the seed ID,
#' ortholog ID, sequence length, feature names, start and end positions,
#' feature weights (optional) and the status to determine if that feature is
#' important for comparison the architecture between 2 proteins* (e.g. seed
#' protein vs ortholog) (optional)
#' @param orthoDf domain dataframe for ortholog domains (same format as seedDf)
#' @param minStart the smallest start position of all domains
#' @param maxEnd the highest stop position of all domains
#' @param labelSize lable size. Default = 12
#' @param titleSize title size. Default = 12
#' @param legendSize legend size. Default = 12
#' @param showScore show/hide E-values and Bit-scores. Default = NULL (hide)
#' @param showWeight Show/hide feature weights. Default = NULL (hide)
#' @param namePosition list of positions for domain names, choose from "plot",
#' "legend" or "axis". Default: "plot"
#' @param firstDist distance of the first domain to plot title. Default = 0.5
#' @param nameType type of domain names, either "Texts" or "Labels" (default)
#' @param nameSize Size of domain names. Default = 3
#' @param segmentSize Height of domain segment. Default = 5
#' @param nameColor color of domain names (for Texts only). Default = "black"
#' @param labelPos position of domain names (for Labels only). Choose from
#' "Above" (default), "Below" or "Inside" the domain bar
#' @param colorPalette color pallete. Default = Paired"
#' @param font font of text. Default = Arial"
#' @return Domain plot of a pair proteins as a arrangeGrob object.
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @importFrom gridExtra grid.arrange
#' @seealso \code{\link{singleDomainPlotting}}, \code{\link{sortDomains}},
#' \code{\link{parseDomainInput}}
#' @examples
#' seed <- "101621at6656"
#' ortho <- "101621at6656|AGRPL@224129@0|224129_0:001955|1"
#' ortho <- gsub("\\|", ":", ortho)
#' grepID <- paste(seed, "#", ortho, sep = "")
#' domainFile <- system.file(
#'     "extdata", "domainFiles/101621at6656.domains",
#'     package = "PhyloProfile", mustWork = TRUE
#' )
#' domainDf <- parseDomainInput(seed, domainFile, "file")
#' domainDf$feature_id_mod <- domainDf$feature_id
#' subdomainDf <- domainDf[grep(grepID, domainDf$seedID), ]
#' subdomainDf$feature <- as.character(subdomainDf$feature)
#' orthoDf <- subdomainDf[subdomainDf$orthoID == ortho,]
#' seedDf <- subdomainDf[subdomainDf$orthoID != ortho,]
#' minStart <- min(subdomainDf$start)
#' maxEnd <- max(c(subdomainDf$end, subdomainDf$length))
#' # resolve overlapping domains
#' seedDf <- PhyloProfile:::resolveOverlapFeatures(seedDf)
#' orthoDf <- PhyloProfile:::resolveOverlapFeatures(orthoDf)
#' # add feature colors
#' featureColorDf <- PhyloProfile:::addFeatureColors(seedDf, orthoDf)
#' seedDf <- featureColorDf[[1]]
#' orthoDf <- featureColorDf[[2]]
#' # do plot
#' g <- PhyloProfile:::pairDomainPlotting(
#'     seed,ortho,seedDf,orthoDf,minStart,maxEnd, font = "sans"
#' )
#' grid::grid.draw(g)

pairDomainPlotting <- function(
        seed = NULL, ortho = NULL, seedDf = NULL, orthoDf = NULL,
        minStart = 0, maxEnd = 999, labelSize = 12, titleSize = 12,
        legendSize = 12, showScore = NULL, showWeight = NULL, 
        namePosition = "plot", firstDist = 0.5, nameType = "Labels", 
        nameSize = 3, segmentSize = 5, nameColor = "#000000",  
        labelPos = "Above", colorPalette = "Paired", font = "Arial"
) {
    if(is.null(seed) | is.null(ortho) | is.null(seedDf) | is.null(orthoDf))
        stop("Seed/Ortho ID or domain dataframe is NULL!")

    sep <- "|"
    plotSeed <- singleDomainPlotting(
        seedDf, seed, sep, labelSize, titleSize, minStart, maxEnd, colorPalette,
        showScore, showWeight, namePosition, firstDist, nameType, nameSize,
        segmentSize, nameColor, labelPos, font
    )
    if (ortho == seed) {
        g <- plotSeed
    } else {
        plotOrtho <- singleDomainPlotting(
            orthoDf, ortho, sep, labelSize, titleSize, minStart, maxEnd,
            colorPalette, showScore, showWeight, namePosition, firstDist,
            nameType, nameSize, segmentSize, nameColor, labelPos, font
        )
        if ("legend" %in% namePosition) {
            g <- joinPlotMergeLegends(
                seedDf, orthoDf, plotSeed, plotOrtho, position = "bottom", font,
                legendSize
            )
        } else {
            seedHeight <- length(levels(as.factor(seedDf$feature)))
            orthoHeight <- length(levels(as.factor(orthoDf$feature)))
            g <- gridExtra::grid.arrange(
                plotSeed, plotOrtho, nrow = 2,
                heights = c(seedHeight, orthoHeight)
            )
        }
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
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @examples
#' # get domain data
#' seedID <- "101621at6656"
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

    # convert feature column into a character vector
    orderedOrthoDf$feature <- as.character(orderedOrthoDf$feature)
    # then convert it back into an ordered factor
    # (to keep this order when plotting)
    orderedOrthoDf$feature <- factor(
        orderedOrthoDf$feature, levels = unique(orderedOrthoDf$feature)
    )
    # return sorted df
    return(orderedOrthoDf)
}

#' Sort one domain dataframe based on list of ordered feature types
#' @description Sort domain dataframe of one protein based on a given list of
#' ordered feature types
#' @param domainDf domain dataframe
#' @param featureClassOrder vector of ordered feature classes
#' @return Dataframe contains sorted domain list.
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @importFrom dplyr left_join
#' @examples
#' # get domain data
#' seedID <- "101621at6656"
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
#' featureClassOrder <- c("pfam", "smart", "tmhmm", "coils", "signalp", "seg",
#'     "flps")
#' # sort
#' PhyloProfile:::sortDomainsByList(orthoDf, featureClassOrder)

sortDomainsByList <- function(domainDf = NULL, featureClassOrder = NULL) {
    if (is.null(domainDf) | is.null(featureClassOrder))
        stop("Domain data or feature type order is NULL!")
    featureClassOrder <- rev(featureClassOrder[
        featureClassOrder %in% levels(as.factor(domainDf$feature_type))
    ])
    orderedDomainDf <- dplyr::left_join(
        data.frame(feature_type = featureClassOrder),
        domainDf, by = "feature_type"
    )
    orderedDomainDf$feature <- factor(
        orderedDomainDf$feature, levels = unique(orderedDomainDf$feature)
    )
    return(orderedDomainDf)
}

#' Modify feature names
#' @description Simplify feature names (e.g. TM for transmembrane domain,
#' LCR for low complexity regions, remove tool names from domain name) and add
#' weight to feature names (if available)
#' @param domainDf domain data as a dataframe object
#' @return Dataframe contains simlified domain names in yLabel column
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @examples
#' domainFile <- system.file(
#'     "extdata", "domainFiles/101621at6656.domains",
#'     package = "PhyloProfile", mustWork = TRUE
#' )
#' seedID <- "101621at6656"
#' domainDf <- parseDomainInput(seedID, domainFile, "file")
#' PhyloProfile:::modifyFeatureName(domainDf)

modifyFeatureName <- function(domainDf = NULL) {
    if (is.null(domainDf)) stop("Domain data cannot be NULL!")
    domainDf$yLabel <- domainDf$feature_id
    domainDf$yLabel[domainDf$yLabel == "transmembrane"] <- "TM"
    domainDf$yLabel[domainDf$yLabel == "low complexity regions"] <- "LCR"
    domainDf$yLabel[domainDf$yLabel == "low_complexity_regions"] <- "LCR"
    return(domainDf)
}


#' Join multiple plots and merge legends
#' @usage joinPlotMergeLegends(df1, df2, plot1, plot2, position, font, 
#'     legendSize)
#' @param df1 Data frame for plot 1
#' @param df2 Data frame for plot 2
#' @param plot1 ggplot object of plot 1
#' @param plot2 ggplot object of plot 2
#' @param position position of legend (bottom or right)
#' @param font font of text
#' @param legendSize font size
#' @return joined plots with merged legend as a grid object
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @importFrom gridExtra grid.arrange
#' @importFrom grid unit grobHeight convertHeight
#' @examples
#' seed <- "101621at6656"
#' ortho <- "101621at6656|AGRPL@224129@0|224129_0:001955|1"
#' ortho <- gsub("\\|", ":", ortho)
#' grepID <- paste(seed, "#", ortho, sep = "")
#' domainFile <- system.file(
#'     "extdata", "domainFiles/101621at6656.domains",
#'     package = "PhyloProfile", mustWork = TRUE
#' )
#' domainDf <- parseDomainInput(seed, domainFile, "file")
#' domainDf$feature_id_mod <- domainDf$feature_id
#' subdomainDf <- domainDf[grep(grepID, domainDf$seedID), ]
#' subdomainDf$feature <- as.character(subdomainDf$feature)
#' orthoDf <- subdomainDf[subdomainDf$orthoID == ortho,]
#' seedDf <- subdomainDf[subdomainDf$orthoID != ortho,]
#' minStart <- min(subdomainDf$start)
#' maxEnd <- max(c(subdomainDf$end, subdomainDf$length))
#' # resolve overlapping domains
#' seedDf <- PhyloProfile:::resolveOverlapFeatures(seedDf)
#' orthoDf <- PhyloProfile:::resolveOverlapFeatures(orthoDf)
#' # add feature colors
#' featureColorDf <- PhyloProfile:::addFeatureColors(seedDf, orthoDf)
#' seedDf <- featureColorDf[[1]]
#' orthoDf <- featureColorDf[[2]]
#' # generate plots
#' plotSeed <- PhyloProfile:::singleDomainPlotting(
#'     seedDf, seed, minStart = minStart, maxEnd = maxEnd, font = "sans"
#' )
#' plotOrtho <- PhyloProfile:::singleDomainPlotting(
#'     orthoDf, ortho, minStart = minStart, maxEnd = maxEnd, font = "sans"
#' )
#' # merge plots
#' PhyloProfile:::joinPlotMergeLegends(
#'     seedDf, orthoDf, plotSeed, plotOrtho, "bottom", font = "sans", 
#'     legendSize = 12)

joinPlotMergeLegends <- function(
    df1 = NULL, df2 = NULL, plot1 = NULL, plot2 = NULL, 
    position = c("bottom", "right"), font = "Arial", legendSize = 12)
{
    if (is.null(plot1) | is.null(df1)) stop("No plot data given!")
    if (is.null(plot2) | is.null(df2))
        return(plot1 + theme(legend.position = position))
    featureOri <- NULL
    # remove legend of the original plots
    plot1 <- plot1 + theme(legend.position = "none")
    plot2 <- plot2 + theme(legend.position = "none")
    # create a temp plot that contains all features
    mergedDf <- rbind(df1, df2)
    colorScheme <- structure(
        mergedDf$color, .Names = mergedDf$featureOri
    )
    tmpPlot <- ggplot() +
        geom_segment(
            data = mergedDf,
            aes(
                x = 0, xend = 1, y = 1, yend = 1,
                color = featureOri
            ),
            linewidth = 5, lineend = "round", alpha = 0.7
        ) +
        scale_color_manual(values = colorScheme) +
        labs(color = "Feature") +
        theme_minimal() +
        theme(
            legend.position = position, 
            text = element_text(family = font, size = legendSize)
        )
    # extract legend from the temp plot above
    getOnlyLegend <- function(plot) {
        plotTable <- ggplot_gtable(ggplot_build(plot))
        legendPlot <- which(
            vapply(
                plotTable$grobs, function(x) x$name, character(1)
            ) == "guide-box"
        )
        legend <- plotTable$grobs[[legendPlot]]
        return(legend)
    }
    legend <- getOnlyLegend(tmpPlot)
    legend_height <- grid::convertHeight(
        grid::grobHeight(legend), "cm", valueOnly = TRUE
    )
    # Combine all three vertically, using measured heights
    combined <- gridExtra::grid.arrange(
        plot1,
        plot2,
        legend,
        ncol = 1, nrow = 3,
        heights = grid::unit.c(
            grid::unit(length(unique(df1$feature)), "null"),
            grid::unit(length(unique(df2$feature)), "null"),
            grid::unit(legend_height, "cm")
        )
    )
    return(combined)
}

#' Identify feature type(s) containing overlapped domains/features
#' @param domainDf input domain dataframe
#' @return List of feature types that have overlapped domains
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @importFrom dplyr group_by add_count
#' @examples
#' # get domain data
#' seedID <- "101621at6656"
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
#' # check overlap features
#' PhyloProfile:::checkOverlapDomains(orthoDf)

checkOverlapDomains <- function(domainDf) {
    if (is.null(domainDf)) stop("Domain data cannot be NULL!")
    feature_type <- NULL
    df <- domainDf[, c("feature", "feature_type", "start", "end")]
    df <- data.frame(
        df %>% dplyr::group_by(feature_type) %>% dplyr::add_count(feature_type)
    )
    df$tmp <- paste(df$feature_type, df$end, sep = "_")
    overlappedType <- lapply(
        df$tmp[df$n > 1],
        function (x) {
            if (length(df$feature[df$tmp == x]) > 1)
                return(df$feature_type[df$tmp == x])

            ed <- df$end[df$tmp == x][1]
            st <- df$start[df$tmp == x][1]
            type <- df$feature_type[df$tmp == x][1]

            subDf <- df[
                !(df$tmp == x) & df$feature_type == type &
                    df$start < ed & df$end > st,
            ]
            if (nrow(subDf) > 0) {
                if(!(subDf$feature[1] == df$feature[df$tmp == x]))
                    return(subDf$feature_type)
            }
        }
    )
    return(unique(unlist(overlappedType)))
}

#' Modify domain dataframe to resolve overlapped domains/features
#' @param domainDf input domain dataframe
#' @return Domain dataframe with modified feature names that join multiple
#' domains of the same type that are not overlapped
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @examples
#' # get domain data
#' seedID <- "101621at6656"
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
#' # resolve overlapped featuers
#' PhyloProfile:::resolveOverlapFeatures(orthoDf)

resolveOverlapFeatures <- function(domainDf) {
    if (is.null(domainDf)) stop("Domain data cannot be NULL!")
    overlappedType <- checkOverlapDomains(domainDf)
    domainDf$featureOri <- domainDf$feature
    domainDf$featureOri <- as.character(domainDf$featureOri)
    if (length(overlappedType) > 0) {
        domainDf$feature <- ifelse(
            !(domainDf$feature_type %in% overlappedType),
            domainDf$feature_type,
            domainDf$feature
        )
    } else {
        domainDf$feature <- domainDf$feature_type
    }
    return(domainDf)
}

#' Linearize PFAM/SMART annotations by best e-value/bitscore
#' @export
#' @param domainDf input domain dataframe
#' @param orthoID ID of protein that needs to be linearized
#' @param value type of values that will be used for linearized, either evalue
#' (default) or bitscore
#' @return Domain dataframe of the selected protein after linearization
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @importFrom dplyr arrange
#' @examples
#' demoDomainDf <- data.frame(
#'     orthoID = rep("protID", 4),
#'     start = c(1, 5, 100, 80),
#'     end = c(30, 40, 130, 110),
#'     evalue = c(0.001, 0.0005, 0.2, 0.004),
#'     feature_type = c(rep("pfam", 2), rep("smart", 2)),
#'     feature_id = c("pf1", "pf2", "sm1", "sm2")
#' )
#' linearizeArchitecture(demoDomainDf, "protID", "evalue")

linearizeArchitecture <- function(
        domainDf = NULL, orthoID = NULL, value = "evalue"
) {
    if (is.null(domainDf) | is.null(orthoID)) stop("Input data is NULL!")
    evalue <- bitscore <- start <- end <- NULL
    # Sort the dataframe by start position and then by evalue or bitscore
    if (value == "evalue") {
        domainDf <- domainDf %>% dplyr::arrange(start, evalue)
    } else if (value == "bitscore") {
        domainDf <- domainDf %>% dplyr::arrange(start, bitscore)
    } else stop("Incorrect value specified! Either 'evalue' or 'bitscore'")

    # Get lines that need to be excluded
    pfamRows <- rownames(domainDf[domainDf$orthoID == orthoID &
                    domainDf$feature_type %in% c("pfam","smart"),])
    exclude_lines <- vapply(
        seq_len(length(pfamRows)-1),
        function(i) {
            if (domainDf[pfamRows[i],]$end >= domainDf[pfamRows[i+1],]$start) {
                # Exclude the row with the higher evalue / lower bitscore
                if (value == "evalue") {
                    if (
                        domainDf[pfamRows[i],]$evalue >
                        domainDf[pfamRows[i+1],]$evalue
                    ) {
                        return((pfamRows[i]))
                    } else {
                        return((pfamRows[i+1]))
                    }
                } else {
                    if (
                        domainDf[pfamRows[i],]$bitscore <
                        domainDf[pfamRows[i+1],]$bitscore
                    ) {
                        return((pfamRows[i]))
                    } else {
                        return((pfamRows[i+1]))
                    }
                }

            } else {
                return("0")
            }
        },
        character(1)
    )
    # return domainDf after removing overlapped features with higher e-values
    outDf <- domainDf[!(row.names(domainDf) %in% exclude_lines), ]
    return(outDf[outDf$orthoID == orthoID,])
}

#' Add colors for each feature/domain
#' @usage addFeatureColors(seedDf, orthoDf, colorType, colorPalette,
#'     ignoreInstanceNo)
#' @description Add colors to features/domains of 2 domain dataframes. Users can
#' choose to color only the shared features, unique features, all features
#' (default) or based on feature types. Default color pallete is "Paired", but
#' it can be changed.
#' @param seedDf Domain dataframe of seed protein (protein 1)
#' @param orthoDf Domain dataframe of orthologs protein (protein 2)
#' @param colorType Choose to color "all", "shared", "unique" features or color
#' by "Feature type". Default: "all"
#' @param colorPalette Choose between "Paired", "Set1", "Set2", "Set3",
#' "Accent", "Dark2" for the color pallete
#' @param ignoreInstanceNo Ignore number of feature instances while identifying
#' shared or unique features. Default: FALSE
#' @return 2 dataframes (seedDf and orthoDf) with an additional column for the
#' assigned color to each feature instance
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @importFrom RColorBrewer brewer.pal
#' @importFrom utils head
#' @importFrom dplyr count
#' @examples
#' # get domain data
#' seedID <- "101621at6656"
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
#' # add colors to features
#' PhyloProfile:::addFeatureColors(seedDf, orthoDf)

addFeatureColors <- function(
        seedDf = NULL, orthoDf = NULL, colorType = "all",
        colorPalette = "Paired", ignoreInstanceNo = FALSE
) {
    if (is.null(seedDf) | is.null(orthoDf)) stop("Domain Df cannot be null!")
    feature <- NULL
    featureSeedCount <- seedDf %>% dplyr::count(feature)
    featureOrthoCount <- orthoDf %>% dplyr::count(feature)
    featureCount <- merge(
        featureSeedCount, featureOrthoCount, by = "feature", all = TRUE
    )
    featureCount$type <- "unique"
    featureCount$type[featureCount$n.x == featureCount$n.y] <- "shared"
    if (ignoreInstanceNo == TRUE) {
        featureCount$type[
            !is.na(featureCount$n.x) & !is.na(featureCount$n.y)] <- "shared"
    }
    featureCount$feature <- as.character(featureCount$feature)
    sharedFeatures <-unique(featureCount$feature[featureCount$type == "shared"])
    uniqueFeatures <-unique(featureCount$feature[featureCount$type == "unique"])
    allFeatures <- c(sharedFeatures, uniqueFeatures)

    if (colorType == "Unique" & length(uniqueFeatures) == 0)
        colorType <- "All"
    if (colorType == "Shared" & length(sharedFeatures) == 0)
        colorType <- "All"

    if (colorType == "Unique") {
        sharedFeaturesColors <- rep("#C9C9C9", length(sharedFeatures))
        uniqueFeaturesColors <- getQualColForVector(uniqueFeatures)
        if (checkColorPalette(uniqueFeatures, colorPalette) == TRUE) {
            uniqueFeaturesColors <-
                suppressWarnings(head(
                    RColorBrewer::brewer.pal(
                        length(uniqueFeatures), colorPalette),
                    length(uniqueFeatures)
                ))
        }
        allColors <- c(sharedFeaturesColors, uniqueFeaturesColors)
        colorScheme <- data.frame(color = allColors, feature = allFeatures)
    } else if (colorType == "Shared") {
        sharedFeaturesColors <- getQualColForVector(sharedFeatures)
        uniqueFeaturesColors <- rep("#C9C9C9", length(uniqueFeatures))
        if (checkColorPalette(sharedFeatures, colorPalette) == TRUE) {
            sharedFeaturesColors <-
                suppressWarnings(head(
                    RColorBrewer::brewer.pal(
                        length(sharedFeatures), colorPalette),
                    length(sharedFeatures)
                ))
        }
        allColors <- c(sharedFeaturesColors, uniqueFeaturesColors)
        colorScheme <- data.frame(color = allColors, feature = allFeatures)
    } else if (colorType == "All") {
        allColors <- getQualColForVector(allFeatures)
        if (checkColorPalette(allFeatures, colorPalette) == TRUE) {
            allColors <-
                suppressWarnings(head(
                    RColorBrewer::brewer.pal(length(allFeatures), colorPalette),
                    length(allFeatures)
                ))
        }
        colorScheme <- data.frame(color = allColors, feature = allFeatures)
    } else {
        tmpDf <- data.frame(str_split_fixed(allFeatures, "_", 2))
        tmpDf$name <- paste(tmpDf$X1, tmpDf$X2, sep = "_")
        tmpDf$name[tmpDf$X2 == ""] <- tmpDf$X1[tmpDf$X2 == ""]
        tmpDf$X2[tmpDf$X2 == ""] <- tmpDf$X1[tmpDf$X2 == ""]
        typeColorDf <- data.frame(
            colors = head(
                suppressWarnings(RColorBrewer::brewer.pal(
                    nlevels(as.factor(tmpDf$X1)), colorPalette)),
                nlevels(as.factor(tmpDf$X1))
            ),
            X1 = levels(as.factor(tmpDf$X1))
        )
        tmpDf <- merge(tmpDf, typeColorDf, all.x = TRUE)
        colorScheme <- data.frame(color = tmpDf$colors, feature = tmpDf$name)
    }

    # add color to seedDf and orthoDf
    seedDf <- merge(seedDf, colorScheme, by = "feature", all.x = TRUE)
    orthoDf <- merge(orthoDf, colorScheme, by = "feature", all.x = TRUE)

    return(list(seedDf, orthoDf))
}
