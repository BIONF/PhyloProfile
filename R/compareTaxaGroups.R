#' Get all taxa that share a common ancestor
#' @description Identify the common ancestor for a selected taxa and return a 
#' list of all taxa that have that common ancestor from an large input taxa set.
#' @export
#' @param inputTaxa ID list of all input taxa (e.g. "ncbi12345")
#' @param inGroup ID list of selected taxa used for identify the common ancestor
#' (e.g.: "ncbi55555")
#' @return A list containing the taxonomy rank and name of the common ancestor,
#' together with a dataframe storing the full taxonomy info of all taxa that
#' share that corresponding common ancestor.
#' @author Vinh Tran (tran@bio.uni-frankfurt.de)
#' @examples
#' inputTaxa <- c("ncbi103449", "ncbi1288291", "ncbi278021", "ncbi27973",
#'     "ncbi31281", "ncbi40302", "ncbi586133", "ncbi6035", "ncbi70536",
#'     "ncbi876142", "ncbi993615")
#' inGroup <-  c("ncbi27973", "ncbi6035")
#' getCommonAncestor(inputTaxa, inGroup)

getCommonAncestor <- function(inputTaxa = NULL, inGroup = NULL) {
    if (is.null(inputTaxa) | is.null(inGroup)) return()
    
    # get list of pre-calculated taxonomy info
    taxMatrix <- getTaxonomyMatrix(TRUE, inputTaxa)
    # get subset taxonomy info for selected in-group taxa
    selectedTaxMatrix <- taxMatrix[
        taxMatrix$abbrName %in% inGroup, 
        which(!duplicated(t(taxMatrix)))
    ]
    # identify common ancestor
    # identify common ancestor
    V1 <- vapply(
        selectedTaxMatrix[,c(seq(4, ncol(selectedTaxMatrix)))], max,
        FUN.VALUE = numeric(1)
    )
    V2 <- vapply(
        selectedTaxMatrix[,c(seq(4, ncol(selectedTaxMatrix)))],
        function (x) as.integer(sum(x)/length(x)),
        FUN.VALUE = numeric(1)
    )
    checkDf <- t(rbind(V1, V2))
    commonRank <- rownames(checkDf[V1 == V2,])[1]
    commonID <- checkDf[V1 == V2,][1]
    commonTaxa <- taxMatrix[taxMatrix[, commonRank] == commonID, ]
    
    return(list(commonRank, commonID, commonTaxa))
}

#' Compare the score distributions between 2 taxon groups
#' @description Given the phylogenetic profiles that contains up to 2 additional
#' variables besides the presence/absence information of the orthologous 
#' proteins. This function will compare the distribution of those variables 
#' between 2 different taxon groups (e.g. parasitic species vs non-parasitic 
#' species), which are defined as in-group and out-group. In-group is identified
#' by the user. Out-group contains all taxa in the input phylogenetic profiles
#' that are not part of the in-group.
#' @usage compareTaxonGroups(data, inGroup, useCommonAncestor, variable, 
#'     significanceLevel)
#' @export
#' @param data input phylogenetic profile in long format (see ?mainLongRaw and 
#' ?createLongMatrix)
#' @param inGroup ID list of in-group taxa (e.g. "ncbi1234")
#' @param useCommonAncestor TRUE/FALSE if using all taxa that share the same 
#' common ancestor with the pre-selected in-group as the in-group taxa. 
#' Default = TRUE.
#' @param variable name of the variable that need to be compared
#' @param significanceLevel significant cutoff for the statistic test (between 
#' 0 and 1). Default = 0.05.
#' @return list of genes that have a significant difference in the variable 
#' distributions between the in-group and out-group taxa and their corresponding
#' p-values.
#' @author Vinh Tran (tran@bio.uni-frankfurt.de)
#' @examples
#' data("mainLongRaw", package="PhyloProfile")
#' data <- mainLongRaw
#' inGroup <- c("ncbi876142", "ncbi586133")
#' variable <- colnames(data)[4]
#' compareTaxonGroups(data, inGroup, TRUE, variable, 0.05)

compareTaxonGroups <- function(
    data = NULL, 
    inGroup = NULL, 
    useCommonAncestor = TRUE, 
    variable = NULL, 
    significanceLevel = 0.05
) {
    if (is.null(data) | is.null(inGroup) | is.null(variable)) return()
    if (!(variable %in% colnames(data))) {
        message("Invalid variable")
        return()
    }
    
    # add other taxa that share a common ancestor with the given in-group
    commonTaxa <- getCommonAncestor(levels(as.factor(data$ncbiID)), inGroup)
    if (useCommonAncestor == TRUE) {
        inGroup <- as.character(commonTaxa[[3]]$abbrName)
    }
    
    # perform distribution comparison test
    data$geneID <- as.character(data$geneID)
    pvalues <- vapply(
        unique(data$geneID),
        function (x) {
            varIn <- data[data$geneID == x & data$ncbiID %in% inGroup, variable]
            varOut <- data[
                data$geneID == x & !(data$ncbiID %in% inGroup), variable
                ]
            pvalue <- distributionTest(varIn, varOut, significanceLevel)
            if (!is.null(pvalue))
                return(distributionTest(varIn, varOut, significanceLevel))
            else return(999)
        },
        FUN.VALUE = numeric(1)
    )
    return(sort(unlist(pvalues)))
}

#' Compare the distribution of 2 numeric vectors
#' @description This function tests the difference between the distributions of 
#' two input numeric samples using the statistical tess. First the 
#' Kolmogorov-Smirnov is used to check if 2 samples have the same distribution.
#' If yes, Wilcoxon-Mann-Whitney will be used to compare the distribution 
#' difference.
#' @usage distributionTest(varIn, varOut, significanceLevel)
#' @param varIn first numeric vector
#' @param varOut second numeric vector
#' @param significanceLevel significant cutoff of the Kolmogorov-Smirnov test. 
#' Default = 0.05.
#' @return p-value of the comparison test.
#' @importFrom stats ks.test
#' @importFrom stats wilcox.test
#' @author Carla Mölbert (carla.moelbert@gmx.de)

distributionTest <- function(
    varIn = NULL, varOut = NULL, significanceLevel = 0.05
){
    if (is.null(varIn) | is.null(varOut)) return()
    # remove NA values
    varIn <- varIn[!is.na(varIn)]
    varOut <- varOut[!is.na(varOut)]
    
    # if there is no data in one of the groups the p-value is NULL
    if (length(varIn) == 0 | length(varOut) == 0) return()
    else {
        # * Kolmogorov-Smirnov Test
        # H0 : The two samples have the same distribution
        ks <- suppressWarnings(
            ks.test(unique(varIn), unique(varOut), exact = FALSE)
        )
        
        if (ks$p.value <= significanceLevel) return(ks$p.value)
        else {
            # * Wilcoxon-Mann-Whitney Test
            # H0: the samples have the same location parameters
            wilcox <- suppressWarnings(
                wilcox.test(
                    varIn, varOut, alternative = "two.sided", paired = FALSE
                )
            )
            return(wilcox$p.value)
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
        # return(perm$p.value)
    }
}

#' Compare the mean values of a variable between 2 taxon groups
#' @description Given the phylogenetic profiles that contains up to 2 additional
#' variables besides the presence/absence information of the orthologous 
#' proteins. This function will compare the mean scores of those variables 
#' between 2 different taxon groups (e.g. parasitic species vs non-parasitic 
#' species), which are defined as in-group and out-group. In-group is identified
#' by the user. Out-group contains all taxa in the input phylogenetic profiles
#' that are not part of the in-group.
#' @usage compareMeansTaxonGroups(data, inGroup, useCommonAncestor, variable)
#' @export
#' @param data input phylogenetic profile in long format (see ?mainLongRaw and 
#' ?createLongMatrix)
#' @param inGroup ID list of in-group taxa (e.g. "ncbi1234")
#' @param useCommonAncestor TRUE/FALSE if using all taxa that share the same 
#' common ancestor with the pre-selected in-group as the in-group taxa. 
#' Default = TRUE.
#' @param variable name of the variable that need to be compared
#' @return List of genes that have a difference in the variable's mean scores
#' between the in-group and out-group taxa and their corresponding delta-mean.
#' @author Vinh Tran (tran@bio.uni-frankfurt.de)
#' @examples
#' data("mainLongRaw", package="PhyloProfile")
#' data <- mainLongRaw
#' inGroup <- c("ncbi876142", "ncbi586133")
#' variable <- colnames(data)[4]
#' compareMeansTaxonGroups(data, inGroup, TRUE, variable)

compareMeansTaxonGroups <- function(
    data = NULL, inGroup = NULL, useCommonAncestor = TRUE, variable = NULL
) {
    if (is.null(data) | is.null(inGroup) | is.null(variable)) return()
    if (!(variable %in% colnames(data))) {
        message("Invalid variable")
        return()
    }
    
    # add other taxa that share a common ancestor with the given in-group
    commonTaxa <- getCommonAncestor(levels(as.factor(data$ncbiID)), inGroup)
    if (useCommonAncestor == TRUE) {
        inGroup <- as.character(commonTaxa[[3]]$abbrName)
    }
    
    # return delta-mean scores for two taxa groups
    data$geneID <- as.character(data$geneID)
    deltaMean <- vapply(
        unique(data$geneID),
        function (x) {
            varIn <- data[data$geneID == x & data$ncbiID %in% inGroup, variable]
            varOut <- data[
                data$geneID == x & !(data$ncbiID %in% inGroup), variable]
            return(
                abs(mean(varIn[!is.na(varIn)]) - mean(varOut[!is.na(varOut)]))
            )
        },
        FUN.VALUE = numeric(1)
    )
    return(sort(unlist(deltaMean)))
}


#' Create data for variable distribution comparison plot
#' @description Create data for plotting the distribution comparison between 2 
#' groups of taxa for a selected gene.
#' @usage dataVarDistTaxGroup(data, inGroup, gene, variable)
#' @export
#' @param data input phylogenetic profile in long format (see ?mainLongRaw and 
#' ?createLongMatrix)
#' @param inGroup ID list of in-group taxa (e.g. "ncbi1234")
#' @param gene ID of gene that need to be plotted the distribution comparison 
#' between in- and out-group taxa.
#' @param variable var1 or c(var1, var2)
#' @return Dataframe containing list of values for all available variables for 
#' the selected genes in in-group and out-group taxa (max. 3 columns).
#' @author Vinh Tran (tran@bio.uni-frankfurt.de)
#' @seealso \code{\link{createLongMatrix}}
#' @examples
#' data("mainLongRaw", package="PhyloProfile")
#' data <- mainLongRaw
#' inGroup <- c("ncbi876142", "ncbi586133")
#' variable <- colnames(data)[c(4, 5)]
#' dataVarDistTaxGroup(data, inGroup, "OG_1017", variable)

dataVarDistTaxGroup <- function(
    data = NULL,
    inGroup = NULL,
    gene = NULL,
    variable = NULL
) {
    if (is.null(data) | is.null(inGroup) | is.null(gene) | is.null(variable)) 
        return()
    # remove "empty" variable (char "")
    variable <- variable[unlist(lapply(variable, function (x) x != ""))]
    # get 2 lists of values for in-group and out-group
    varIn <- data.frame(
        data[data$geneID == gene & data$ncbiID %in% inGroup, ][, variable]
    )
    colnames(varIn) <- variable
    varOut <- data.frame(
        data[data$geneID == gene & !(data$ncbiID %in% inGroup), ][, variable]
    )
    colnames(varOut) <- variable
    if (nrow(varIn) == 0 & nrow(varOut) == 0) return()
    
    varIn$type <- "In-group"
    varOut$type <- "Out-group"
    out <- rbind(varIn[complete.cases(varIn),], varOut[complete.cases(varOut),])
    return(out[, c(variable, "type")])
}

#' Create variable distribution comparison plot 
#' @description Create variable distribution plots between 2 groups of taxa for
#' a selected gene.
#' @export
#' @param data dataframe for plotting. Last column 
#' indicates what type of taxon group (in- or out-group). The first (or first 2)
#' column contains values of the variables. See ?dataVarDistTaxGroup
#' @param plotParameters plot parameters, including size of x-axis, y-axis, 
#' legend and title; position of legend ("right", "bottom" or "none"); 
#' mean/median point; names of in-group and out-group; and plot title.
#' NOTE: Leave blank or NULL to use default values.
#' @return Distribution plots as a grob (gtable) object. Use grid.draw to plot.
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @seealso \code{\link{dataVarDistTaxGroup}}
#' @importFrom ggplot2 geom_violin
#' @importFrom ggplot2 position_dodge
#' @importFrom ggplot2 geom_boxplot
#' @importFrom ggplot2 scale_x_discrete
#' @importFrom ggplot2 stat_summary
#' @importFrom ggplot2 ggplotGrob
#' @importFrom grid unit.c
#' @importFrom grid unit
#' @importFrom grid textGrob
#' @importFrom grid gpar
#' @examples
#' data("mainLongRaw", package="PhyloProfile")
#' data <- mainLongRaw
#' inGroup <- c("ncbi876142", "ncbi586133")
#' variable <- colnames(data)[c(4, 5)]
#' plotDf <- dataVarDistTaxGroup(data, inGroup, "OG_1017", variable)
#' plotParameters <- list(
#'     "xSize" = 12,
#'     "ySize" = 12,
#'     "titleSize" = 15,
#'     "legendSize" = 12,
#'     "legendPosition" = "right",
#'     "mValue" = "mean",
#'     "inGroupName" = "In-group",
#'     "outGroupName" = "Out-group",
#'     "title" = "OG_1017"
#' )
#' g <- varDistTaxPlot(plotDf, plotParameters)
#' grid::grid.draw(g)

varDistTaxPlot <- function(data, plotParameters) {
    if (is.null(data)) return()
    if (missing(plotParameters)) return()
    
    # rename in-group and out-group
    data$type[data$type == "In-group"] <- plotParameters$inGroupName
    data$type[data$type == "Out-group"] <- plotParameters$outGroupName
    
    # function for plotting a single plot
    generatePlot <- function(plotDf, parameters, variable) {
        type <- NULL
        .data <- NULL
        xNames <- c(
            paste(
                parameters$inGroupName, " \n n = ", 
                nrow(plotDf[plotDf$type == parameters$inGroupName,]), sep = ""
            ),
            paste(
                parameters$outGroupName, " \n n = ", 
                nrow(plotDf[plotDf$type == parameters$outGroupName,]), sep = ""
            )
        )
        
        plot <- ggplot(plotDf, aes(x = factor(type), y = .data[[variable]])) +
            geom_violin(
                aes(fill = factor(type)), position = position_dodge(), 
                scale = "width", alpha = .5) +
            geom_boxplot(width = 0.1) +
            scale_x_discrete(labels = xNames) +
            # add mean/median point and color of that point
            stat_summary(
                aes(colour = parameters$mValue), fun.y = parameters$mValue,
                geom = "point", size = 3, show.legend = TRUE, shape = 8
            ) +
            scale_color_manual("", values = c("red")) +
            # add title and theme
            labs(x = element_blank(), y = variable) +
            theme_minimal() +
            theme(
                axis.text.x = element_text(size = parameters$xSize, hjust = 1),
                axis.text.y = element_text(size = parameters$ySize),
                axis.title.y = element_text(size = parameters$ySize),
                legend.position = parameters$legendPosition,
                legend.text = element_text(size = parameters$legendSize),
                legend.title = element_blank()
            )
        return(plot)
    }
    
    # adapted from http://rpubs.com/sjackman/grid_arrange_shared_legend
    gridArrangeSharedLegend <- function(
        ...,  ncol = length(list(...)), nrow = 1, 
        position = c("bottom", "right"), title = NA
    ) {
        plots <- list(...)
        position <- match.arg(position)
        g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
        legend <- g[[
            which(vapply(g, function(x) x$name, FUN.VALUE = "") == "guide-box")
            ]]
        lheight <- sum(legend$height)
        lwidth <- sum(legend$width)
        gl <- lapply(plots, function(x) x + theme(legend.position="none"))
        gl <- c(gl, ncol = ncol, nrow = nrow)
        
        combined <- switch(
            position,
            "bottom" = arrangeGrob(
                do.call(arrangeGrob, gl),
                legend,
                ncol = 1,
                heights = unit.c(unit(1, "npc") - lheight, lheight),
                top = textGrob(
                    title, 
                    vjust = 1, 
                    gp = gpar(
                        fontface = "bold", fontsize = plotParameters$titleSize
                    )
                )
            ),
            "right" = arrangeGrob(
                do.call(arrangeGrob, gl),
                legend,
                ncol = 2,
                widths = unit.c(unit(1, "npc") - lwidth, lwidth),
                top = textGrob(
                    title, 
                    vjust = 1, 
                    gp = gpar(
                        fontface = "bold", fontsize = plotParameters$titleSize
                    )
                )
            )
        )
        return(combined)
    }
    
    # return plot(s)
    if (ncol(data) == 2) {
        plotVar1 <- generatePlot(data, plotParameters, colnames(data)[1])
        return(
            arrangeGrob(
                plotVar1,
                top = textGrob(
                    plotParameters$title, vjust = 1, 
                    gp = gpar(
                        fontface = "bold", fontsize = plotParameters$titleSize
                    )
                )
            )
        )
    } else {
        plotVar1 <- generatePlot(data, plotParameters, colnames(data)[1])
        plotVar2 <- generatePlot(data, plotParameters, colnames(data)[2])
        if (plotParameters$legendPosition == "none") {
            return(
                arrangeGrob(
                    plotVar1, plotVar2,
                    nrow = 1,
                    top = textGrob(
                        plotParameters$title, vjust = 1, 
                        gp = gpar(
                            fontface = "bold", 
                            fontsize = plotParameters$titleSize
                        )
                    )
                )
            )
        } else {
            return(
                gridArrangeSharedLegend(
                    plotVar1, plotVar2, 
                    position = plotParameters$legendPosition, 
                    title = plotParameters$title
                )
            )
        }
    }
}

#' Create data for feature distribution comparison plot
#' @description Create data for plotting the distribution of the protein domain 
#' features between 2 group of taxa for a selected gene.
#' @usage dataFeatureTaxGroup(mainDf, domainDf, inGroup, gene, featureThreshold)
#' @export
#' @param mainDf input phylogenetic profile in long format (see ?mainLongRaw 
#' and ?createLongMatrix)
#' @param domainDf dataframe contains domain info for the seed and ortholog. 
#' This including the seed ID, orthologs IDs, sequence lengths, feature names, 
#' start and end positions, feature weights (optional) and the status to 
#' determine if that feature is important for comparison the architecture 
#' between 2 proteins* (e.g. seed protein vs ortholog) (optional). (see 
#' ?parseDomainInput)
#' @param inGroup ID list of in-group taxa (e.g. "ncbi1234")
#' @param gene ID of gene that need to be plotted the feature distribution 
#' comparison between in- and out-group taxa.
#' @param featureThreshold cutoff for occurence frequencies of the protein 
#' feature instances in each taxon group (in percentage). Default = 0.
#' @return Dataframe containing all feature names, their frequencies (absolute 
#' count and percentage) in each taxon group and the corresponding taxa group.
#' @author Vinh Tran (tran@bio.uni-frankfurt.de)
#' @seealso \code{\link{createLongMatrix}}, \code{\link{parseDomainInput}}
#' @examples
#' data("mainLongRaw", package="PhyloProfile")
#' data <- mainLongRaw
#' gene <- "OG_1017"
#' inputFile <- system.file(
#'     "extdata", "domainFiles/OG_1017.domains",
#'     package = "PhyloProfile", mustWork = TRUE
#' )
#' type <- "file"
#' domainDf <- parseDomainInput(gene, inputFile, type)
#' inGroup <- c("ncbi876142", "ncbi586133")
#' featureThreshold <- 15
#' dataFeatureTaxGroup(data, domainDf, inGroup, gene, featureThreshold)

dataFeatureTaxGroup <- function(
    mainDf = NULL,
    domainDf = NULL,
    inGroup = NULL,
    gene = NULL,
    featureThreshold = 0
) {
    if (is.null(mainDf) | is.null(inGroup) | is.null(gene) | is.null(domainDf)) 
        return()
    
    mainDf$orthoID <- gsub("\\|", ":", mainDf$orthoID)
    domainDfSub <- merge(
        domainDf[grep(gene, domainDf$seedID),],
        mainDf[, c("orthoID", "ncbiID")],
        by = "orthoID", all.x = TRUE
    )
    domainDfSub <- domainDfSub[complete.cases(domainDfSub),]
    
    domainIn <- domainDfSub$feature[domainDfSub$ncbiID %in% inGroup]
    domainOut <- domainDfSub$feature[!(domainDfSub$ncbiID %in% inGroup)]
    
    countDomainIn <- data.frame(table(unlist(as.character(domainIn))))
    countDomainIn$type <- "In-group"
    countDomainIn$percentage <- (countDomainIn$Freq/sum(countDomainIn$Freq))*100
    
    countDomainOut <- data.frame(table(unlist(as.character(domainOut))))
    countDomainOut$type <- "Out-group"
    countDomainOut$percentage <- 
        (countDomainOut$Freq/sum(countDomainOut$Freq))*100
    
    outDf <- rbind(countDomainIn, countDomainOut)
    colnames(outDf) <- c("Feature", "Count", "Taxon_group", "Percentage")
    return(outDf[outDf$Percentage >= featureThreshold,])
}

#' Create feature distribution comparison plot
#' @description Create protein feature distribution plots between 2 groups of 
#' taxa for a selected gene.
#' @export
#' @param data dataframe for plotting (see ?dataFeatureTaxGroup)
#' @param plotParameters plot parameters, including size of x-axis, y-axis, 
#' legend and title; position of legend ("right", "bottom" or "none"); names of
#' in-group and out-group; flip the plot coordinate ("Yes" or "No").
#' NOTE: Leave blank or NULL to use default values.
#' @return Distribution plots as a ggplot2 object.
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @seealso \code{\link{dataFeatureTaxGroup}}
#' @examples
#' data("mainLongRaw", package="PhyloProfile")
#' data <- mainLongRaw
#' gene <- "OG_1017"
#' inputFile <- system.file(
#'     "extdata", "domainFiles/OG_1017.domains",
#'     package = "PhyloProfile", mustWork = TRUE
#' )
#' type <- "file"
#' domainDf <- parseDomainInput(gene, inputFile, type)
#' inGroup <- c("ncbi876142", "ncbi586133")
#' plotDf <- dataFeatureTaxGroup(data, domainDf, inGroup, gene, 15)
#' plotParameters <- list(
#'     "xSize" = 12,
#'     "ySize" = 12,
#'     "angle" = 15,
#'     "legendSize" = 12,
#'     "inGroupName" = "In-group",
#'     "outGroupName" = "Out-group",
#'     "flipPlot" = "No"
#' )
#' featureDistTaxPlot(plotDf, plotParameters)

featureDistTaxPlot <- function(data, plotParameters) {
    Feature <- NULL
    Percentage <- NULL
    Taxon_group <- NULL
    if (is.null(data)) return()
    if (missing(plotParameters)) return()
    
    data$Taxon_group[data$Taxon_group == "In-group"] <- 
        plotParameters$inGroupName
    data$Taxon_group[data$Taxon_group == "Out-group"] <- 
        plotParameters$outGroupName
    
    plot <- ggplot(data, aes(x = Feature, y = Percentage, fill = Taxon_group)) +
        geom_bar(stat="identity", width=.5, position = "dodge") +
        theme_minimal() +
        theme(
            axis.title.y = element_text(size = plotParameters$ySize),
            axis.text.y = element_text(size = plotParameters$ySize),
            axis.title.x = element_text(size = plotParameters$xSize),
            axis.text.x = element_text(
                size = plotParameters$xSize, 
                angle = plotParameters$angle, 
                hjust = 1
            ),
            legend.text = element_text(size = plotParameters$legendSize),
            legend.title = element_blank()
        )
    if (plotParameters$flipPlot == "Yes")
        plot <- plot + coord_flip()
    
    return(plot)
}



#### WHAT IS p.adjust() ????? ################
# significantGenesDf$pvalues <- {
#     p.adjust(
#         significantGenesDf$pvalues, method = "holm",
#         n = length(significantGenesDf$pvalues)
#     )
# }