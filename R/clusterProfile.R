#' Get data for calculating the distance matrix
#' @export
#' @usage getDataClustering(data, profileType, var1AggregateBy,
#'     var2AggregateBy)
#' @param data processed profile data
#' @param profileType type of data used for calculating the distance matrix.
#' Either "binary" (consider only the presence/absence status of orthlogs), or
#' "var1"/"var2" for taking values of the additional variables into account.
#' @param var1AggregateBy aggregate method for VAR1 (min, max, mean or median)
#' @param var2AggregateBy aggregate method for VAR2 (min, max, mean or median)
#' @return A wide dataframe contains values for calculating distance matrix.
#' @author Carla Mölbert (carla.moelbert@gmx.de)
#' @note Documented by Vinh Tran (tran@bio.uni-frankfurt.de)
#' @seealso \code{\link{fromInputToProfile}}
#' @examples
#' data("fullProcessedProfile", package="PhyloProfile")
#' data <- fullProcessedProfile
#' profileType <- "binary"
#' var1AggregateBy <- "max"
#' var2AggregateBy <- "mean"
#' getDataClustering(data, profileType, var1AggregateBy, var2AggregateBy)

getDataClustering <- function(
    data, profileType, var1AggregateBy, var2AggregateBy
) {
    supertaxon <- NULL
    presSpec <- NULL

    # remove lines where there is no found ortholog
    subDataHeat <- subset(data, data$presSpec > 0)

    # transform data into wide matrix
    if (profileType == "binary") {
        subDataHeat <- subDataHeat[, c("geneID", "supertaxon", "presSpec")]
        subDataHeat$presSpec[subDataHeat$presSpec > 0] <- 1
        subDataHeat <- subDataHeat[!duplicated(subDataHeat), ]
        wideData <- spread(subDataHeat, supertaxon, presSpec)
    } else {
        var <- profileType

        subDataHeat <- subDataHeat[, c("geneID", "supertaxon", var)]
        subDataHeat <- subDataHeat[!duplicated(subDataHeat), ]

        # aggreagte the values by the selected method
        if (var == "var1") aggregateBy <- var1AggregateBy
        else aggregateBy <- var2AggregateBy

        subDataHeat <- aggregate(
            subDataHeat[, var],
            list(subDataHeat$geneID, subDataHeat$supertaxon),
            FUN = aggregateBy
        )

        colnames(subDataHeat) <- c("geneID", "supertaxon", var)
        wideData <- spread(subDataHeat, supertaxon, var)
    }

    # set name for wide matrix as gene IDs
    dat <- wideData[, 2:ncol(wideData)]
    rownames(dat) <- wideData[, 1]
    dat[is.na(dat)] <- 0

    return(dat)
}

#' Calculate the distance matrix
#' @export
#' @param profiles profile data for distance calculating
#' @param method distance calculation method ("euclidean", "maximum",
#' "manhattan", "canberra", "binary", "distanceCorrelation",
#' "mutualInformation" or "pearson" for binary data; "distanceCorrelation" or
#' "mutualInformation" for non-binary data)
#' @return A distance matrix for input phylogenetic profiles.
#' @author Carla Mölbert (carla.moelbert@gmx.de)
#' @note Documented by Vinh Tran (tran@bio.uni-frankfurt.de)
#' @seealso \code{\link{getDataClustering}}
#' @examples
#' data("fullProcessedProfileLarge", package="PhyloProfile")
#' data <- fullProcessedProfileLarge
#' profileType <- "binary"
#' profiles <- getDataClustering(
#'     data, profileType, var1AggregateBy, var2AggregateBy)
#' method <- "mutualInformation"
#' getDistanceMatrix(profiles, method)

getDistanceMatrix <- function(profiles, method) {

    profiles <-  profiles[, colSums(profiles != 0) > 0]
    profiles <-  profiles[rowSums(profiles != 0) > 0, ]

    distMethods <- c("euclidean", "maximum", "manhattan", "canberra", "binary")
    if (method %in% distMethods) {
        distanceMatrix <- dist(profiles, method = method)
    } else if (method == "distanceCorrelation") {
        matrix <- data.frame()
        for (i in seq_len(nrow(profiles))) { # rows
            for (j in seq_len(nrow(profiles))) { # columns
                if (i == j) {
                    matrix[i,i] = 1 # if this cell NA as.dist not work probably
                    break
                }
                dist <- dcor(unlist(profiles[i,]), unlist(profiles[j,]))
                # Swich the value so that the profiles with a high correlation
                # are clustered together
                matrix[i,j] <- 1 - dist
            }
        }

        profileNames <- rownames(profiles)
        colnames(matrix) <- profileNames[seq_len(length(profileNames)) - 1]
        rownames(matrix) <- profileNames
        distanceMatrix <- as.dist(matrix)
    } else if (method == "mutualInformation") {
        distanceMatrix <- mutualInfo(as.matrix(profiles))
        distanceMatrix <- max(distanceMatrix, na.rm = TRUE) - distanceMatrix
    } else if (method == "pearson") {
        distanceMatrix <-  cor.dist(as.matrix(profiles))
    }

    return(distanceMatrix)
}

#' Create a dendrogram tree from the distance matrix
#' @export
#' @param distanceMatrix calculated distance matrix
#' @param clusterMethod clustering method ("single", "complete",
#' "average" for UPGMA, "mcquitty" for WPGMA, "median" for WPGMC,
#' or "centroid" for UPGMC)
#' @return A dendrogram tree object
#' @importFrom stats as.dendrogram
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @seealso \code{\link{getDataClustering}},
#' \code{\link{getDistanceMatrix}}, \code{\link{hclust}}
#' @examples
#' data("fullProcessedProfileLarge", package="PhyloProfile")
#' data <- fullProcessedProfileLarge
#' profileType <- "binary"
#' profiles <- getDataClustering(
#'     data, profileType, var1AggregateBy, var2AggregateBy)
#' distMethod <- "mutualInformation"
#' distanceMatrix <- getDistanceMatrix(profiles, distMethod)
#' clusterMethod <- "complete"
#' clusterDataDend(distanceMatrix, clusterMethod)

clusterDataDend <- function(distanceMatrix, clusterMethod) {
    if (is.null(distanceMatrix)) return()
    dd.col <- as.dendrogram(hclust(distanceMatrix, method = clusterMethod))
    return(dd.col)
}

#' Plot dendrogram tree
#' @export
#' @param dd dendrogram object
#' @return A dendrogram plot
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @seealso \code{\link{clusterDataDend}}
#' @examples
#' data("fullProcessedProfileLarge", package="PhyloProfile")
#' data <- fullProcessedProfileLarge
#' profileType <- "binary"
#' profiles <- getDataClustering(
#'     data, profileType, var1AggregateBy, var2AggregateBy)
#' distMethod <- "mutualInformation"
#' distanceMatrix <- getDistanceMatrix(profiles, distMethod)
#' clusterMethod <- "complete"
#' dd <- clusterDataDend(distanceMatrix, clusterMethod)
#' getDendrogram(dd)

getDendrogram <- function(dd) {
    if (is.null(dd)) return()
    py <- dendextend::as.ggdend(dd)
    p <- ggplot(py, horiz = TRUE, theme = theme_minimal()) +
        theme(axis.title = element_blank(), axis.text.y = element_blank())
    return(p)
}
