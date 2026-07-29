#' Get data for calculating distance matrix from phylogenetic profiles
#' @export
#' @usage getDataClustering(data, profileType = "binary", var1AggBy = "max",
#'     var2AggBy = "max")
#' @param data a data frame contains processed and filtered profiles (see
#' ?fullProcessedProfile and ?filterProfileData, ?fromInputToProfile)
#' @param profileType type of data used for calculating the distance matrix.
#' Either "binary" (consider only the presence/absence status of orthlogs),
#' "orthoID" (consider ortholog IDs as values for clustering),
#' "var1"/"var2" for taking values of the additional variables into account.
#' Default = "binary".
#' @param var1AggBy aggregate method for VAR1 (min, max, mean or median).
#' Default = "max".
#' @param var2AggBy aggregate method for VAR2 (min, max, mean or median).
#' Default = "max".
#' @return A wide dataframe contains values for calculating distance matrix.
#' @author Carla Mölbert (carla.moelbert@gmx.de), Vinh Tran
#' (tran@bio.uni-frankfurt.de)
#' @importFrom data.table dcast setDT
#' @seealso \code{\link{fromInputToProfile}}
#' @examples
#' data("finalProcessedProfile", package="PhyloProfile")
#' data <- finalProcessedProfile
#' profileType <- "binary"
#' var1AggregateBy <- "max"
#' var2AggregateBy <- "mean"
#' getDataClustering(data, profileType, var1AggregateBy, var2AggregateBy)

getDataClustering <- function(
    data = NULL, profileType = "binary", var1AggBy = "max", var2AggBy = "max"
) {
    if (is.null(data)) stop("Input data cannot be NULL!")
    presSpec <- geneID <- supertaxon <- orthoID <- NULL
    # Filter out rows where presSpec is 0
    subDataHeat <- data %>% dplyr::filter(presSpec > 0)
    
    # Predict orthoID format (default is "other")
    subDataHeat$orthoID <- as.character(subDataHeat$orthoID)
    idFormat <- ifelse(
        length(subDataHeat$orthoID) > 0 &&
            grepl(
                subDataHeat$supertaxonID[1],
                strsplit(subDataHeat$orthoID[1], '\\|', fixed = TRUE)[[1]][2]
            ), "bionf", "other"
    )
    # Process based on profileType
    if (profileType == "binary") {
        subDataHeat <- subDataHeat %>%
            dplyr::select(geneID, supertaxon, presSpec) %>%
            dplyr::mutate(presSpec = ifelse(presSpec > 0, 1, 0)) %>%
            dplyr::distinct()
        wideData <- data.table::dcast(
            data.table::setDT(subDataHeat),
            geneID ~ supertaxon, value.var = "presSpec"
        )
    } else if (profileType == "orthoID") {
        subDataHeat <- subDataHeat %>%
            dplyr::select(geneID, supertaxon, orthoID) %>%
            dplyr::mutate(
                orthoID = if (idFormat == "bionf") {
                    paste(
                        vapply(
                            strsplit(as.character(orthoID), "\\|"),
                            `[`, character(1), 2
                        ),
                        vapply(
                            strsplit(as.character(orthoID), "\\|"),
                            `[`, character(1), 3
                        ),
                        sep = "#"
                    )
                } else {
                    paste(supertaxon, orthoID, sep = "#")
                }
            ) %>% dplyr::distinct()
        # Assign numeric factor to orthoID
        subDataHeat <- subDataHeat %>%
            dplyr::mutate(orthoID = as.numeric(factor(orthoID)))
        wideData <- data.table::dcast(
            data.table::setDT(subDataHeat),
            geneID ~ supertaxon, value.var = "orthoID", fun.aggregate = sum
        )
    } else {
        var <- profileType
        subDataHeat <- subDataHeat %>%
            dplyr::select(geneID, supertaxon, !!sym(var)) %>%
            dplyr::distinct()
        # Aggregate values by selected method
        aggregateBy <- ifelse(var == "var1", var1AggBy, var2AggBy)
        subDataHeat <- subDataHeat %>%
            dplyr::group_by(geneID, supertaxon) %>%
            dplyr::summarise(
                !!sym(var) := match.fun(aggregateBy)(!!sym(var)),
                .groups = 'drop'
            )
        wideData <- data.table::dcast(
            data.table::setDT(subDataHeat),
            geneID ~ supertaxon, value.var = var
        )
    }
    # Convert to data frame and set geneID as row names
    wideData <- as.data.frame(wideData)
    dat <- wideData[, -1]
    rownames(dat) <- wideData$geneID
    # Replace NA with 0
    dat[is.na(dat)] <- 0
    return(dat)
}

#' Calculate the distance matrix
#' @export
#' @param profiles dataframe contains profile data for distance calculating
#' (see ?getDataClustering)
#' @param method distance calculation method ("euclidean", "maximum",
#' "manhattan", "canberra", "binary", "distanceCorrelation",
#' "mutualInformation" or "pearson" for binary data; "distanceCorrelation" or
#' "mutualInformation" for non-binary data). Default = "mutualInformation".
#' @return A calculated distance matrix for input phylogenetic profiles.
#' @importFrom bioDist mutualInfo
#' @importFrom bioDist cor.dist
#' @importFrom stats dist as.dist
#' @importFrom energy dcor
#' @importFrom stats dist
#' @author Carla Mölbert (carla.moelbert@gmx.de), Vinh Tran
#' (tran@bio.uni-frankfurt.de)
#' @seealso \code{\link{getDataClustering}}
#' @examples
#' data("finalProcessedProfile", package="PhyloProfile")
#' data <- finalProcessedProfile
#' profileType <- "binary"
#' profiles <- getDataClustering(
#'     data, profileType, var1AggregateBy, var2AggregateBy)
#' method <- "mutualInformation"
#' getDistanceMatrix(profiles, method)

getDistanceMatrix <- function(profiles = NULL, method = "mutualInformation") {
    if (is.null(profiles)) stop("Profile data cannot be NULL!")
    profiles <-  profiles[, colSums(profiles != 0) > 0]
    profiles <-  profiles[rowSums(profiles != 0) > 0, ]
    distMethods <- c("euclidean", "maximum", "manhattan", "canberra", "binary")
    if (method %in% distMethods) {
        distanceMatrix <- as.matrix(stats::dist(profiles, method = method))
        colnames(distanceMatrix)<-rownames(distanceMatrix)<-rownames(profiles)
    } else if (method == "distanceCorrelation") {
        n <- seq_len(nrow(profiles))
        matrix <- matrix(0L, nrow = nrow(profiles), ncol = nrow(profiles))
        for (i in seq_len(nrow(profiles))) { # rows
            p_i <- unlist(profiles[i,])
            for (j in seq_len(nrow(profiles))) { # columns
                if (i == j) break
                matrix[i, j] <- energy::dcor(p_i, unlist(profiles[j,]))
            }
        }
        # Swich the value so that the profiles with a high correlation
        # are clustered together
        matrix <- 1 - matrix
        matrix <- as.data.frame(matrix)
        profileNames <- rownames(profiles)
        colnames(matrix) <- profileNames[seq_len(length(profileNames)) - 1]
        rownames(matrix) <- profileNames
        distanceMatrix <- stats::as.dist(matrix)
    } else if (method == "mutualInformation") {
        distanceMatrix <- bioDist::mutualInfo(as.matrix(profiles))
        distanceMatrix <- max(distanceMatrix, na.rm = TRUE) - distanceMatrix
    } else if (method == "pearson") {
        distanceMatrix <-  bioDist::cor.dist(as.matrix(profiles))
    }
    return(as.matrix(distanceMatrix))
}

#' Create a hclust object from the distance matrix
#' @export
#' @param distanceMatrix calculated distance matrix as dist object
#' @param clusterMethod clustering method ("single", "complete",
#' "average" for UPGMA, "mcquitty" for WPGMA, "median" for WPGMC,
#' or "centroid" for UPGMC). Default = "complete".
#' @return An object class hclust generated based on input distance matrix and
#' a selected clustering method.
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @importFrom fastcluster hclust
#' @seealso \code{\link{getDataClustering}},
#' \code{\link{getDistanceMatrix}}, \code{\link{hclust}}
#' @examples
#' data("finalProcessedProfile", package="PhyloProfile")
#' data <- finalProcessedProfile
#' profileType <- "binary"
#' profiles <- getDataClustering(
#'     data, profileType, var1AggregateBy, var2AggregateBy)
#' distMethod <- "mutualInformation"
#' distanceMatrix <- getDistanceMatrix(profiles, distMethod)
#' clusterMethod <- "complete"
#' clusterDataDend(as.dist(distanceMatrix), clusterMethod)

clusterDataDend <- function(distanceMatrix = NULL, clusterMethod = "complete") {
    if (is.null(distanceMatrix)) stop("Distance matrix cannot be NULL!")
    if (!inherits(distanceMatrix,"dist")) stop("Input must be a 'dist' object!")
    # Validate clustering method
    validMethods <- c(
        "complete", "ward.D", "ward.D2", "single", "average", "mcquitty",
        "median", "centroid"
    )
    if (!clusterMethod %in% validMethods) {
        stop(
            paste(
                "Invalid clustering method. Choose one of:",
                paste(validMethods, collapse = ", ")
            )
        )
    }
    # Perform clustering using fastcluster
    dendrogram <- fastcluster::hclust(distanceMatrix, method = clusterMethod)
    return(dendrogram)
}

#' Plot dendrogram tree
#' @export
#' @param dd dendrogram object (see ?clusterDataDend)
#' @return A dendrogram plot for the genes in the input phylogenetic profiles.
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @importFrom ape plot.phylo
#' @seealso \code{\link{clusterDataDend}}
#' @examples
#' data("finalProcessedProfile", package="PhyloProfile")
#' data <- finalProcessedProfile
#' profileType <- "binary"
#' profiles <- getDataClustering(
#'     data, profileType, var1AggregateBy, var2AggregateBy)
#' distMethod <- "mutualInformation"
#' distanceMatrix <- getDistanceMatrix(profiles, distMethod)
#' clusterMethod <- "complete"
#' dd <- clusterDataDend(as.dist(distanceMatrix), clusterMethod)
#' getDendrogram(dd)

getDendrogram <- function(dd = NULL) {
    if (is.null(dd)) stop("Input dendrogram cannot be NULL!")
    p <- ape::plot.phylo(as.phylo(dd))
    return(p)
}
