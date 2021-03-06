% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clusterProfile.R
\name{getDistanceMatrix}
\alias{getDistanceMatrix}
\title{Calculate the distance matrix}
\usage{
getDistanceMatrix(profiles = NULL, method = "mutualInformation")
}
\arguments{
\item{profiles}{dataframe contains profile data for distance calculating
(see ?getDataClustering)}

\item{method}{distance calculation method ("euclidean", "maximum",
"manhattan", "canberra", "binary", "distanceCorrelation",
"mutualInformation" or "pearson" for binary data; "distanceCorrelation" or
"mutualInformation" for non-binary data). Default = "mutualInformation".}
}
\value{
A calculated distance matrix for input phylogenetic profiles.
}
\description{
Calculate the distance matrix
}
\examples{
data("finalProcessedProfile", package="PhyloProfile")
data <- finalProcessedProfile
profileType <- "binary"
profiles <- getDataClustering(
    data, profileType, var1AggregateBy, var2AggregateBy)
method <- "mutualInformation"
getDistanceMatrix(profiles, method)
}
\seealso{
\code{\link{getDataClustering}}
}
\author{
Carla Mölbert (carla.moelbert@gmx.de), Vinh Tran
(tran@bio.uni-frankfurt.de)
}
