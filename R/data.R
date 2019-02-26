#' An example of a taxonomy tree in newick format.
#'
#' @format A data frame with only one entry
#' \describe{
#'     \item{V1}{ tree in newick format}
#' }
"pp_tree"

#' An example of a taxonomy matrix.
#'
#' @format A data frame with 10 rows and 162 variables:
#' \itemize{
#'     \item{abbrName}{ e.g. "ncbi10090"}
#'     \item{ncbiID}{ e.g. "10090"}
#'     \item{fullName}{ e.g. "Mus musculus"}
#'     \item{strain}{ e.g. "10090"}
#'     ...
#' }
"pp_taxonomy_matrix"

#' An example of a raw long input file.
#'
#' @format A data frame with 20 rows and 5 variables:
#' \itemize{
#'     \item{geneID}{ Seed or ortholog group ID, e.g. "OG_1017"}
#'     \item{ncbiID}{ Taxon ID, e.g. "ncbi176299"}
#'     \item{orthoID}{ Ortholog ID, e.g. "A.fabrum@176299@1582"}
#'     \item{FAS}{ First additional variable}
#'     \item{traceability}{ Second additional variable}
#' }
"main_long_raw"

#' An example of a raw long input file together with taxonomy info.
#'
#' @format A data frame with 20 rows and 12 variables:
#' \itemize{
#'     \item{geneID}{ Seed or ortholog group ID, e.g. "OG_1017"}
#'     \item{ncbiID}{ Taxon ID, e.g. "ncbi176299"}
#'     \item{orthoID}{ Ortholog ID, e.g. "A.fabrum@176299@1582"}
#'     \item{var1}{ First additional variable}
#'     \item{var2}{ Second additional variable}
#'     \item{paralog}{ Number of co-orthologs in the current taxon}
#'     \item{abbrName}{ e.g. "ncbi176299"}
#'     \item{taxonID}{ Taxon ID, e.g. "176299"}
#'     \item{fullName}{ Full taxon name, e.g. "Agrobacterium fabrum str. C58"}
#'     \item{supertaxonID}{ Supertaxon ID (only different than ncbiID in case
#'     working with higher taxonomy rank than input's)}
#'     \item{supertaxon}{ Name of the corresponding supertaxon}
#'     \item{rank}{ Rank of the supertaxon}
#' }
"profile_with_taxonomy"

#' An example of a fully processed phylogenetic profile.
#'
#' @format A data frame with 20 rows and 17 variables:
#' \itemize{
#'     \item{supertaxon}{ Supertaxon name together with its ordered index, e.g.
#'     "1001_Mammalia"}
#'     \item{geneID}{ Seed or ortholog group ID, e.g. "OG_1017"}
#'     \item{ncbiID}{ Taxon ID, e.g. "ncbi10090"}
#'     \item{orthoID}{ Ortholog ID, e.g. "M.musculus@10090@112934"}
#'     \item{var1}{ First additional variable}
#'     \item{var2}{ Second additional variable}
#'     \item{paralog}{ Number of co-orthologs in the current taxon}
#'     \item{abbrName}{ NCBI ID of the ortholog, e.g. "ncbi10090"}
#'     \item{taxonID}{ Taxon ID of the ortholog, in this case: "0"}
#'     \item{fullName}{ Full taxon name of the ortholog, e.g. "Mus musculus"}
#'     \item{supertaxonID}{ Supertaxon ID (only different than ncbiID in case
#'     working with higher taxonomy rank than input's). e.g. "40674"}
#'     \item{rank}{ Rank of the supertaxon, e.g. "class"}
#'     \item{category}{ "cat}
#'     \item{presSpec}{ The percentage of species presenting in each supertaxon}
#'     \item{mVar1}{ Value of the 1. variable after grouping into supertaxon}
#'     \item{mVar2}{ Value of the 2. variable after grouping into supertaxon}
#'     \item{numberSpec}{ Total number of species in each supertaxon}
#' }
"full_processed_profile"

#' An example of a large processed phylogenetic profile.
#'
#' @format A data frame with 80 rows and 9 variables:
#' \itemize{
#'     \item{geneID}{ Seed or ortholog group ID, e.g. "OG_1017"}
#'     \item{supertaxon}{ Supertaxon name together with its ordered index, e.g.
#'     "1001_Mammalia"}
#'     \item{supertaxonID}{ Supertaxon ID (only different than ncbiID in case
#'     working with higher taxonomy rank than input's). e.g. "40674"}
#'     \item{var1}{ First additional variable}
#'     \item{presSpec}{ The percentage of species presenting in each supertaxon}
#'     \item{category}{ "cat}
#'     \item{orthoID}{ Ortholog ID, e.g. "M.musculus@10090@112934"}
#'     \item{var2}{ Second additional variable}
#'     \item{paralog}{ Number of co-orthologs in the current taxon}
#' }
"full_processed_profile_large"
