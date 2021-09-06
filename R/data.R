#' Gene expression dataset of five tissues in three
#' different species
#'
#' An example dataset containing average gene expression values of
#' five tissues in three different species.
#'
#' @docType data
#'
#' @usage avg_counts <- SBF::TissueExprSpecies
#'
#' @format A list containing gene expression data frame from three species.
#'
#' @keywords datasets
#'
#'
#' @examples
#' # load dataset
#' avg_counts <- SBF::TissueExprSpecies
#'
#' # species names
#' names(avg_counts)
#'
#' # dimension of matrices for different species
#' sapply(avg_counts, dim)
#'
#' # names of different tissue types
#' names(avg_counts[[names(avg_counts)[1]]])
"TissueExprSpecies"
