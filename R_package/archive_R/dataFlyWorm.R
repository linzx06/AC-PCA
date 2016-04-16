#' The modENCODE RNA-Seq data, fly and worm, embryonic stage
#'
#' The modENCODE RNA-Seq data, fly and worm, embryonic stage. For the orthologs 
#' in fly that map to multiple orthologs in worm, we took median to get a one to one match, 
#' resulting in 4831 ortholog paris. Data downloaded from https://www.encodeproject.org/comparative/transcriptome/.
#'
#' \itemize{
#'   \item data_fly expression levels for fly.
#'   \item data_worm expression levels for worm.
#'   \item fly_time time windows for fly. In the unit of hours.
#'   \item worm_time time windows for worm In the unit of hours.
#'   \item X data matrix for fly and worm combined.
#'   To gain robustness, we used the rank across samples within the same species. 
#'   The rank matrix was then scaled to have unit variance.
#'   \item X_time labels for the time windows/points
#'   \item X_species labels for the species
#'   \item Y confounder matrix (see user's guide for details). 
#'   \item fly_gene the fly ortholog
#' }
#'
#' @docType data
#' @keywords datasets
#' @name data_fly_worm
#' @format A list with multiple elements
#' @examples
#' load_all()
#' data(data_fly_worm)
NULL

