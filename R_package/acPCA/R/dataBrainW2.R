#' The human brain exon array dataset, time window 2
#'
#' A subset of 1,000 genes are included for demonstration purpose.
#'
#' \itemize{
#'   \item X the N by p data matrix, N is the number of samples and p is the number of genes.
#'   \item Y the N by q confounder matrix (See implementation details in the user's guide).
#'   \item Yid labels for the individuals, left and right hemispheres from the same donors are treated as different individuals (See implementation details in the user's guide).
#'   \item regions labels for the brain regions
#'   \item hemispheres '1' represents left hemisphere, '3' represents right hemisphere
#'   \item donor_labs labels for the donors
#' }
#'
#' @docType data
#' @keywords datasets
#' @name data_brain_w2
#' @format A list with multiple elements
#' @examples
#' load_all()
#' data(data_brain_w2)
NULL

