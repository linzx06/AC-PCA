#' Simulated data example 2
#'
#' In this simulated example, the confounder is a factor with three levels, corresponding to three groups. 
#' Instead of assuming that the confounder affects all variables (example 1), 
#' we assume that it affects only a subset of variables (half of the variables in this example).
#' Here are details for the data:
#'
#' \itemize{
#'   \item X the N by p data matrix, number of samples=30, number of variables=400
#'   \item Y the N by q confounder matrix, q=3, representing three groups of data
#'   \item lab the labels (used in the plots in the user's guide)
#'   \item colors colors of the labels (used in the plots in the user's guide)
#'   \item true_pattern the true underlying latent pattern
#' }
#'
#' @docType data
#' @keywords datasets
#' @name data_example2
#' @format A list with multiple elements
#' @examples
#' load_all()
#' data(data_example2)
NULL