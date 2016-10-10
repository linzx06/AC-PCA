#' Simulated data example 1
#'
#' In this simulated example, the confounder is a factor with three levels, corresponding to three groups. 
#' The confounder contributes globally to all variables(genes).
#' Here are details for the data:
#'
#' \itemize{
#'   \item X the N by p data matrix, number of samples=30, number of variables=400
#'   \item Y the N by q confounder matrix, q=3. Y is chosen such that the penalty term equals the between groups sum of squares.
#'   \item lab the labels (used in the plots in the user's guide)
#'   \item group the group labels
#'   \item true_pattern the true underlying latent pattern
#' }
#'
#' @docType data
#' @keywords datasets
#' @name data_example1
#' @format A list with multiple elements
#' @examples
#' load_all()
#' data(data_example1)
NULL