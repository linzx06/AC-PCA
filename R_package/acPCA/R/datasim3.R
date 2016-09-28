#' Simulated data example 3
#'
#' The confounder is unobserved and we only know the primary variable of interest (the biological conditions). 
#' There are 10 biological conditions, each with 3 replicates. 
#' The variation is shared among replicates for half of the genes and not shared for the other genes.
#' Here are details for the data:
#'
#' \itemize{
#'   \item X the N by p data matrix, number of samples=30, number of variables=400
#'   \item Y the N by q confounder matrix, q=30. Y is designed such that samples within the
#'    same biological condition are shrinked together. Details are provided in the user's guide.
#'   \item lab labels for the biological conditions.
#'   \item true_pattern the true underlying latent pattern
#' }
#'
#' @docType data
#' @keywords datasets
#' @name data_example3
#' @format A list with multiple elements
#' @examples
#' load_all()
#' data(data_example3)
NULL