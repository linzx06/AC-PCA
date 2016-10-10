#' Simulated data example 3
#'
#' The confounder is unobserved and we only know the primary variable of interest (the biological conditions). 
#' There are 10 biological conditions, each with 3 replicates. 
#' The variation is shared among replicates for half of the genes and not shared for the other genes.
#' Here are details for the data:
#'
#' \itemize{
#'   \item X the N by p data matrix, number of samples=30, number of variables=400
#'   \item Y the N by q confounder matrix, q=30. For a biological condition, 
#'   treating the 3 replicates as 3 groups, it can be shown that the penalty term 
#'   equals the summation of the between groups sum of squares over the biological conditions.
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