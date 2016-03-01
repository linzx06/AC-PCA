#' Calculate the kernel matrix
#'
#' @param Y the n by q confounder matrix, where n is the number of samples. Missing values in Y should be labeled as NA. 
#' @param kernel the kernel to use: "linear", "gaussian".
#' @param bandwidth bandwidth h for Gaussian kernel. Optional. 
#' @param scaleY scale the columns in Y to unit standard deviation. Default is False.
#' @return The kernel matrix
#' \item{K}{the n by n kernel matrix for Y} 
#' @export
#' @examples
#' load_all()
#' Y <- data_example1$Y 
#' K1 <- calkernel(Y, kernel="linear") ##linear kernel
#' K2 <- calkernel(Y, kernel="gaussian", bandwidth=1) ##Gaussian kernel
calkernel <- function(Y, kernel, bandwidth, scaleY=F){
  Y <- scale(Y, center = F, scale = scaleY)  
  ####missing data
  Y[is.na(Y)] <- mean(Y, na.rm=T)
  
  if (kernel=="linear"){
    K <- tcrossprod(Y)
  } else if (kernel=="gaussian"){
    if (is.null(bandwidth)==T){
      stop("For gaussian kernel, please specify the bandwidth") 
    } else{
      K <- as.matrix(dist(Y, method = "euclidean"))
      K <- exp(-K^2/2/bandwidth^2)  
    }
  } else {
    stop("Please select a valid kernel, linear kernel or gaussian kernel")
  } 
  return(K)
}