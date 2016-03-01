calkernel <- function(Y, kernel, bandwidth){
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