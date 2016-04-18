#' Perform AC-PCA for simultaneous dimension reduction and adjustment for confounding variation
#'
#' @param X the n by p data matrix, where n is the number of samples, p is the number of variables. Missing values in X should be labeled as NA. If a whole sample in X is missing, it should be removed.
#' @param Y the n by q confounder matrix, where n is the number of samples, q is the number of confounding factors. Missing values in Y should be labeled as NA. 
#' @param lambda the tuning parameter, non-negative.
#' @param centerX center the columns in X. Default is True.
#' @param scaleX scale the columns in X to unit standard deviation. Default is False.
#' @param scaleY scale the columns in Y to unit standard deviation. Default is False.
#' @param nPC number of principal components to compute
#' @param kernel the kernel to use: "linear", "gaussian".
#' @param bandwidth bandwidth h for Gaussian kernel. Optional. 
#' @return The principal components and the projected data
#' \item{v}{the principal components, p by nPC matrix} 
#' \item{Xv}{the projected data, i.e. X times v}
#' \item{...}{Input parameters for the function}
#' @export
#' @examples
#' load_all()
#' data(data_example1)
#' X <- data_example1$X; Y <- data_example1$Y 
#' 
#' ##use linear kernel
#' result_linear <- acPCA(X=X, Y=Y, lambda=0.5, kernel="linear", nPC=2) 
#' 
#' ##use Gaussian kernel
#' result_gaussian <- acPCA(X=X, Y=Y, lambda=0.5, kernel="gaussian", bandwidth=1, nPC=2)
#' 
#' ##first tune lambda, and then use the best lambda
#' result_cv <- acPCAcv(X=X, Y=Y, lambdas=seq(0, 1, 0.05), kernel="linear", nPC=2, plot=T)
#' result <- acPCA(X=X, Y=Y, lambda=result_cv$best_lambda, kernel="linear", nPC=2)
acPCA <- function(X, Y, lambda, centerX=T, scaleX=F, scaleY=F, nPC=2, kernel=c("linear", "gaussian"), bandwidth=NULL){
  ####check whether a whole row in X is missing
  Xmis <- apply(X, 1, function(row){sum(!is.na(row))})
  if (sum(Xmis==0)){
    stop(paste("The following samples in X is missing, please remove them in X and Y: rows ", paste(which(Xmis==0), collapse =  " "), sep=""))
  }
  ####check whether the number of samples in X and Y match
  if (dim(X)[1]!=dim(Y)[1]){
    stop("The numbers of samples in X ( nrow(X) ) and Y ( nrow(Y) ) do not match")
  }
  ####check whether lambda is non-negative
  if (lambda<0){
    stop("lambda should be non-negative")
  }
  nsam <- dim(X)[1] 
  p <- dim(X)[2]
  if (is.null(dim(Y))){
    Y <- matrix(Y, ncol=1)
  }
  X <- scale(X, center = centerX, scale = scaleX)  
  Y <- scale(Y, center = F, scale = scaleY)  
  ####missing data
  X[is.na(X)] <- mean(X, na.rm=T)
  Y[is.na(Y)] <- mean(Y, na.rm=T)
  
  ####calculate kernel matrix for Y
  K <- calkernel(Y, kernel, bandwidth)
  
  ####the rotation, v
  v <- matrix(eigs_sym(calAv, k=nPC, which = "LA", n=p, args=list(X=X, K=K, lambda=lambda))$vectors, ncol=nPC)
  ####the projection, Xv
  Xv <- X%*%v
  return(list(Xv=Xv, v=v, lambda=lambda, kernel=kernel, bandwidth=bandwidth))
}

calAv <- function(v, args){
  X <- args$X
  K <- args$K
  lambda <- args$lambda
  return( crossprod(X, (diag(dim(K)[1])-lambda*K)%*%(X%*%matrix(v, ncol=1))) )
}