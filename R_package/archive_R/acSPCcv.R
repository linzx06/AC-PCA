#' Perform cross-validation to tune the sparsity parameter c2 in function acSPC
#'
#' @param X the n by p data matrix, where n is the number of samples, p is the number of variables. Missing values in X should be labeled as NA. If a whole sample in X is missing, it should be removed.
#' @param Y the n by q confounder matrix, where n is the number of samples, q is the number of confounding factors. Missing values in Y should be labeled as NA.
#' @param X4Y the "X" used to calculate the empirical Hilbert Schmidt criterion. Default is set to X. Optional.
#' @param c1 tuning parameter. Default is set to v'X4Y'KX4Yv. Optional. 
#' @param c2s a vector of tuning parameters controlling sparsity.
#' @param v_ini the initial v. Recommended to be the estimate of the non-sparse version. 
#' @param v_substract the principal components to be subtracted. A p by k matrix, where k is the number of PCs to be substracted. Optional.
#' @param kernel the kernel to use: "linear", "gaussian".
#' @param bandwidth bandwidth h for Gaussian kernel. Optional. 
#' @param centerX center the columns in X. Default is True.
#' @param scaleX scale the columns in X to unit standard deviation. Default is False.
#' @param scaleY scale the columns in Y to unit standard deviation. Default is False.
#' @param fold the fold number for cross-validation. Default is 10.
#' @param plot True or False. plot=T generates the c2 vs. mean squared error(MSE) plot. Default is True.
#' @param quiet True or False. Output the progress of the program. Default is False. 
#' @param ... other parameters
#' @return Results for tuning the sparsity parameter c2
#' \item{mse}{a vector of MSEs for c2s. Mean is taken for the sum of squared errors in each fold.} 
#' \item{c2s}{the input c2s}
#' \item{best_c2}{c2 with the smallest MSE}
#' \item{mse_sd}{a vector of standard deviations for the MSEs. Standard deviation across the folds is calculated.} 
#' @export
#' @examples
#' load_all()
#' data(data_example5)
#' X <- data_example5$X; Y <- data_example5$Y
#' result1cv <- acPCAcv(X=X, Y=Y, lambdas=seq(0, 2, 0.05), kernel="linear", plot=F, quiet=T)
#' result1 <- acPCA(X=X, Y=Y, lambda=result1cv$best_lambda, kernel="linear")
#' v_ini <- result1$v[,1]
#' resultcv_spc1_coarse <- acSPCcv( X=X, Y=Y, c2s=seq(1, 0, -0.1)*sum(abs(v_ini)), 
#'                                  v_ini=v_ini, kernel="linear") 
#' resultcv_spc1_fine <- acSPCcv( X=X, Y=Y, c2s=seq(0.7, 0.4, -0.02)*sum(abs(v_ini)), 
#'                                v_ini=v_ini, kernel="linear") 
#' result_spc1 <- acSPC( X=X, Y=Y, c2=resultcv_spc1_fine$best_c2, v_ini=v_ini, kernel="linear") 
acSPCcv <- function( X, Y, c1=NULL, c2s, v_ini, v_substract=NULL, X4Y=NULL, kernel=c("linear", "gaussian"), bandwidth=NULL, centerX=T, scaleX=F, scaleY=F, maxiter=25, delta=10^-4, fold=10, plot=T, quiet=F){
  ####if Y is a vector, change it to a matrix
  if (is.null(dim(Y))){
    Y <- matrix(Y, ncol=1)
  }
  ####check whether all elements in c2s are non-negative
  if (sum(c2s<0)){
    stop("Values in c2s must be non-negative")
  }
  ####check whether a whole row in X is missing
  Xmis <- apply(X, 1, function(row){sum(!is.na(row))})
  if (sum(Xmis==0)){
    stop(paste("The following samples in X is missing, please remove them in X and Y: rows ", paste(which(Xmis==0), collapse =  " "), sep=""))
  }
  ####check whether the number of samples in X and Y match
  if (dim(X)[1]!=dim(Y)[1]){
    stop("The numbers of samples in X ( nrow(X) ) and Y ( nrow(Y) ) do not match")
  }
  ####check whether the number of samples in X4Y and Y match
  if (!is.null(X4Y)){
    if (dim(X4Y)[1]!=dim(Y)[1]){
      stop("The numbers of samples in X4Y ( nrow(X4Y) ) and Y ( nrow(Y) ) do not match")
    }  
  }
  
  X <- scale(X, center = centerX, scale = scaleX)  
  Y <- scale(Y, center = F, scale = scaleY)  
  ####missing data
  X[is.na(X)] <- mean(X, na.rm=T)
  Y[is.na(Y)] <- mean(Y, na.rm=T)
  
  K <- calkernel(Y, kernel, bandwidth)
  if (is.null(c1)){
    if (is.null(X4Y)){
      Xv <- X%*%as.matrix(v_ini)
      c1 <- as.numeric(crossprod(Xv, K%*%Xv))  
    } else {
      Xv <- X4Y%*%as.matrix(v_ini)
      c1 <- as.numeric(crossprod(Xv, K%*%Xv)) 
    }
  }
  ####labels for the fold
  lab <- matrix(sample(fold, prod(dim(X)), replace=T), nrow=nrow(X), ncol=ncol(X))
  difsqr <- matrix(nrow=fold, ncol=length(c2s))
  for (f in 1:fold){
    if (!quiet){
      print(paste("Running fold", f))  
    }
    X_cv <- X
    X_cv[lab==f] <- mean(X_cv[lab!=f])
    for (i in 1:length(c2s)){
      if (!quiet){
        if (i==length(c2s)){
          cat("\r", round(i/length(c2s)*100), "%", "completed\n")
        } else {
          cat("\r", round(i/length(c2s)*100), "%", "completed")  
        }
        
      }
      c2 <- c2s[i]
      if (c2==0){
        difsqr[f, i] <- sum(X[lab==f]^2)  
      } else {
        if (is.null(X4Y)){
          tmp <- acSPC(X=X_cv, Y=Y, c1=c1, c2=c2, v_ini=v_ini, v_substract=v_substract, X4Y=X, kernel=kernel, bandwidth=bandwidth, centerX=centerX, scaleX=scaleX, scaleY=scaleY, maxiter=maxiter, delta=delta, filter=F)   
        } else {
          tmp <- acSPC(X=X_cv, Y=Y, c1=c1, c2=c2, v_ini=v_ini, v_substract=v_substract, X4Y=X4Y, kernel=kernel, bandwidth=bandwidth, centerX=centerX, scaleX=scaleX, scaleY=scaleY, maxiter=maxiter, delta=delta, filter=F)
        }
        v <- tmp$v; u <- tmp$u
        tmp <- tcrossprod((sum((X_cv%*%v)*u)*u), v)
        difsqr[f, i] <- sum((tmp[lab==f] - X[lab==f])^2)
        gc()  
      }
    }
  }
  mse <- apply(difsqr, 2, mean)
  mse_sd <- apply(difsqr, 2, sd)
  best_c2 <- c2s[which.min(mse)]
  if (plot==T){
    plot(c2s, mse, ylab="mean squared error(MSE)", xlab=expression('Sparsity parameter c'[2]), pch=20, cex=0.5, main="Cross-validation", cex.lab=1.5)
    lines(c2s, mse, col="red", lty=2)
    abline(v=best_c2, col="blue")
    legend( "topright", legend=expression('Best c'[2]), lty=1, col="blue" )  
  }
  return( list(mse=mse, mse_sd=mse_sd, c2s=c2s, best_c2=best_c2) ) 
}
