#' Perform cross-validation to tune lambda
#'
#' @param X the n by p data matrix, where n is the number of samples. Missing values in X should be labeled as NA. If a whole sample in X is missing, it should be removed.
#' @param Y the n by q confounder matrix, where n is the number of samples. Missing values in Y should be labeled as NA. 
#' @param lambdas a vector with the tuning parameters, non-negative values.
#' @param centerX center the columns in X. Default is True.
#' @param scaleX scale the columns in X to unit standard deviation. Default is False.
#' @param scaleY scale the columns in Y to unit standard deviation. Default is False.
#' @param nPC number of principal components to compute
#' @param kernel the kernel to use: "linear", "gaussian".
#' @param bandwidth bandwidth h for Gaussian kernel. Optional. 
#' @param fold the fold number for cross-validation. Default is 10.
#' @param foldlab optional. A vector of labels for cross-validation, should take values from 1 to fold. Length of the vector should match with NROW(X).
#' @param perc the best lambda is defined to be the smallest lambda with loss smaller than (max(loss)-min(loss))*perc + min(loss). 
#' @param plot True or False. plot=T generates the diagnosis plot (lambda vs. loss). Default is True.
#' @param quiet True or False. Output the progress of the program. Default is False. 
#' @return Results for cross-validation
#' \item{loss}{a vector with the loss. Same length as lambdas} 
#' \item{best_lambda}{the best lambda after cross-validation}
#' \item{...}{Input parameters for the function}
#' @export
#' @examples
#' load_all()
#' data(data_example1)
#' X <- data_example1$X; Y <- data_example1$Y 
#'
#' ##first tune lambda, and then use the best lambda. Linear kernel
#' result_cv_linear <- acPCAcv(X=X, Y=Y, lambdas=seq(0, 1, 0.05), kernel="linear", nPC=2, plot=T)
#' result_linear <- acPCA(X=X, Y=Y, lambda=result_cv$best_lambda, kernel="linear", nPC=2)
#' 
#' ##Gaussian kernel
#' result_cv_gaussian <- acPCAcv(X=X, Y=Y, lambdas=seq(0, 1, 0.05), kernel="gaussian", bandwidth=1, nPC=2, plot=T)
#' result_gaussian <- acPCA(X=X, Y=Y, lambda=result_cv$best_lambda, kernel="linear", bandwidth=1, nPC=2)
acPCAcv <- function(X, Y, lambdas, centerX=T, scaleX=F, scaleY=F, nPC=2, kernel=c("linear", "gaussian"), bandwidth=NULL, fold=10, foldlab=NULL, perc=0.05, plot=T, quiet=F){
  ####check whether a whole row in X is missing
  Xmis <- apply(X, 1, function(row){sum(!is.na(row))})
  if (sum(Xmis==0)){
    stop(paste("The following samples in X is missing, please remove them in X and Y: rows ", paste(which(Xmis==0), collapse =  " "), sep=""))
  }
  ####check whether the number of samples in X and Y match
  if (dim(X)[1]!=dim(Y)[1]){
    stop("The numbers of samples in X ( nrow(X) ) and Y ( nrow(Y) ) do not match")
  }
  ####check whether the elements in lambdas are non-negative
  if (sum(lambdas<0)){
    stop("All elements in lambdas should be non-negative")
  }
  nsam <- dim(X)[1] 
  p <- dim(X)[2]
  
  if (is.null(foldlab)){
    if (nsam < fold){
      print("Cannot tune lambda when fold > sample size, setting fold=sample size")
      foldlab <- sample(nsam, nsam)
    } else if (nsam==fold){
      foldlab <- sample(nsam, fold)
    } else {
      foldlab <- rep(NA, nsam)
      labs <- 1:nsam
      num <- floor(nsam/fold)
      #res <- nsam - num*fold
      for (i in 1:fold){
        tmp <- sample(1:length(labs), num)
        foldlab[labs[tmp]] <- i
        labs <- labs[-tmp]
      }
      foldlab[labs] <- sample(fold, length(labs))
    }  
  } 
  ####if Y is a vector, change it to a matrix
  if (is.null(dim(Y))){
    Y <- matrix(Y, ncol=1)
  }
  X <- scale(X, center = centerX, scale = scaleX)  
  Y <- scale(Y, center = F, scale = scaleY)  
  ####missing data
  X[is.na(X)] <- mean(X, na.rm=T)
  Y[is.na(Y)] <- mean(Y, na.rm=T)
  Xv <- array( dim=c(length(lambdas), nsam, nPC) )
  
  ####variable with the largest magnitude of loading
  K <- calkernel(Y, kernel, bandwidth)
  maxlab <- maxlab_sign <- matrix(nrow=length(lambdas), ncol=nPC)
  for (llab in 1:length(lambdas)){
    lambda <- lambdas[llab]
    #vA <- matrix(eigs_sym(calAv, k=nPC, which = "LA", n=p, args=list(X=X, K=K, lambda=lambda))$vectors, ncol=nPC)
    vA <- try(matrix(eigs_sym(calAv, k=nPC, which = "LA", n=p, args=list(X=X, K=K, lambda=lambda))$vectors, ncol=nPC), silent=T)
    if (class(vA)=="try-error"){
      vA <- try(matrix(eigs_sym(calAv1, k=nPC, which = "LA", n=p, args=list(X=X, K=K, lambda=lambda))$vectors, ncol=nPC), silent=T)
      if (class(vA)=="try-error"){
        vA <- try(matrix(eigs_sym(calAv2, k=nPC, which = "LA", n=p, args=list(X=X, K=K, lambda=lambda) )$vectors, ncol=nPC), silent=T)
        if (class(vA)=="try-error"){
          vA <- try(matrix(eigs_sym(crossprod(X, (diag(dim(K)[1])-lambda*K)%*%X), k=nPC, which = "LA")$vectors, ncol=nPC))
        }
      } 
    }
    maxlab[llab,] <- apply(abs(vA), 2, which.max)
    if (nPC==1){
      maxlab_sign[llab,] <- sign(vA[maxlab[llab,],])
    } else {
      maxlab_sign[llab,] <- sign(diag(vA[maxlab[llab,],]))  
    }
  }
  
  for (skip in 1:fold){
    if (!quiet){
      print(paste("Running fold", skip)) 
    }
    Xtmp <- X[foldlab!=skip,]
    Ytmp <- as.matrix(Y[foldlab!=skip,])
    ####calculate kernel matrix for Y
    Ktmp <- calkernel(Ytmp, kernel, bandwidth)
    
    for (llab in 1:length(lambdas)){
      lambda <- lambdas[llab]
      #print(lambda)
      v <- try(matrix(eigs_sym(calAv, k=nPC, which = "LA", n=p, args=list(X=Xtmp, K=Ktmp, lambda=lambda))$vectors, ncol=nPC), silent=T)
      if (class(v)=="try-error"){
        v <- try(matrix(eigs_sym(calAv1, k=nPC, which = "LA", n=p, args=list(X=Xtmp, K=Ktmp, lambda=lambda))$vectors, ncol=nPC), silent=T)
        if (class(v)=="try-error"){
          v <- try(matrix(eigs_sym(crossprod(Xtmp, (diag(dim(Ktmp)[1])-lambda*Ktmp)%*%Xtmp), k=nPC, which = "LA")$vectors, ncol=nPC), silent=T)
        } 
      }
      if (class(v)!="try-error"){
        if (nPC==1){
          v <- v*sign(v[maxlab[llab,],])*maxlab_sign[llab,]
          Xv[llab, which(foldlab==skip),] <- X[which(foldlab==skip),]%*%v  
        } else {
          v <- v%*%diag(sign(diag(v[maxlab[llab,],]))*maxlab_sign[llab,])
          Xv[llab, which(foldlab==skip),] <- X[which(foldlab==skip),]%*%v   
        }
      } 
    }
  }
  ###exclude the lambdas with numerical issues
  labna <- apply(Xv, 1, function(xv){sum(is.na(xv))})!=0
  Xv <- Xv[!labna,,]
  lambdas <- lambdas[!labna]
  ###calculate the loss
  loss <- apply(Xv, 1, function(Xv1, K){sum(diag(crossprod(Xv1, K%*%Xv1)))}, calkernel(Y, kernel, bandwidth))
  thres <- (max(loss)-min(loss))*perc + min(loss)
  best_lambda <- lambdas[min(which(loss<=thres))]
  if (plot==T){
    plot(lambdas, loss, ylab="loss", xlab=expression(lambda), pch=20, cex=0.5, main="Cross-validation", cex.lab=1.5)
    lines(lambdas, loss, col="red", lty=2)
    abline(v=best_lambda, col="blue")
    legend( "topright", legend=expression(paste("Best ", lambda)), lty=1, col="blue" ) 
  }
  return(list(loss=loss, best_lambda=best_lambda, lambdas=lambdas, kernel=kernel, bandwidth=bandwidth))
}
