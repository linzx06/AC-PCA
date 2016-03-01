acSPCcv <- function( X, Y, c1=NULL, c2s, v_ini, v_substract=NULL, X4Y=NULL, kernel=c("linear", "gaussian"), bandwidth=NULL, centerX=T, scaleX=F, scaleY=F, maxiter=25, delta=10^-4, fold=5, plot=T, quiet=F){
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
          tmp <- acSPCM(X=X_cv, Y=Y, c1=c1, c2=c2, v_ini=v_ini, v_substract=v_substract, X4Y=X, kernel=kernel, bandwidth=bandwidth, centerX=centerX, scaleX=scaleX, scaleY=scaleY, maxiter=maxiter, delta=delta, filter=F)   
        } else {
          tmp <- acSPCM(X=X_cv, Y=Y, c1=c1, c2=c2, v_ini=v_ini, v_substract=v_substract, X4Y=X4Y, kernel=kernel, bandwidth=bandwidth, centerX=centerX, scaleX=scaleX, scaleY=scaleY, maxiter=maxiter, delta=delta, filter=F)
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
