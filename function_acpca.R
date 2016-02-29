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

calAv <- function(v, args){
  X <- args$X
  K <- args$K
  lambda <- args$lambda
  return( crossprod(X, (diag(dim(K)[1])-lambda*K)%*%(X%*%matrix(v, ncol=1))) )
  #return( crossprod(X, (diag(dim(K)[1])-lambda*K)%*%X%*%matrix(v, ncol=1)) )
}

calAv1 <- function(v, args){
  X <- args$X
  K <- args$K
  lambda <- args$lambda
  #return( crossprod(X, (diag(dim(K)[1])-lambda*K)%*%(X%*%matrix(v, ncol=1))) )
  return( crossprod(X, (diag(dim(K)[1])-lambda*K)%*%X%*%matrix(v, ncol=1)) )
}

calAv2 <- function(v, args) {
  X <- args$X
  K <- args$K
  lambda <- args$lambda
  ## X'(I - lambda * K)Xv
  Xv = X %*% v
  as.numeric(crossprod(X, Xv - lambda * (K %*% Xv)))
}

acPCA<- function(X, Y, lambda, centerX=T, scaleX=F, scaleY=F, nPC=2, kernel=c("linear", "gaussian"), bandwidth=NULL){
  ####check whether a whole row in X is missing
  Xmis <- apply(X, 1, function(row){sum(!is.na(row))})
  if (sum(Xmis==0)){
    stop(paste("The following samples in X is missing, please remove them in X and Y: rows ", paste(which(Xmis==0), collapse =  " "), sep=""))
  }
  ####check whether the number of samples in X and Y match
  if (dim(X)[1]!=dim(Y)[1]){
    stop("The numbers of samples in X ( nrow(X) ) and Y ( nrow(Y) ) do not match")
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
  v <- try(matrix(eigs_sym(calAv, k=nPC, which = "LA", n=p, args=list(X=X, K=K, lambda=lambda))$vectors, ncol=nPC), silent=T)
  if (class(v)=="try-error"){
    v <- try(matrix(eigs_sym(calAv1, k=nPC, which = "LA", n=p, args=list(X=X, K=K, lambda=lambda))$vectors, ncol=nPC), silent=T)
    if (class(v)=="try-error"){
      v <- try(matrix(eigs_sym(calAv, k=nPC, which = "LA", n=p, args=list(X=X, K=K, lambda=lambda*(1+10^-8)))$vectors, ncol=nPC), silent=T) 
      if (class(v)=="try-error"){
        v <- try(matrix(eigs_sym(crossprod(X, (diag(dim(K)[1])-lambda*K)%*%X), k=nPC, which = "LA")$vectors, ncol=nPC), silent=T)
        if (class(v)=="try-error"){
          stop("Numerical issue on lambda, try another lambda close to the input lambda")
        }
      }
    } 
  }
  ####the projection, Xv
  Xv <- X%*%v
  
  return(list(Xv=Xv, v=v))
}

acPCAcv <- function(X, Y, lambdas, centerX=T, scaleX=F, scaleY=F, nPC=2, fold=5, foldlab=NULL, kernel=c("linear", "gaussian"), bandwidth=NULL, perc=0.05, plot=T, quiet=F){
  ####check whether a whole row in X is missing
  Xmis <- apply(X, 1, function(row){sum(!is.na(row))})
  if (sum(Xmis==0)){
    stop(paste("The following samples in X is missing, please remove them in X and Y: rows ", paste(which(Xmis==0), collapse =  " "), sep=""))
  }
  ####check whether the number of samples in X and Y match
  if (dim(X)[1]!=dim(Y)[1]){
    stop("The numbers of samples in X ( nrow(X) ) and Y ( nrow(Y) ) do not match")
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
  return(list(loss=loss, lambdas=lambdas, kernel=kernel, best_lambda=best_lambda))
}

acSPC <- function( X, Y, c1, c2, v_ini, X4Y=NULL, kernel=c("linear", "gaussian"), bandwidth=NULL, centerX=T, scaleX=F, scaleY=F, maxiter=50, delta=10^-8, filter=T){
  ####if Y is a vector, change it to a matrix
  if (is.null(dim(Y))){
    Y <- matrix(Y, ncol=1)
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
  ####Decomposition on the kernel matrix to get the M matrix
  if (kernel=="linear" & dim(Y)[2]<=dim(Y)[1]){
    if (is.null(X4Y)){
      M <- crossprod(Y, X)  
    } else {
      M <- crossprod(Y, X4Y) 
    }
  } else {
    K <- calkernel(Y, kernel, bandwidth)
    tmp <- svd(K)
    U <- tmp$u; S <- tmp$d
    posS <- which(S>0)
    Delta <- U[,posS]%*%diag(sqrt(S[posS]))
    if (is.null(X4Y)){
      M <- crossprod(Delta, X)  
    } else {
      M <- crossprod(Delta, X4Y) 
    }
  }
  ####svd on M
  tmp = svd(M);
  S <- tmp$d; V <- tmp$v
  ####initialization
  v <- v_ini;
  valA <- c();
  devia <- c();
  iter <- 1; converge <- 0;
  val <- -Inf;
  ####iteration
  while (converge==0 & iter<= maxiter){
    ####update u
    tmp <- X%*%v;
    u <- tmp/getnorm2(tmp);
    ####update v
    utx <- t(u)%*%X;
    tmin <- 0; tmax <- getnorm2(utx);
    tmpproj <- proj( M, S, V, utx, c1, c2, v, tmin, tmax, iter)
    vnew <- tmpproj$v;
    #vnew[abs(vnew)<=max(abs(vnew))/10^5] <- 0
    #check convergence
    devia <- c(devia, getnorm2(v-vnew))
    valnew <- tmpproj$t;
    if (valnew - val <= delta){
      converge <- 1
    }
    iter <- iter + 1;
    v <- vnew
    val <- valnew
    valA <- c(valA, val)    
  }
  ####Because of the bioconvexity of the problem and numerical issues, sometimes there are a lot of small entries in v very close to 0
  if (filter==T){
    v <- as.matrix(topthresfilter(v, top=0.01, alpha=0.5*10^-4))
  }
  return(list(v=v, u=u, objec=valA, devia=devia, converge=converge))
}

acSPCcv <- function( X, Y, c1, c2s, v_ini, X4Y=NULL, kernel=c("linear", "gaussian"), bandwidth=NULL, centerX=T, scaleX=F, scaleY=F, maxiter=25, delta=10^-4, fold=5, plot=T, quiet=F){
  ####if Y is a vector, change it to a matrix
  if (is.null(dim(Y))){
    Y <- matrix(Y, ncol=1)
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
      c2 <- c2s[i]
      write.table(c(f, c2), file="tunec2.tmp")
      if (is.null(X4Y)){
        tmp <- acSPC(X=X_cv, Y=Y, c1=c1, c2=c2, v_ini=v_ini, X4Y=X, kernel=kernel, bandwidth=bandwidth, centerX=centerX, scaleX=scaleX, scaleY=scaleY, maxiter=maxiter, delta=delta, filter=F)   
      } else {
        tmp <- acSPC(X=X_cv, Y=Y, c1=c1, c2=c2, v_ini=v_ini, X4Y=X4Y, kernel=kernel, bandwidth=bandwidth, centerX=centerX, scaleX=scaleX, scaleY=scaleY, maxiter=maxiter, delta=delta, filter=F)
      }
      v <- tmp$v; u <- tmp$u
      tmp <- tcrossprod((sum((X_cv%*%v)*u)*u), v)
      difsqr[f, i] <- sum((tmp[lab==f] - X[lab==f])^2)
      gc()
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

acSPCcvM <- function( X, Y, c1=NULL, c2s, v_ini, v_substract=NULL, X4Y=NULL, kernel=c("linear", "gaussian"), bandwidth=NULL, centerX=T, scaleX=F, scaleY=F, maxiter=25, delta=10^-4, fold=5, plot=T, quiet=F){
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

acSPCM <- function( X, Y, c1=NULL, c2, v_ini, v_substract=NULL, X4Y=NULL, kernel=c("linear", "gaussian"), bandwidth=NULL, centerX=T, scaleX=F, scaleY=F, maxiter=50, delta=10^-8, filter=T){
  ####if Y is a vector, change it to a matrix
  if (is.null(dim(Y))){
    Y <- matrix(Y, ncol=1)
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
  ####Decomposition on the kernel matrix to get the M matrix
  if (kernel=="linear" & dim(Y)[2]<=dim(Y)[1]){
    if (is.null(X4Y)){
      M <- crossprod(Y, X)  
    } else {
      M <- crossprod(Y, X4Y) 
    }
  } else {
    tmp <- svd(K)
    U <- tmp$u; S <- tmp$d
    posS <- which(S>0)
    Delta <- U[,posS]%*%diag(sqrt(S[posS]))
    if (is.null(X4Y)){
      M <- crossprod(Delta, X)  
    } else {
      M <- crossprod(Delta, X4Y) 
    }
  }
  ####svd on M
  tmp = svd(M);
  S <- tmp$d; V <- tmp$v
  
  ####substract the first several PCs (i.e. v_substract), if v_substract is provided
  if (!is.null(v_substract)){
    if(is.null(dim(v_substract))){
      v_substract <- as.matrix(v_substract, ncol=1)
    }
    ###substract the PCs in v_substract one by one
    npc <- dim(v_substract)[2]
    for (pc in 1:npc){
      X <- X - tcrossprod( X%*%matrix(v_substract[,pc], ncol=1), matrix(v_substract[,pc], ncol=1) )   
    }
  }
  
  ####initialization
  v <- v_ini;
  valA <- c();
  devia <- c();
  iter <- 1; converge <- 0;
  val <- -Inf;
  ####iteration
  while (converge==0 & iter<= maxiter){
    ####update u
    tmp <- X%*%v;
    u <- tmp/getnorm2(tmp);
    ####update v
    utx <- t(u)%*%X;
    tmin <- 0; tmax <- getnorm2(utx);
    tmpproj <- proj( M, S, V, utx, c1, c2, v, tmin, tmax, iter)
    vnew <- tmpproj$v;
    #vnew[abs(vnew)<=max(abs(vnew))/10^5] <- 0
    #check convergence
    devia <- c(devia, getnorm2(v-vnew))
    valnew <- tmpproj$t;
    if (valnew - val <= delta){
      converge <- 1
    }
    iter <- iter + 1;
    v <- vnew
    val <- valnew
    valA <- c(valA, val)    
  }
  ####Because of the bioconvexity of the problem and numerical issues, sometimes there are a lot of small entries in v very close to 0
  if (filter==T){
    v <- as.matrix(topthresfilter(v, top=0.01, alpha=0.5*10^-4))
  }
  return(list(v=v, u=u, objec=valA, devia=devia, converge=converge))
}

proj <- function( Y, S, V, utx, c1, c2, v_ini, tmin, tmax, iternum){
  epsilon <- 10^-3/iternum; epsilon1 <- 10^-3/iternum; epsilon2 <- 10^-4/iternum; maxiter <- round(50*sqrt(iternum))
  tlow <- tmin; tup <- tmax; v <- v_ini
  while (tup - tlow >= tup*epsilon){
    t <- (tup + tlow)/2;
    iter <- 1; converge <- 0; 
    l <- 0;  
    while (converge==0 & iter<=maxiter){
      #project to l1 ball
      if (sum(abs(v)) > c2){
        v <- proj_l1(v, c2);  
      }
      #project to l2 ball
      if (getnorm2(v) > 1){
        v <- proj_l2( v, 1);  
      }
      #project to hyperplane
      if (t - utx%*%v > 0){
        v <- proj_hp( v, t(utx), t)
      }
      #project to elliptic
      P <- crossprod(V, v)^2; 
      if (getnorm2(Y%*%v) > sqrt(c1)){
        l <- proj_l2g_getl( S, P, c1)
        #v <- v - (V%*%diag(1/(1 + 1/(l*S^2))))%*%crossprod(V,v)
        ####when there is only one non-zero singular value
        if (length(S)==1){
          v <- v - V*(1/(1 + 1/(l*S^2)))*as.numeric(crossprod(V,v))
        } else {
          v <- v - V%*%(diag(1/(1 + 1/(l*S^2)))%*%crossprod(V,v))  
        }
        
      }
      #check constrain
      if (t - utx%*%v <= t*epsilon & sum(abs(v)) <= c2*(1 + epsilon1) & getnorm2(v) <= 1 + epsilon2){
        converge <- 1  
      }
      iter <- iter + 1
    }
    
    if (converge == 1){
      tlow <- utx%*%v
      vcon <- v
    } else {
      tup <- t  
    }
  }
  return(list(v = vcon, t = utx%*%vcon))
}

proj_l2 <- function(v, c){
  return(v/getnorm2(v)*c)    
}

proj_hp <- function( v, w, t){
  t <- as.numeric(t)
  return( v - w*((sum(w*v) - t)/getnorm2(w)^2) )
}

proj_l2g_getl <- function( S, P, c){
  S2 <- S^2;
  epsilon <- 10^-10; maxiter <- 50;
  l <- 10; lnew <- 0
  iter <- 1;
  while (abs(l - lnew)>=epsilon && iter <= maxiter){
    l <- lnew  
    d1 <- sum(P*S2/((1 + l*S2)^2)) - c
    d2 <- -2*sum( P*S^4/ ( (1 + l*S2 )^3 ) )
    lnew <- l - d1/d2
  }
  return(lnew)  
}

proj_l1 <- function(v, c){
  if (sum(abs(v)) < c){
    return(v)
  }
  u <- sort(abs(v), decreasing=T)
  sv <- cumsum(u)
  rho <- max(which(u > (sv - c) / (1:length(u))))
  theta <- max(0, (sv[rho] - c) / rho)
  return( sign(v) * pmax(abs(v) - theta, 0) )
}

getnorm2 <- function(vector){
  return(sqrt(sum(vector^2)))
}

topthresfilter <- function(v, top, alpha){
  thres <- mean(sort(abs(v), decreasing = T)[1:round(length(v)*top)])*alpha
  v[abs(v)<=thres] <- 0
  return(v)
}

calK <- function(X, h){
  K <- diag(length(X))
  for (i in 1:length(X)){
    for (j in 1:length(X)){
      K[i, j] <- exp(-1/h^2*(X[i]-X[j])^2)
    }    
  }
  return(K)
}