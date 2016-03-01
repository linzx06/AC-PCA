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

calAv <- function(v, args){
  X <- args$X
  K <- args$K
  lambda <- args$lambda
  return( crossprod(X, (diag(dim(K)[1])-lambda*K)%*%(X%*%matrix(v, ncol=1))) )
}

calAv1 <- function(v, args){
  X <- args$X
  K <- args$K
  lambda <- args$lambda
  return( crossprod(X, (diag(dim(K)[1])-lambda*K)%*%X%*%matrix(v, ncol=1)) )
}

calAv2 <- function(v, args) {
  X <- args$X
  K <- args$K
  lambda <- args$lambda
  Xv = X %*% v
  as.numeric(crossprod(X, Xv - lambda * (K %*% Xv)))
}