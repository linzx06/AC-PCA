#' Perform sparse AC-PCA for variable selection
#'
#' @param X the n by p data matrix, where n is the number of samples. Missing values in X should be labeled as NA. If a whole sample in X is missing, it should be removed.
#' @param Y the n by q confounder matrix, where n is the number of samples. Missing values in Y should be labeled as NA. 
#' @param X4Y the "X" used to calculate the empirical Hilbert Schmidt criterion. Default is set to X. Optional.
#' @param c1 tuning parameter. Default is set to v'X4Y'KX4Yv. Optional. 
#' @param c2 tuning parameter controlling sparsity.
#' @param v_ini the initial v. Recommended to be the estimate of the non-sparse version. 
#' @param v_substract the principal components to be subtracted. A p by k matrix, where k is the number of PCs to be substracted. Optional.
#' @param kernel the kernel to use: "linear", "gaussian".
#' @param bandwidth bandwidth h for Gaussian kernel. Optional. 
#' @param centerX center the columns in X. Default is True.
#' @param scaleX scale the columns in X to unit standard deviation. Default is False.
#' @param scaleY scale the columns in Y to unit standard deviation. Default is False.
#' @param ... other parameters
#' @return Results for sparse AC-PCA
#' \item{v}{the sparse principal component} 
#' \item{u}{the u vector}
#' \item{converge}{whether the algorithm converged}
#' @export
#' @examples
#' load_all()
#' data(data_example5)
#' X <- data_example5$X ###the data matrix
#' Y <- data_example5$Y
#'
#' result1 <- acPCA(X=X, Y=Y, lambda=1, kernel="linear")
#' v_ini <- result1$v[,1]
#' result_spc1 <- acSPC( X=X, Y=Y, c2=0.5*sum(abs(v_ini)), 
#'                       v_ini=v_ini, kernel="linear")
#' ###examples with more details are provided in the function acSPCcv
acSPC <- function( X, Y, X4Y=NULL, c1=NULL, c2, v_ini, v_substract=NULL, kernel=c("linear", "gaussian"), bandwidth=NULL, centerX=T, scaleX=F, scaleY=F, maxiter=50, delta=10^-8, filter=T){
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
  return(list(v=v, u=u, converge=converge))
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
