#' Perform AC-PCA for simultaneous dimension reduction and adjustment for confounding variation
#'
#' @param X the n by p data matrix, where n is the number of samples, p is the number of variables. Missing values in X should be labeled as NA. If a whole sample in X is missing, it should be removed.
#' @param Y the n by q confounder matrix, where n is the number of samples, q is the number of confounding factors. Missing values in Y should be labeled as NA. 
#' @param lambda the tuning parameter, non-negative.
#' @param nPC number of principal components to compute
#' @param eval True or False. eval=T evaluates the significance of the PCs. Default is F. 
#' @param numPerm the number of permutations to evaluate the significance of the PCs. Default is 100. 
#' @param alpha the significance level. Default is 0.05. If the kth PC is not significant, we don't consider the PCs after it.
#' If the eigenvalue and variance explained by the PCs give inconsistent results, we choose the maximum number of significant PCs.
#' @param plot True or False. plot=T generates the plots. Default is True.
#' @param centerX center the columns in X. Default is True.
#' @param centerY center the columns in Y. Default is True.
#' @param scaleX scale the columns in X to unit standard deviation. Default is False.
#' @param scaleY scale the columns in Y to unit standard deviation. Default is False.
#' @param kernel the kernel to use: "linear", "gaussian".
#' @param bandwidth bandwidth h for Gaussian kernel. Optional. 
#' @return The principal components and the projected data
#' \item{v}{the principal components, p by nPC matrix} 
#' \item{Xv}{the projected data, i.e. X times v}
#' \item{eigenvalueX}{eigenvalues for the PCs}
#' \item{varianceX}{variance explained by the PCs}
#' \item{varianceX_perc}{percentage of total variance in X explained by the PCs. If eval=F, NA is returned.}
#' \item{eigenvalueXperm}{eigenvalues for the PCs, permutation. If eval=F, NA is returned.}
#' \item{varianceXperm}{variance explained by the PCs, permutation. If eval=F, NA is returned.}
#' \item{sigPC}{the significant PCs. If eval=F, NA is returned.}
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
#' # linear kernel, the penalty term has the ANOVA interpretation
#' result_tune_linear <- acPCAtuneLambda(X=X, Y=Y, nPC=2, lambdas=seq(0, 20, 0.05), kernel="linear", anov=T)
#' result_linear <- acPCA(X=X, Y=Y, nPC=2, lambda=result_tune_linear$best_lambda, nPC=2, kernel="linear")
#'
acPCA <- function(X, Y, nPC, lambda, eval=F, numPerm=100, alpha=0.05, plot=T, centerX=T, centerY=T, scaleX=F, scaleY=F, kernel=c("linear", "gaussian"), bandwidth=NULL){
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
  ####check alpha
  if (alpha<=0 | alpha>1){
    stop("alpha must be between 0 and 1")
  }
  nsam <- dim(X)[1] 
  p <- dim(X)[2]
  if (is.null(dim(Y))){
    Y <- matrix(Y, ncol=1)
  }
  
  ####missing data
  X[is.na(X)] <- mean(X, na.rm=T)
  Y[is.na(Y)] <- mean(Y, na.rm=T)
  X <- scale(X, center = centerX, scale = scaleX)  
  Y <- scale(Y, center = centerY, scale = scaleY)  

  ####calculate kernel matrix for Y
  K <- calkernel(Y, kernel, bandwidth)
  
  ####AC-PCA
  result_acpca <- eigs_sym(calAv, k=nPC, which = "LA", n=p, args=list(X=X, K=K, lambda=lambda))
  v <- matrix(result_acpca$vectors, ncol=nPC)
  ####eigenvalues
  eigenX <- result_acpca$values
  ####the projection, Xv
  Xv <- X%*%v
  ####variance in Xv
  varX <- apply(Xv, 2, var)
  ####total variance in X
  totvar <- sum(apply(X, 2, var))
  ####percentage
  varX_perc <- varX/totvar
  
  if (eval){
    eigenXperm <- c()
    varXperm <- c()
    varXperm_perc <- c()
    for (i in 1:numPerm){
      Xperm <- X[sample(nrow(X)),]
      tmp <- eigs_sym(calAv, k = nPC, which = "LA", n = p, 
                      args = list(X = Xperm, K = K, lambda = lambda))
      eigenXperm <- rbind(eigenXperm, tmp$values)
      varPC <- apply(Xperm%*%tmp$vectors, 2, var)
      varXperm <- rbind(varXperm, varPC)
      varXperm_perc <- rbind(varXperm_perc, varPC/totvar)
    }
    if (plot){
      par(mfrow=c(1, 2))
      ylimrange <- c(min(c(as.numeric(eigenXperm), eigenX)),
                     max(c(as.numeric(eigenXperm), eigenX)) )
      boxplot(eigenXperm, ylim=ylimrange, xlab="PC", ylab="Eigenvalue", main="")
      points(1:nPC, eigenX, col="red", pch = 4, cex = 1.5, lwd=1.5)
      
      ylimrange <- c(min(c(as.numeric(varXperm), varX)),
                     max(c(as.numeric(varXperm), varX)) )
      boxplot(varXperm, ylim=ylimrange, xlab="PC", ylab="Variance", main="")
      points(1:nPC, varX, col="red", pch = 4, cex = 1.5, lwd=1.5) 
      
      mtext("Data: red cross; Permutation: boxplot", side = 3, outer = TRUE, line=-2)
    }
    labs1 <- which(eigenX < apply(eigenXperm, 2, function(x){quantile(x, 1-alpha)}))
    if (length(labs1)==0){
      sigPC1 <- 0
    } else {
      sigPC1 <- min(labs1)-1
    }
    labs2 <- which(varX < apply(varXperm, 2, function(x){quantile(x, 1-alpha)}))
    if (length(labs2)==0){
      sigPC2 <- 0
    } else {
      sigPC2 <- min(labs2)-1
    }
    sigPC <- max(sigPC1, sigPC2)
  } else {
    sigPC <- NA
    eigenXperm <- NA
    varXperm <- NA
    varXperm_perc <- NA
  }
  return(list(Xv=Xv, v=v, sigPC=sigPC, eigenX=eigenX, 
              varX=varX, varX_perc=varX_perc,
              eigenXperm=eigenXperm, varXperm=varXperm, varXperm_perc=varXperm_perc,
              lambda=lambda, kernel=kernel, bandwidth=bandwidth))
}

calAv <- function(v, args){
  X <- args$X
  K <- args$K
  lambda <- args$lambda
  return( crossprod(X, (diag(dim(K)[1])-lambda*K)%*%(X%*%matrix(v, ncol=1))) )
}