function [ obj ] = acPCA( X, Y, lambda, nPC, kernel, bandwidth, opts)
%% 
% acPCA -- simultaneous dimension reduction and adjustment for confounding 
% variation.

% input 
% X: a the n by p data matrix, where n is the number of samples, p is the 
% number of variables
% Missing values in X should be labeled as NaN. If a whole sample in X is
% missing, it should be removed.
% Y: the n by q confounder matrix, where n is the number of samples, q is 
% the number of confounding factors.  
% Missing values in Y should be labeled as NaN.
% lambda: tuning parameter, non-negative
% nPC: number of principal components to compute
% kernel: the kernel to use, should be either 'linear', 'gaussian'.
% bandwidth: bandwidth for gaussian kernel. Provide any number for 'linear'
% kernel, it won't affect the result.
% opts: some other options:
% opts.centerX: center the columns in X. Default is 1(True).
% opts.scaleX: scale the columns in X to unit standard deviation. Default is 0(False).
% opts.scaleY: scale the columns in Y to unit standard deviation. Default is 0(False).

% output
% obj:
% obj.v: the principal components, p by nPC matrix
% obj.Xv: the projected data, i.e. X times v
% ...: input parameters for the function
%%
strcmp(kernel,'linear')

% check whether a whole row in X is missing
Xmis = sum(~isnan(X), 2);
if (sum(Xmis==0))
    error('some samples in X is missing');
end
  ####check whether the number of samples in X and Y match
  if (dim(X)[1]!=dim(Y)[1]){
    stop("The numbers of samples in X ( nrow(X) ) and Y ( nrow(Y) ) do not match")
  }
  ####check whether lambda is non-negative
  if (lambda<0){
    stop("lambda should be non-negative")
  }
