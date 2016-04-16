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
% opts: optional, some other options:
% opts.centerX: center the columns in X. Default is 1(True).
% opts.scaleX: scale the columns in X to unit standard deviation. Default is 0(False).
% opts.scaleY: scale the columns in Y to unit standard deviation. Default is 0(False).

% output
% obj:
% obj.v: the principal components, p by nPC matrix
% obj.Xv: the projected data, i.e. X times v
% ...: input parameters for the function
%%

if nargin < 7
    opts = [];  
    opts.centerX = 1; 
    opts.scaleX = 0;
    opts.scaleY = 0; 
end
% check whether a whole row in X is missing
Xmis = sum(~isnan(X), 2);
if (sum(Xmis==0))
    error('some samples in X is missing');
end
[nX, p] = size(X);
[nY, ~] = size(Y);
% check whether the number of samples in X and Y match
if (nX~=nY)
    error('The numbers of samples in X and Y do not match')
end
% check whether lambda is non-negative
if (lambda<0)
    error('lambda should be non-negative')
end
% center the X matrix
if (opts.centerX)
    X = X-repmat(mean(X,'omitnan'),nX,1);
end
% scale the X matrix
if (opts.scaleX)
    Xsd = std(X,'omitnan');
    Xsd(Xsd==0) = 1;
    X = X./repmat(Xsd,nX,1);
end
% scale the Y matrix
if (opts.scaleY)
    Ysd = std(Y,'omitnan');
    Ysd(Ysd==0) = 1;
    Y = Y./repmat(Ysd,nY,1);
end
% input the missing values in X and Y with the mean
X(isnan(X)) = mean(mean(X, 'omitnan'));
Y(isnan(Y)) = mean(mean(Y, 'omitnan'));
K = calkernel(Y, kernel, bandwidth);
calAv = @(v) (X' - lambda*X'*K)*(X*v);
optseigs.isreal = 1; 
optseigs.issym = 1;
[V,~] = eigs(calAv,p,nPC,'LA', optseigs);
% output
obj.Xv = X*V;
obj.v = V;
obj.lambda = lambda;
obj.kernel = kernel;
if (strcmp(kernel,'gaussian'))
    obj.bandwidth = bandwidth;
end