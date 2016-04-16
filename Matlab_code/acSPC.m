function [ obj ] = acSPC( X, Y, X4Y, c1, c2, v_ini, v_substract, kernel, bandwidth, opts)
%% 
% acSPC -- acPCA with sparse loadings

% input 
% X: a the n by p data matrix, where n is the number of samples, p is the 
% number of variables
% Missing values in X should be labeled as NaN. If a whole sample in X is
% missing, it should be removed.
% Y: the n by q confounder matrix, where n is the number of samples, q is 
% the number of confounding factors.  
% Missing values in Y should be labeled as NaN.
% X4Y: 'd' or a matrix of the same dimension as X, this matrix is used to 
% calculate the empirical Hilbert Schmidt criterion. The option 'd' means 
% default and uses X.
% c1: 'd' or a non-negative scalar. Tuning parameter for v'X'KXv. The option 'd' means 
% default and uses v_ini'X4Y'KX4Yv_ini. 
% c2: tuning parameter controlling sparsity, non-negative.
% v_ini: the initial v. Recommended to be the estimate of the non-sparse 
% version, for a "warm" start.
% v_substract: 'd' or a p by k matrix, the principal components 
% to be subtracted, where k is the number of PCs to be substracted. The 
% option 'd' means default and does not substract anything, i.e. calculate 
% the first principal component.
% kernel: the kernel to use, should be either 'linear', 'gaussian'.
% bandwidth: bandwidth for gaussian kernel. Provide any number for 'linear'
% kernel, it won't affect the result.
% opts: optional, some other options:
% opts.centerX: center the columns in X. Default is 1(True).
% opts.scaleX: scale the columns in X to unit standard deviation. Default is 0(False).
% opts.scaleY: scale the columns in Y to unit standard deviation. Default is 0(False).
% opts...: other parameters for convergence.

% output
% obj:
% obj.v: the sparse principal component.
% obj.u: the u vector.
% obj.converge: whether the algorithm converged.
%%

if nargin < 10
    opts = [];  
    opts.centerX = 1; 
    opts.scaleX = 0;
    opts.scaleY = 0; 
    opts.maxiter = 50;
    opts.delta = 10^-8;
    opts.filter = 1;
end

% check whether a whole row in X is missing
Xmis = sum(~isnan(X), 2);
if (sum(Xmis==0))
    error('some samples in X is missing');
end
[nX, p] = size(X);
[nY, q] = size(Y);

% check whether the number of samples in X and Y match
if (nX~=nY)
    error('The numbers of samples in X and Y do not match')
end
if (~strcmp(X4Y,'d'))
    % check whether the number of samples in X4Y and Y match
    [nX4Y, ~] = size(X4Y);
    if (nX4Y~=nY)
        error('The numbers of samples in X4Y and Y do not match')
    end
end
[nY, ~] = size(X4Y);
% check whether lambdas is non-negative
if (sum(lambdas<0))
    error('lambdas should be non-negative')
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
if (strcmp(c1,'d'))
    if (strcmp(X4Y,'d'))
        c1 = v_ini'*X'*K*X*v_ini;
    else
        X4Y(isnan(X4Y)) = mean(mean(X4Y, 'omitnan'));
        c1 = v_ini'*X4Y'*K*X4Y*v_ini;
    end
end

% Decomposition on the kernel matrix to get the M matrix
if (strcmp(kernel,'linear') && q<=nY)
    if (strcmp(X4Y,'d'))
      M = Y'*X;  
    else
      M = Y'*X4Y;  
    end
else
    [U,S,~] = svd(K, 'econ');
    S = diag(S);
    Delta = U(:,S>0)*diag(sqrt(S(S>0)));
    if (strcmp(X4Y,'d'))
      M = Delta'*X;  
    else
      M = Delta'*X4Y; 
    end
end

% svd on M
[~,S,V] = svd(M, 'econ');

% substract the first several PCs (i.e. v_substract), if v_substract is provided
if (strcmp(v_substract,'d'))
    % substract the PCs in v_substract one by one
    [~, k] = size(v_substract);
    for pc = 1:k
      X = X - (X*v_substract(:,pc))*v_substract(:,pc)';   
    end
end
