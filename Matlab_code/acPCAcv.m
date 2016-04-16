function [ obj ] = acPCAcv(X, Y, lambdas, nPC, kernel, bandwidth, opts)
%% 
% acPCAcv -- performs cross-validation for acPCA to select the lambda
% parameter

% input 
% X: a the n by p data matrix, where n is the number of samples, p is the 
% number of variables
% Missing values in X should be labeled as NaN. If a whole sample in X is
% missing, it should be removed.
% Y: the n by q confounder matrix, where n is the number of samples, q is 
% the number of confounding factors.  
% Missing values in Y should be labeled as NaN.
% lambdas: a vector of tuning parameters, non-negative
% nPC: number of principal components to compute
% kernel: the kernel to use, should be either 'linear', 'gaussian'.
% bandwidth: bandwidth for gaussian kernel. Provide any number for 'linear'
% kernel, it won't affect the result.
% opts: optional, some other options:
% opts.centerX: center the columns in X. Default is 1(True).
% opts.scaleX: scale the columns in X to unit standard deviation. Default is 0(False).
% opts.scaleY: scale the columns in Y to unit standard deviation. Default is 0(False).
% opts.fold: # of folds for cross validation. Default is 10. 
% opts.perc: a cut-off to select the elbow point of the loss. Default is 0.05.  
% opts.plot: 1 or 0. 1 generates the diagnosis plot (lambda vs. loss). 
% Default is 1.
% opts.quiet: 1 or 0. Output the progress of the program. Default is 0.
% Default is 0.

% output
% obj:
% obj.loss: a vector with the loss. Same length as lambdas.
% obj.best_lambda: the best lambda after cross-validation
% ...: input parameters for the function
%%

if nargin < 7
    opts = [];  
    opts.centerX = 1; 
    opts.scaleX = 0;
    opts.scaleY = 0; 
    opts.fold = 10; 
    opts.perc = 0.05;
	opts.plot = 1; 
    opts.quiet = 0;
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

% generate labels for CV
if (nX < opts.fold)
    warning('Cannot tune lambda when fold > sample size, setting fold=sample size');
    foldlab = randsample(nX,nX);
end

if (nX == opts.fold)
    foldlab = randsample(nX,nX);
end

if (nX > opts.fold)
    num = floor(nX/opts.fold);
    foldlab = repmat(1:opts.fold, 1, num)';
    foldlab = [foldlab;randsample(opts.fold, nX-num*opts.fold)];
    foldlab = foldlab(randsample(nX,nX));
end
  
K = calkernel(Y, kernel, bandwidth);

% variable with the largest magnitude of loading, make sure signs match
maxlab = nan(length(lambdas), nPC);
maxlab_sign = maxlab;

optseigs.isreal = 1; 
optseigs.issym = 1;
for lab = 1:length(lambdas)
    lambda = lambdas(lab);
    calAv = @(v) (X' - lambda*X'*K)*(X*v);
    [V,~] = eigs(calAv,p,nPC,'LA', optseigs);
    [~,I] = max(abs(V));
    maxlab(lab,:) = I;
    maxlab_sign(lab,:) = diag( sign(V(I,:)))'; 
end
Xv = nan(length(lambda), nX, nPC);
for skip = 1:opts.fold
    if (~opts.quiet)
        disp(['Running fold ' num2str(skip)]);
    end
    % CV datasets
    Xtmp = X(foldlab~=skip,:);
    Ytmp = Y(foldlab~=skip,:);
    Ktmp = calkernel(Ytmp, kernel, bandwidth);
    for lab = 1:length(lambdas)
        lambda = lambdas(lab);
        calAvtmp = @(v) (Xtmp' - lambda*Xtmp'*Ktmp)*(Xtmp*v);
        [Vtmp,~] = eigs(calAvtmp,p,nPC,'LA', optseigs);
        Vtmp = Vtmp * diag(sign(diag(Vtmp(maxlab(lab,:),:)))) * diag(maxlab_sign(lab,:));
        Xv(lab, foldlab==skip,:) = X(foldlab==skip,:)*Vtmp;
    end
end

% calculate the loss
loss = nan(length(lambdas),1);
for lab = 1:length(lambdas)
    loss(lab) = sum(diag(squeeze(Xv(lab,:,:))'*K*squeeze(Xv(lab,:,:))));
end
thres = (max(loss)-min(loss))*opts.perc + min(loss);
best_lambda = min(lambdas(loss<=thres));

% the lambdas vs. loss plot
if (opts.plot) 
    hax=axes; 
    plot(lambdas, loss, 'bo--');
    hold on 
    line([best_lambda best_lambda],get(hax,'YLim'),'Color',[1 0 0])
    xlabel('lambdas') 
    ylabel('loss') 
end


obj.loss = loss;
obj.best_lambda = best_lambda;
obj.kernel = kernel;
if (strcmp(kernel,'gaussian'))
    obj.bandwidth = bandwidth;
end