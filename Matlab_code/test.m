%%%test the matlab functions
cd /Users/zhixianglin/AC-PCA/Matlab_code/

%% test eigs function
% test on random data
Xtmp = rand(5,10);
Ytmp = rand(5,2);
Ktmp = calkernel( Ytmp, 'linear', 1);

lambda = 1;
calAvtmp = @(v) (Xtmp' - lambda*Xtmp'*Ktmp)*(Xtmp*v);
opts.isreal = 1; 
opts.issym = 1;
[Vtmp,Dtmp] = eigs(calAvtmp,10,2,'LA', opts);
diag(Dtmp)
Atmp = (Xtmp' - lambda*Xtmp'*Ktmp)*Xtmp;
Atmp = (Atmp+Atmp')/2;
eigs(Atmp, 2, 'LA')

% test on simulated data
K = calkernel( Y, 'linear', 1);
[~,p] = size(X);
lambda = 1;
calAv = @(v) (X' - lambda*X'*K)*(X*v);
opts.isreal = 1; 
opts.issym = 1;
[Vtmp,Dtmp] = eigs(calAv,p,2,'LA', opts);
diag(Dtmp)
Atmp = (X' - lambda*X'*K)*X;
Atmp = (Atmp+Atmp')/2;
eigs(Atmp, 2, 'LA')
%% test acPCAcv and acPCA
tic;
resultcv = acPCAcv(X, Y,  0:0.1:2, 2, 'linear', 1);
toc;
result = acPCA(X, Y, resultcv.best_lambda, 2, 'linear', 1);
Xv = result.Xv;
c = [repmat(2,12,1);repmat(3,24,1)];
scatter(Xv(:,1), Xv(:,2),[],c,'filled');
%% test svd
K = calkernel( Y, 'linear', 1);
[U,S,~] = svd(K, 'econ');
S = diag(S);
Delta = U(:,S>0)*diag(sqrt(S(S>0)));

M = Delta'*X; 
tmp1 = M'*M;
tmp2 = X'*K*X;