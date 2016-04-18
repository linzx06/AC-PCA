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

%% test acSPCcv and acSPC
cd /Users/zhixianglin/AC-PCA/Matlab_code/
result1cv = acPCAcv(X, Y,  0:0.05:2, 2, 'linear', 1);
result1cv.best_lambda
result1 = acPCA(X, Y, result1cv.best_lambda, 2, 'linear', 1);
v_ini = result1.v;
v_ini = v_ini(:,1);
result2 = acSPC( X, Y, 'd', 'd', 0.5*sum(abs(v_ini)), v_ini, 'd', 'linear', 1);
tic;
result2cv = acSPCcv( X, Y, 'd', 'd', (1:(-0.1):0)*sum(abs(v_ini)), v_ini, 'd', 'linear', 1);
toc;
