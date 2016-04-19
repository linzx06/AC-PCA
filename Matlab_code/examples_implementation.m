%%%%%%%%%%%%%%%%%%%%%% Implementation examples %%%%%%%%%%%%%%%%%%%%%
% More details for the input and output of the functions are included the 
% function files ('.m')

%% Implementation of acPCAcv and acPCA 
dataA = load('data_example1.mat');
% input for acPCAcv and acPCA
X = dataA.X;
Y = dataA.Y;
lambdas = 0:0.05:2;
nPC = 2; 
kernel = 'linear';
bandwidth = 1; %bandwidth for gaussian kernel. Provide any number for 'linear'
% kernel, it won't affect the result.

result1cv = acPCAcv(X, Y, lambdas, nPC, kernel, bandwidth);
lambda = result1cv.best_lambda;
result1 = acPCA(X, Y, lambda, nPC, kernel, bandwidth);
result1.Xv % the principal components, p by nPC matrix
result1.v % the projected data, i.e. X times v
%%

%% Implementation of acSPCcv and acSPC, we use a warm start from acPCA. PC1 
dataA = load('data_example5.mat');
% input for acPCAcv and acPCA
X = dataA.X;
Y = dataA.Y;
lambdas = 0:0.05:2;
nPC = 1; 
kernel = 'linear';
bandwidth = 1; %bandwidth for gaussian kernel. Provide any number for 'linear'
% kernel, it won't affect the result.

result2cv = acPCAcv(X, Y, lambdas, nPC, kernel, bandwidth);
result2 = acPCA(X, Y, result2cv.best_lambda, nPC, kernel, bandwidth);

% input for acSPCcv and acSPC
v_ini = result2.v;
v_substract = 'd';
X4Y = 'd';
c1 = 'd';
kernel = 'linear'; 
bandwidth = 1; 

% a coarse search first
c2s_coarse = (1:(-0.1):0)*sum(abs(v_ini));
result3cv_coarse = acSPCcv( X, Y, X4Y, c1, c2s_coarse, v_ini, v_substract, kernel, bandwidth);

% a finer search 
c2s_fine = (0.7:(-0.02):0.4)*sum(abs(v_ini));
result3cv_fine = acSPCcv( X, Y, X4Y, c1, c2s_fine, v_ini, v_substract, kernel, bandwidth);

% use the best c2 in the finer search
c2 = result3cv_fine.best_c2;
result3 = acSPC( X, Y, X4Y, c1, c2, v_ini, v_substract, kernel, bandwidth);
result3.v
%%

%% Implementation of acSPCcv and acSPC, we use a warm start from acPCA. PC2
dataA = load('data_brain_w2.mat');
X = dataA.X;
Y = dataA.Y;
lambdas = 0:0.05:2;
nPC = 2; 
kernel = 'linear';
bandwidth = 1; %bandwidth for gaussian kernel. Provide any number for 'linear'
result4cv = acPCAcv(X, Y, lambdas, nPC, kernel, bandwidth);
result4 = acPCA(X, Y, result4cv.best_lambda, nPC, kernel, bandwidth);

% input for acSPCcv and acSPC
v_ini = result4.v(:,2);
v_substract = result4.v(:,1);
X4Y = 'd';
c1 = 'd';
kernel = 'linear'; 
bandwidth = 1; 

c2s_coarse = (1:(-0.1):0)*sum(abs(v_ini));
result5cv_coarse = acSPCcv( X, Y, X4Y, c1, c2s_coarse, v_ini, v_substract, kernel, bandwidth); 
%%