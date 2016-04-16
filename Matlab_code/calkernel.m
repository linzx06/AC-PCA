function [ K ] = calkernel( Y, kernel, bandwidth, scaleY)
%% 
% calkernel -- calculate the kernel matrix

% input 
% Y: the n by q confounder matrix, where n is the number of samples, q is 
% the number of confounding factors.  
% Missing values in Y should be labeled as NaN.
% kernel: the kernel to use, should be either 'linear', 'gaussian'.
% bandwidth: bandwidth for gaussian kernel. Provide any number for 'linear'
% kernel, it won't affect the result.
% scaleY: optional. whether columns in Y should have unit standard deviation. Default is 0(False).

% output
% K: n by n kernel matrix
%%

if nargin < 4 
    scaleY = 0; 
end
% scale the Y matrix
if (scaleY)
    Ysd = std(Y,'omitnan');
    Ysd(Ysd==0) = 1;
    Y = Y./repmat(Ysd,nY,1);
end
Y(isnan(Y)) = mean(mean(Y, 'omitnan'));

if (strcmp(kernel,'linear'))
    K = Y*Y';
end
if (strcmp(kernel,'gaussian'))
     K = dist(Y');
     K = exp(-K.^2/2/bandwidth^2);
end

if (strcmp(kernel,'linear') && strcmp(kernel,'gaussian'))
    error('Please select a valid kernel, linear kernel or gaussian kernel')    
end
