function [ obj ] =  topthresfilter(v, top, alpha) 
vsort = sort(abs(v), 'descend');
thres = mean( vsort(1:ceil(top*length(v))) ) * alpha;
v( abs(v) <= thres ) = 0;
obj = v;