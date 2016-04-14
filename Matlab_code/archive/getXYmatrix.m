function [ obj ] = getXYmatrix( data )
%%
[nid, ~, ng] = size(data);
X = [];
Y = [];
%X
for i = 1:nid
    dtmp = data(i,:,:);
    dtmp = squeeze(dtmp);
    dtmp = dtmp(~any(isnan(dtmp),2),:);
    X = [X; dtmp];
end
%Y
for l = 1:(nid-1)
    for k = (l+1):nid
        ddtmp = data(l,:,:) - data(k,:,:);
        ddtmp = squeeze(ddtmp);
        ddtmp = ddtmp(~any(isnan(ddtmp),2),:);
        Y = [Y; ddtmp];    
    end
end

obj.X = X;
obj.Y = Y;

end