function [ obj ] = calQmatrix( data, lambda )
%%
[nid, nbr, ng] = size(data);
Qmatrix = zeros(ng);
%individual contribution
for i = 1:nid
    dtmp = data(i,:,:);
    dtmp = squeeze(dtmp);
    dtmp = dtmp(~any(isnan(dtmp),2),:);
    Qmatrix = Qmatrix + dtmp'*dtmp;
end
%cross-individual contribution
for br = 1:nbr
    dtmp = data(:,br,:);
    dtmp = squeeze(dtmp);
    dtmp = dtmp(~any(isnan(dtmp),2),:);
    [nid, ~] = size(dtmp);
    if nid<=1
      
    else 
        for l = 1:(nid-1)
            for k = (l+1):nid
                ddtmp = dtmp(l,:) - dtmp(k,:);
                Qmatrix = Qmatrix - lambda*(ddtmp'*ddtmp);        
            end
        end
    end
end  

obj = Qmatrix;

end


