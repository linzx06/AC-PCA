function [ obj ] = center_na( X )
    [nrow,~] = size(X);
    lab = isnan(X);
    numnna = nrow - sum(lab);
    X(lab==1) = 0;
    X = X - repmat(sum(X)./numnna, [nrow, 1]);
    X(lab==1) = NaN;
    obj = X;
end