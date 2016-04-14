function [ obj ] = center_na2( X )
    obj = X - mean(X(~isnan(X)));
end