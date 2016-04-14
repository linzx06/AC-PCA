function [ obj ] = v_l2g( S, P, c1, l)
    obj = sum(l*S.^2./(1 + l*S.^2).*P) - l*c1^2;
end