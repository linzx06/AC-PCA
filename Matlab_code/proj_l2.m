function [ obj ] = proj_l2( v, c)
    obj = v/norm(v)*c;
end