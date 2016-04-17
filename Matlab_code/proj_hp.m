function [ obj ] = proj_hp( v, w, t)
    obj = v - ((w'*v - t)/norm(w)^2)*w;
end