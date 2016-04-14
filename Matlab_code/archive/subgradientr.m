function [ obj ] = subgradientr( muX, Y, c1, c2, v_ini, alpha, maxiter, d)
    valA = [];
    %%update v
    v = v_ini;
    iter = 1; 
    v = v - alpha/iter*muX';
    valA = [valA muX*v];
    l2muX = norm(muX);
    seqA = unidrnd(4,1,maxiter);
    while iter<= maxiter
        seq = seqA(iter);
        if seq == 1
            if norm( v, 1 ) > c2;
                v = v - alpha/(iter^d)/c2*sign(v);
            end
        end
        if seq == 2
            if norm( Y*v ) > c1;
                tmp = Y'*(Y*v);
                v = v - alpha/(iter^d)/c1/sqrt(v'*tmp)*tmp;
            end
        end
        if seq == 3
            if norm( v ) > 1;
                v = v - alpha/(iter^d)/norm(v)*v;
            end
        end
        if seq == 4
            v = v - alpha/(iter^d)/l2muX*muX';
            valA = [valA muX*v];
        end   
        iter = iter + 1;
    end
    obj.v = v;
    obj.objec_v = valA;
end