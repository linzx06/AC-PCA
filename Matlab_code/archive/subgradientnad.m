function [ obj ] = subgradientnad( muX, Y, c1, c2, v_ini, alpha, maxiter)
    valA = [];
    %%update v
    v = v_ini;
    iter = 1; 
    v = v - alpha/iter*muX';
    valA = [valA muX*v];
    l2muX = norm(muX);
    while iter<= maxiter
        flag = 1;
        %%check constrain
        if norm( v, 1 ) > c2
            tmp = sign(v);
            v = v - alpha/sqrt(iter)/c2*tmp;
            flag = 0;
        end
        if norm( v ) > 1 && flag==1
            v = v - alpha/sqrt(iter)*v;
            flag = 0;
        end
        if norm( Y*v ) > c1
            v = v - alpha/sqrt(iter)/(c1^2)*Y'*(Y*v);
            flag = 0;
        end
        if flag == 1
            v = v - alpha/sqrt(iter)/l2muX*muX';
            valA = [valA muX*v];
        end    
        iter = iter + 1;
        if mod(iter, 2000)==0
            alpha = alpha*1;
        end
    end
    obj.v = v;
    obj.objec_trace = valA;
    obj.objec = muX*v;
end