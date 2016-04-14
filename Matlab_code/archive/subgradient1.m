function [ obj ] = subgradient1( muX, Y, c1, c2, v_ini, alpha, opts)
    if nargin < 7
        opts = [];
        opts.maxiter = 1000;
    end
    valA = [];
    %%update v
    v = v_ini;
    iter = 1; 
    v = v - alpha/iter*muX';
    valA = [valA muX*v];
    while iter<= opts.maxiter
        flag = 1;
        %%check constrain
        if norm( v, 1 ) > c2;
            vo = v;
            tmp = sign(v);
            v = v - alpha/iter*tmp;
            v = v - alpha/iter*muX';
            valA = [valA muX*v];
            flag = 0;
        end
        if norm( Y*v ) > c1;
            vo = v;
            tmp = Y'*(Y*v);
            v = v - alpha/iter/sqrt(v'*tmp)*tmp;
            v = v - alpha/iter*muX';
            valA = [valA muX*v];
            flag = 0;
        end
        if norm( v ) > 1;
            vo = v;
            v = v - alpha/iter/norm(v)*v;
            v = v - alpha/iter*muX';
            valA = [valA muX*v];
            flag = 0;
        end
        if flag == 1
            v = v - alpha/iter*muX';
            %display(muX*v);
            valA = [valA muX*v];
        end    
        iter = iter + 1;
    end
    obj.v = v;
    obj.objec_v = valA;
end