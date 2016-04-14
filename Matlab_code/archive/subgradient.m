function [ obj ] = subgradient( muX, Y, c1, c2, v_ini, opts)
    if nargin < 6
        opts = [];
        opts.maxiter = 2000;
        %opts.delta = 10^-8;
        opts.alpha = 0.01;
    end
    alpha = opts.alpha;
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
            v = v - alpha/sqrt(iter)*tmp;
            %if sum(isnan(v))>0
                %display('check c2');
            %end
            flag = 0;
        end
        if norm( Y*v ) > c1;
            vo = v;
            tmp = Y'*(Y*v);
            %v = v - alpha/iter*(tmp/norm(tmp));
            v = v - alpha/sqrt(iter)/sqrt(v'*tmp)*tmp;
            %if sum(isnan(v))>0
                %display('check c1');
            %end
            flag = 0;
        end
        if norm( v ) > 1;
            vo = v;
            v = v - alpha/sqrt(iter)/norm(v)*v;
            %if sum(isnan(v))>0
                %display('check 1');
            %end
            flag = 0;
        end
        if flag == 1
            v = v - alpha/sqrt(iter)*muX';
            %display(muX*v);
            valA = [valA muX*v];
        end    
        iter = iter + 1;
    end
    obj.v = v;
    obj.objec_v = valA;
end