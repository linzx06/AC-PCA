function [ obj ] = pwdsubgrd( X, Y, c1, c2, v_ini, opts)
    if nargin < 6
        opts = [];
        opts.maxiter = 20;
        opts.delta = 10^-8;
        opts.center = 1;
        opts.alpha = 0.1;
    end
    %%center X
    if opts.center==1
        X = X - repmat(mean(X),size(X,1),1);
    end
    v = v_ini;
    vlen = length(v);
    devia = [];
    valA = [];
    iter = 1; converge = 0;
    val = -Inf;
    while converge==0 && iter<= opts.maxiter
        display(iter)
        %%update u
        tmp = X*v;
        u = tmp./norm(tmp);
        %%update v
        tmp = -u'*X;     
        subtmp = subgradientn( tmp, Y, c1, c2, v, opts.alpha/sqrt(iter), 20000);
        vnew = subtmp.v;
        %%check convergence
        devia = [devia, norm(v-vnew)];
        valnew = u'*X*vnew;
        if valnew - val <= opts.delta && valnew - val >= 0
            converge = 1;
        end
        iter = iter + 1;
        v = vnew;
        val = valnew;
        valA = [valA val];
        display(val);
    end
    obj.v = v;
    obj.devia = devia;
    obj.objec_v = valA;
    obj.converge = converge;
end