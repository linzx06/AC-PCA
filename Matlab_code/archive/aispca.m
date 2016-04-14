function [ obj ] = aispca( X, Y, c1, c2, v_ini, mode, opts)
    if nargin < 7
        opts = [];
        opts.maxiter = 50;
        opts.delta = 10^-8;
        opts.center = 1;
    end
    %%center X
    if opts.center==1
        X = X - repmat(mean(X),size(X,1),1);
    end
    %%svd on Y
    [~,S,V] = svd(Y, 'econ');
    S = diag(S);
    %%initialization
    v = v_ini;
    valA = [];
    devia = [];
    iter = 1; converge = 0;
    val = -Inf;
    %%iteration
    while converge==0 && iter<= opts.maxiter
        %%update u
        tmp = X*v;
        u = tmp./norm(tmp);
        %%update v
        utx = u'*X;
        tmin = 0; tmax = norm(utx);
        tmpproj = proj2( Y, S, V, utx, c1, c2, v, tmin, tmax, iter, mode);
        vnew = tmpproj.v;
        %%check convergence
        devia = [devia, norm(v-vnew)];
        valnew = tmpproj.t;
        if valnew - val <= opts.delta
            converge = 1;
        end
        iter = iter + 1;
        v = vnew;
        val = valnew;
        valA = [valA val];
    end
    obj.v = v;
    obj.u = u;
    obj.objec_v = valA;
    obj.devia = devia;
    obj.converge = converge;
end