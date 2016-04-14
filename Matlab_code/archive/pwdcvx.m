function [ obj ] = pwdcvx( X, Y, c1, c2, v_ini, opts)
    if nargin < 6
        opts = [];
        opts.maxiter = 50;
        opts.delta = 10^-8;
        opts.center = 1;
    end
    %%center X
    if opts.center==1
        X = X - repmat(mean(X),size(X,1),1);
    end
    v = v_ini;
    vlen = length(v);
    devia = [];
    valA = [];
    iter = 0; converge = 0;
    val = -Inf;
    while converge==0 && iter<= opts.maxiter
        %%update u
        tmp = X*v;
        u = tmp./norm(tmp);
        %%update v
        tmp = -u'*X; 
        cvx_begin quiet
            variable vnew(vlen)
            minimize( tmp*vnew )
            subject to
            norm( Y*vnew ) <= c1;
            norm( vnew, 1 ) <= c2;
            norm( vnew ) <= 1;
        cvx_end
        %%check convergence
        devia = [devia, norm(v-vnew)];
        valnew = u'*X*vnew;
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
    obj.devia = devia;
    obj.objec_v = valA;
    obj.converge = converge;
end