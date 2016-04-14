function [ obj ] = pwdcvxCV( X, Y, c1, c2s, v_ini, opts)
    if nargin < 6
        opts = [];
        opts.maxiter = 25;
        opts.delta = 10^-4;
        opts.fold = 5;
        opts.center = 1;
    end
    vlen = length(v_ini);
    lab = unidrnd(opts.fold, size(X));
    mse = [];
    mse_sd = [];
    for c2 = c2s
        display(c2)
        X_hat = zeros(size(X));
        difsqr = [];
        for f = 1:opts.fold
            %display(f);
            X_cv = X; 
            if opts.center==1
                X_cv(lab==f) = NaN; 
                X_cv = center_na2(X_cv); 
            end
            X_cv(lab==f) = 0; 
            iter = 0;
            converge = 0;
            v = v_ini;
            val = -Inf;
            while converge==0 && iter<= opts.maxiter
                %%update u
                tmp = X_cv*v;
                u = tmp./norm(tmp);
                %%update v
                tmp = -u'*X_cv; 
                cvx_begin quiet
                    variable vnew(vlen)
                    minimize( tmp*vnew )
                    subject to
                    norm( Y*vnew ) <= c1;
                    norm( vnew, 1 ) <= c2;
                    norm( vnew ) <= 1;
                cvx_end
                %%check convergence
                valnew = u'*X_cv*vnew;
                if valnew - val <= opts.delta
                    converge = 1;
                end
                iter = iter + 1;
                v = vnew;
                val = valnew;
                %display(val);
            end
            tmp = u'*X_cv*v*u*v';
            difsqr = [difsqr sum((tmp(lab==f) - X(lab==f)).^2)];
        end
        mse = [mse mean(difsqr)];
        mse_sd = [mse_sd std(difsqr)];
    end   
    obj.mse = mse;
    obj.mse_sd = mse_sd;
    obj.c2s = c2s;
end