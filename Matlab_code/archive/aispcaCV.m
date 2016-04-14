function [ obj ] = aispcaCV( X, Y, c1, c2s, v_ini, opts)
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
    v = v_ini;
    for c2 = c2s
        display(c2)
        X_hat = zeros(size(X));
        for f = 1:opts.fold
            display(f);
            X_cv = X; 
            if opts.center==1
                X_cv(lab==f) = NaN; 
                X_cv = center_na(X_cv); 
            end
            X_cv(lab==f) = 0; 
            %v = v_ini;
   
            tmp = aispca( X_cv, Y, c1, c2, v, 'fast');
            v = tmp.v;
            u = tmp.u;
            tmp = u'*X_cv*v*u*v';
            X_hat(lab==f) = tmp(lab==f);
        end
        mse = [mse mean(mean((X_hat-X).^2))];
        %v_ini = v;
    end   
    obj.mse = mse;
    obj.c2s = c2s;
end