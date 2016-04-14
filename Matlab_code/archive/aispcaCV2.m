function [ obj ] = aispcaCV2( X, Y, c1, c2s, v_ini, mode, opts)
    if nargin < 7
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
    %v = v_ini;
    for c2 = c2s
        display(c2)
        X_hat = zeros(size(X));
        difsqr = [];
        for f = 1:opts.fold
            %display(f);
            X_cv = X; 
            if opts.center==1
                X_cv(lab==f) = NaN; 
                %X_cv = center_na2(X_cv); 
            end
            X_cv(lab==f) = mean(X_cv(lab~=f)); 
            %v = v_ini;
   
            tmp = aispca( X_cv, Y, c1, c2, v_ini, mode);
            v = tmp.v;
            u = tmp.u;
            tmp = u'*X_cv*v*u*v';
            difsqr = [difsqr sum((tmp(lab==f) - X(lab==f)).^2)];
        end
        mse = [mse mean(difsqr)];
        mse_sd = [mse_sd std(difsqr)];
        %v_ini = v;
    end   
    obj.mse = mse;
    obj.mse_sd = mse_sd;
    obj.c2s = c2s;
end