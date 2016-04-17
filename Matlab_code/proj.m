function [ obj ] = proj( M, S, V, utx, c1, c2, v_ini, tmin, tmax, iternum)
epsilon = 10^-3/iternum; 
epsilon1 = 10^-3/iternum; 
epsilon2 = 10^-4/iternum; 
maxiter = round(50*sqrt(iternum));

tlow = tmin; tup = tmax; v = v_ini;
while tup - tlow >= tup*epsilon
    t = (tup + tlow)/2;
    iter = 1; converge = 0; 
    while converge==0 && iter<=maxiter
        %%project to l1 ball
        if norm(v, 1) > c2
            v = proj_l1(v, c2);
        end 
        %%project to l2 ball
        if norm(v) > 1
            v = proj_l2( v, 1);
        end 
        %%project to hyperplane
        if t - utx*v > 0
            v = proj_hp( v, utx', t);
        end
        %%project to elliptic
        P = (V'*v).^2; 
        if norm(M*v) > sqrt(c1)      
            l = proj_l2g_getl( S, P, c1);
            v = v - V*(diag(1./(1+ 1./(l*S.^2)))*(V'*v));
        end
        %%check constrain
        if t - utx*v <= t*epsilon && norm(v, 1) <= c2*(1 + epsilon1) && norm(v) <= 1 + epsilon2
           converge=1; 
        end
        iter = iter + 1; 
    end
    if converge == 1
        tlow = utx*v;
        vcon = v;
    else
        tup = t;
    end
end   


obj.v = vcon;
obj.t = utx*vcon;

end