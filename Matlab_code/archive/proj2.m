function [ obj ] = proj2( Y, S, V, utx, c1, c2, v_ini, tmin, tmax, iternum, mode)
if strcmp(mode,'normal')
    epsilon = 10^-3/iternum; epsilon1 = 10^-3/iternum; epsilon2 = 10^-4/iternum; maxiter = round(50*sqrt(iternum));
end
if strcmp(mode,'fast')
    epsilon = 10^-3; epsilon1 = 10^-3; epsilon2 = 10^-4; maxiter = 50;
end
tlow = tmin; tup = tmax; v = v_ini;
while tup - tlow >= tup*epsilon
    t = (tup + tlow)/2;
    iter = 1; converge = 0; 
    l = 0;
    while converge==0 && iter<=maxiter
        %display(iter);
        %%project to l1 ball
        if norm(v, 1) <= c2
        else
            v = proj_l1(v, c2);
        end 
        %%project to l2 ball
        if norm(v) <= 1
        else
            v = proj_l2( v, 1);
        end 
        %%project to hyperplane
        if t - utx*v <= 0
        else
            v = proj_hp( v, utx', t);
        end
        %%project to elliptic
        P = (V'*v).^2; 
        if norm(Y*v)<=c1      
        else    
            l = proj_l2g_getl( S, P, c1);
            %display(l);
            %if l < 0
                %display('error l is negative');
            %end
            %tmp = (eye(p) + l*(Y'*Y));
            %v = tmp\v;
            %v = v - (V*diag(1./(1+ 1./(l*S.^2))))*(V'*v);
            v = v - V*(diag(1./(1+ 1./(l*S.^2)))*(V'*v));
        end
        %%check constrain
        if t - utx*v <= t*epsilon && norm(v, 1) <= c2*(1 + epsilon1) && norm(v) <= 1 + epsilon2
           converge=1; 
        end
        iter = iter + 1; 
        %display([iter l utx*v-t norm(v, 1)-c2 norm(v)-1]);
    end
    if converge == 1
        tlow = utx*v;
        vcon = v;
    else
        tup = t;
    end
    %display([tlow tup]);
end   


obj.v = vcon;
obj.t = utx*vcon;

end