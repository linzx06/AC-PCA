function [ obj ] = proj( Y, S, V, utx, c1, c2, v_ini, tmin, tmax)
epsilon = 10^-2; maxiter = 500;
tlow = tmin; tup = tmax; v = v_ini;
[~, p] = size(Y);
while tup - tlow >= epsilon
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
        if abs(utx*v-t) <= epsilon
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
            v = (eye(p) + l*(Y'*Y))\v;
        end
        %%check constrain
        if abs(utx*v-t) <= epsilon && norm(v, 1) <= c2 && norm(v) <= 1
           converge=1; 
        end
        iter = iter + 1; 
        display([iter l utx*v-t norm(v, 1)-c2 norm(v)-1]);
    end
    if converge == 1
        tlow = t;
    else
        tup = t;
    end
    display([tlow tup]);
end   


obj.v = v;
obj.t = t;

end