function [ obj ] = proj_l2g_getl( S, P, c1)
%S: a vector with singular values, P is a vector with V'v first few
%elements, same length as S
    S2 = S.^2;
    epsilon = 10^-10; maxiter = 50;
    l = 10; lnew = 0; %lt = []; it is important that the initial value for l is 0, otherwise it may go negative and get messy
    iter = 1;
    while abs(l - lnew)>=epsilon && iter <= maxiter
        l = lnew;
        d1 = sum(P.*S2./((1 + l*S2).^2)) - c1^2;
        d2 = -2*sum( P.*S.^4./ ( (1 + l*S2 ).^3 ) );
        lnew = l - d1/d2;
        %lt = [lt lnew];
    end
    obj = lnew;
    %obj.lnew = lnew;
    %obj.lt = lt;
    %obj.iter = lnew;
end