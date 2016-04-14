
tmp = -u'*X; 
tic;
        cvx_begin quiet
            variable vnew(vlen)
            minimize( tmp*vnew )
            subject to
            norm( Y*vnew ) <= c1;
            norm( vnew, 1 ) <= c2;
            norm( vnew ) <= 1;
        cvx_end
        toc;
tic;   
test = subgradient( tmp, Y, c1, c2, v_ini);
toc;

tic;   
test = subgradient( tmp, Y, c1, c2, v_ini);
toc;

tic;   
test = subgradient1( tmp, Y, c1, c2, v_ini);
toc;

tic;   
test = subgradientn( tmp, Y, c1, c2, v_ini);
toc;

tic;   
test = subgradientn( tmp, Y, c1, c2, v_ini, 0.08);
toc;

tic;   
test = subgradientr( tmp, Y, c1, c2, v_ini, 0.1, 20000, 0.3);
toc;
test.objec_v

tic;   
test = subgradientn( tmp, Y, c1, c2, v_ini, 0.1, 10000);
toc;

lab = 10000:26000;
plot(lab, test.objec_trace(lab))

tic;   
test = subgradientn( tmp, Y, c1, c2, test.v, 0.001, 10000);
toc;
plot(length(test.objec_trace), test.objec_trace);

tic;   
test = subgradientn( tmp, Y, c1, c2, v_ini, 0.1, 1000);
toc;

tic;   
test = subgradientnad( tmp, Y, c1, c2, v_ini, 0.1, 10000);
toc;

tic;   
test = subgradientn( tmp, Y, c1, c2, v_ini, 0.1, 2000);
toc;

c1 = 10;
tic;   
test10 = subgradientn( tmp, Y, c1, c2, v_ini, 0.1, 2000);
toc;

tic;
        cvx_begin quiet
            variable vnew(vlen)
            minimize( tmp*vnew )
            subject to
            norm( Y*vnew ) <= c1;
            norm( vnew, 1 ) <= c2;
            norm( vnew ) <= 1;
        cvx_end
        toc;
tic;


c1 = c2;
tic;
test1 = pwdcvx( X, Y, c1, c2, v_ini);
toc;
tic;
test2 = pwdsubgrd( X, Y, c1, c2, v_ini);
toc;