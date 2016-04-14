cd ~/Dropbox/brain_gradient/code/matlab_code/ 
%%test each projection function
v = rand(1000,1);
norm(v,1)
norm(v,2)

vl1 = proj_l1(v, 100);
vl2 = proj_l2(v, 1);

Y = rand(100, 1000);
utx = rand(1,1000);
t = utx*v;
[U,S,V] = svd(Y, 'econ');
S = diag(S);

P = (V'*v).^2; P = P(1:100);
c1 = 1000;   
l = proj_l2g_getl( S, P, c1);

vs = [];
ls = 1:1:1000;
ls = 0:0.001:1;
for l = ls
    vs = [vs v_l2g( S, P, c1, l)];
end
plot(ls, vs)

vhp = proj_hp( v, utx', t/2);

%%test when Y is identity
Y = eye(1000);
[U,S,V] = svd(Y);
S = diag(S);

P = (V'*v).^2; 
c1 = 10;   
l = proj_l2g_getl( S, P, c1);

%%check on the real data
Qm =  calQmatrix( data,  2);
%test = eigs(Qm, 2);
[V,D] = eig(Qm);
[D,I] = sort(diag(D));
V = V(:, I);
d = size(data, 3);
V = V(:,d:(-1):(d-1));

XY = getXYmatrix(data);
Xm = XY.X;
Ym = XY.Y;

c1 = norm(Ym * V(:,1));
%%test cvx
c2 = c1;
v_ini = V(:,1);
utx = v_ini'*Xm'*Xm;


[U1,S1,V1] = svd(Ym);
[U,S,V] = svd(Ym, 'econ');
S = diag(S);

tmin = 0;
tmax = norm(utx);
test = proj( Ym, S, V, utx, c1, c2, v_ini, tmin, tmax);
tic;
test = proj2( Ym, S, V, utx, c1, c2, v_ini, tmin, tmax);
toc;

tic;
cvx_begin quiet
            variable vnew(vlen)
            minimize( -utx*vnew )
            subject to
            norm( Ym*vnew ) <= c1;
            norm( vnew, 1 ) <= c2;
            norm( vnew ) <= 1;
cvx_end
toc;

%%%vary c2
tic;
cvx_begin quiet
            variable vnew(vlen)
            minimize( -utx*vnew )
            subject to
            norm( Ym*vnew ) <= c1;
            norm( vnew, 1 ) <= c2/2;
            norm( vnew ) <= 1;
cvx_end
toc;

tic;
cvx_begin quiet
            variable vnew(vlen)
            minimize( -utx*vnew )
            subject to
            norm( Ym*vnew ) <= c1;
            norm( vnew, 1 ) <= c2/4;
            norm( vnew ) <= 1;
cvx_end
toc;

%%test the whole function
tic;
test1 = aispca( Xm, Ym, c1, c2/2, v_ini);
toc;

tic;
test2 = pwdcvx( Xm, Ym, c1, c2/2, v_ini);
toc;

tic;
test3 = aispca( Xm, Ym, c1, c2/2, v_ini);
toc;

%%check trace
c2s = norm(v_ini, 1)*[0.95:-0.05:0.05];
v = v_ini;
vA = [];
for c2 = c2s
    display(c2);
    tmp = aispca( Xm, Ym, c1, c2, v);
    v = tmp.v;
    vA = [vA, v];
end

tt = find(vA(:,19)~=0);
plot(c2s, vA(tt,:))


test1 = pwdcvx( Xm, Ym, norm(Ym*v,2), norm(v,1), v);
tic;
test2 = pwdcvx( Xm, Ym, norm(Ym*v,2), norm(v,1), v_ini);
toc;
v = v_ini;
tic;
tmp = aispca( Xm, Ym, c1, c2, v);
toc;
tic;
v = v_ini;
test1 = pwdcvx( Xm, Ym, c1, c2, v);
toc;

tic;
tmpCV = aispcaCV( Xm, Ym, c1, norm(v_ini, 1).*[0.9:-0.05:0.05], v_ini);
toc;

tic;
tmpCVf = aispcaCV( Xm, Ym, c1, norm(v_ini, 1).*[0.9:-0.05:0.05], v_ini);
toc;

tic;
tmpCVf2 = aispcaCV( Xm, Ym, c1, norm(v_ini, 1).*[0.9:-0.05:0.05], v_ini);
toc;

tic;
tmpCVf2 = aispcaCV2( Xm, Ym, c1, norm(v_ini, 1).*[0.9:-0.05:0.05], v_ini);
toc;

tic;
tmpCVf2 = aispcaCV2( Xm, Ym.*0, c1, norm(v_ini, 1).*[0.9:-0.05:0.05], v_ini);
toc;

tic;
tmpCVf3 = aispcaCV( Xm, Ym, c1, norm(v_ini, 1).*[0.98:-0.02:0.02], v_ini);
toc;

tic;
tmpCVcvx= pwdcvxCV( Xm, Ym, c1, norm(v_ini, 1).*[0.9:-0.1:0.1], v_ini);
toc;

%%%check on variance explained
%total variation
tot_var = 0;
for i = 1:6
    dtmp = squeeze(data(i,:,:));
    dtmp = dtmp(~isnan(dtmp(:,1)),:);
    display(size(dtmp))
    tot_var = tot_var + sum(var(dtmp));
end

v = tmp1.v;
aivar = 0;
for i = 1:6
    dtmp = squeeze(data(i,:,:));
    dtmp = dtmp(~isnan(dtmp(:,1)),:);
    tmp = dtmp*v;
    aivar = aivar + var(tmp);
end

aispca( Xm, Ym, c1, c2, v, 'fast');

v = v_ini;
varAI = [];
for c2 = norm(v_ini, 1).*[0.9:-0.05:0.05]
    display(c2);
    tmp = aispca( Xm, Ym, c1, c2, v_ini, 'fast');
    v = tmp.v;
    aivar = 0;
    for i = 1:6
        dtmp = squeeze(data(i,:,:));
        dtmp = dtmp(~isnan(dtmp(:,1)),:);
        tmp = dtmp*v;
        aivar = aivar + var(tmp);
    end
    varAI = [varAI aivar];
end

vA = [];
for c2 = norm(v_ini, 1).*[0.9:-0.05:0.05]
    display(c2);
    tmp = aispca( Xm, Ym, c1, c2, v_ini, 'fast');
    vA = [vA tmp.v];
end
[v_ini, D] = eigs(datatest'*datatest, 1);
re1 = aispca( datatest, datatest.*0, 30, 30, -v_ini, 'normal');
vre1 = re1.v;
re2 = pwdcvx( datatest, datatest.*0, 30, 30, -v_ini);
vre2 = re2.v;

tic;
tmpCVcvx= pwdcvxCV( datatest, datatest.*0, 30, [10 20 30 35 45], v_ini);
toc;
tic;
tmpCVf3 = aispcaCV2( datatest, datatest.*0, 30, [10 20 30 35 45], v_ini, 'fast');
toc;
tic;
tmpCVf4 = aispcaCV2( datatest, datatest.*0, 30, [10 20 30 35 45], v_ini, 'normal');
toc;


[v_ini, D] = eigs(Xm'*Xm, 1);
tic;
tmpCVcvx= pwdcvxCV( Xm, Xm.*0, 30, [10 20 30 35 45], v_ini);
toc;

tic;
tmpCVf3 = aispcaCV2( Xm, Xm.*0, 30, [10 20 30 35 45], v_ini, 'normal');
toc;

tic;
tmpCVf4 = aispcaCV2( Xm, Xm.*0, 30, [10 20 30 35 45], v_ini, 'fast');
toc;

tic;
tmpCVcvx= pwdcvxCV( Xm, Ym, c1, norm(v_ini, 1).*[0.95:-0.1:0.15], v_ini);
toc;

tic;
tmpCVf5 = aispcaCV2( Xm, Ym, c1, norm(v_ini, 1).*[0.95:-0.1:0.15], v_ini, 'fast');
toc;

tic;
tmpCVf6 = aispcaCV2( Xm, Ym, c1, norm(v_ini, 1).*[0.95:-0.1:0.15], v_ini, 'normal');
toc;

tic;
tmpCVf7 = aispcaCV3( Xm, Ym, c1, norm(v_ini, 1).*[0.95:-0.1:0.15], v_ini, 'fast');
toc;

tic;
tmpCVf8 = aispcaCV3( Xm, Ym, c1, norm(v_ini, 1).*[0.95:-0.1:0.15], v_ini, 'normal');
toc;

tic;
tmpCVf9 = aispcaCV4( Xm, Ym, c1, norm(v_ini, 1).*[0.95:-0.1:0.15], v_ini, 'fast');
toc;

tic;
tmpCVf10 = aispcaCV4( Xm, Ym, c1, norm(v_ini, 1).*[0.95:-0.1:0.15], v_ini, 'normal');
toc;

result1 = aispca( Xm, Ym, c1, 6.3904 , v_ini, 'normal'); vest1 = result1.v;
result2 = aispca( Xm, Ym, c1, 3.7, v_ini, 'normal'); vest2 = result2.v;


result = aispca( Xm, Ym, c1, 7, v_ini, 'normal'); vest = result.v; sum(abs(vest)>=10^-5)

tic;
result2 = aispca( Xm, Ym, c1, 5, v_ini, 'normal');
toc;

