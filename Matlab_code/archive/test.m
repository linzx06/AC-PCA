cd ~/Dropbox/brain_gradient/code/matlab_code/

Qm =  calQmatrix( data, 1 );
test = eigs(Qm, 2);
[V,D] = eigs(Qm);

XY = getXYmatrix(data);
Xm = XY.X;
Ym = XY.Y;

c1 = norm(Ym * V(:,1));
%%test cvx
c2 = c1;
v_ini = V(:,1);
[sv iv] = sort(abs(v_ini), 'descend');

test = pwdcvx( X, Y, c1, c2, v_ini);
v_sp = test.v;
%the non-zero elements
sum(abs(v_sp)>=10^-3)

%a = 1 - 2.*rand(100,1);
%norm(a, 1)
%norm(a, 2)
%norm(a, 3)

%%the trace of soft thresholding
ls = [0.1:0.1:1 1.5:0.5:5 10];
v_tr = [];
cvs = [];
for l=ls
    c2 = c1*l;
    test = pwdcvx( X, Y, c1, c2, v_ini);
    v_tr = [v_tr  test.v];
    cvs = [cvs test.converge];
end

v_trnew = [];
for i = 1:19
    lab1 = (i-1)*5000 + 1;
    lab2 = lab1 + 4999;
    v_trnew = [v_trnew v_tr(lab1:lab2)];
end
plot(ls, v_trnew([1:10 iv(1:10)'],:))

plot(ls, v_trnew(iv(1:10)',:))

plot(ls, v_trnew(iv(1:20)',:))

[tmp iv5] = sort(abs(v5), 'descend');
v5 = v_trnew(:, 5);
[tmp iv5] = sort(abs(v5), 'descend');

c2s = c1*[0.2:0.2:1 1.5 2 5 10];
test_cv = pwdcvxCV( X, Y, c1, c2s, v_ini);

c2s = c1*[0.1 0.4 1 1.5 2 5 10];
test_cv = pwdcvxCV( X, Y, c1, c2s, v_ini);

c2s = 10000;
test_cv_l = pwdcvxCV( X, Y, c1, c2s, v_ini);

c2s = 23:2:45;
test_cv_f = pwdcvxCV( X, Y, c1, c2s, v_ini);