function w = proj_l1(v, c2)
% PROJECTONTOL1BALL Projects point onto L1 ball of specified radius.
%
% w = ProjectOntoL1Ball(v, b) returns the vector w which is the solution
%   to the following constrained minimization problem:
%
%    min   ||w - v||_2
%    s.t.  ||w||_1 <= b.
%
%   That is, performs Euclidean projection of v to the 1-norm ball of radius
%   b.
%
% Author: John Duchi (jduchi@cs.berkeley.edu)

if (c2 < 0)
  error('Radius of L1 ball is negative: %2.3f\n', c2);
end
if (norm(v, 1) < c2)
  w = v;
  return;
end
u = sort(abs(v),'descend');
sv = cumsum(u);
rho = find(u > (sv - c2) ./ (1:length(u))', 1, 'last');
theta = max(0, (sv(rho) - c2) / rho);
w = sign(v) .* max(abs(v) - theta, 0);