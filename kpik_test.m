%  Code to run kpik.m for Poisson with 63 unknowns
nh = 100;
A = -kr_pois(nh); 
n = nh^2;
E = eye(n);
LE = chol(E,'lower');

xtemp = linspace(0,1,nh);
x = repmat(xtemp, 1, nh);
y = reshape(repmat(xtemp, length(xtemp), 1), 1, length(xtemp)^2);
b = @(x, y) sin(pi.*x).*cos(pi.*y);
B = b(x, y)';

m = 100;
tol = 1e-9;
tolY = 1e-12;

[Z,r]=kpik(A,E,LE,B,m,tol,tolY);

fprintf('final true absolute residual norm: \n')
disp(norm(A*Z*Z'*E+E*Z*Z'*A'+B*B'))    %this matrix should never be formed for n large 