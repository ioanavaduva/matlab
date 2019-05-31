%  Code to run kpik.m for Poisson with 63 unknowns
n = 500;
h = 1/n;
T = (-(diag(2*ones(n, 1)) + diag (-1*ones(n-1, 1), 1) + diag (-1*ones(n-1, 1), -1)))/h^2;
E = eye(n);
LE = chol(E,'lower');

%x = linspace(0,1,n);
%y = linspace(1, 0, n);
%b = @(x, y) sin(pi.*x).*cos(pi.*y);
%B = b(x, y)';
B = ones(n, 1);

m = 100;
tol = 1e-9;
tolY = 1e-12;

[Z,r]=kpik(T,E,LE,B,m,tol,tolY);

fprintf('final true absolute residual norm: \n')
disp(norm(T*Z*Z'*E+E*Z*Z'*T'+B*B'))    %this matrix should never be formed for n large 