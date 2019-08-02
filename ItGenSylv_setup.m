% Set-up for running ItGenSylv.m code. 
% We want to solve a problem of the form 
%       AX + XB + M1XN1 + M2XN2 + ... = C, 
% with N = [N1, N2, ...] and M = [M1, M2, ...].
% Here, A and B are the matrices coming from the finite differences
% discretisation of the Poisson equation

n = 2000; 
h = 1/n;

% A and B
eps = 0.0083;
A = eps*(diag(2*ones(n, 1)) + diag (-1*ones(n-1, 1), 1) + diag (-1*ones(n-1, 1), -1))/h^2;
B = A;

% for N and M
x = linspace(0,1,n);
y = x;
ph1 = @(x) 1-x.^2;
ph2 = @(x) -2*x;
ps1 = @(y) 2*y;
ps2 = @(y) 1-y.^2;

v = ones(n, 1);
B2 = spdiags([v -v], [1, -1], n, n);

% Create diagonal matrices
px1 = ph1(x);
px2 = ph2(x);
py1 = ps1(y);
py2 = ps2(y);

Phi1 = spdiags(px1(:), 0, n, n);
Phi2 = spdiags(px2(:), 0, n, n);
Psi1 = spdiags(py1(:), 0, n, n);
Psi2 = spdiags(py2(:), 0, n, n);

M1 = Phi1 * B2;
M2 = Phi2;
N1 = Psi1;
N2 = B2' * Psi2;

M = {M1, M2};
N = {N1, N2};

% rhs set up
b = @(x, y) sin(pi*x).*cos(pi*y);
rhs1 = b(x, y)';
rhs2 = rhs1;
C = rhs1*rhs2';

% X0 initial guess
X0 = zeros(n, n);