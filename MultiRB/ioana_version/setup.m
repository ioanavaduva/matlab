% Setup for MultiRB. 
% We want to solve a problem of the form 
%       AXI + M0XN0 + M1XN1 + M2XN2 = C, 
% with N = [I, N0, N1, N2] and M = [A, M0, M1, M2].

n = 64; 
h = 1/n;

% A and B
eps = 1;
A = eps*(diag(2*ones(n, 1)) + diag (-1*ones(n-1, 1), 1) + diag (-1*ones(n-1, 1), -1))/h^2;
I = speye(n);

% for N and M
x = linspace(0,1,n);
y = x;
ph1 = @(x) 1-x.^2;
ph2 = @(x) -2*x;
ps1 = @(y) 2*y;
ps2 = @(y) 1-y.^2;

v = ones(n, 1);
B2 = spdiags([v -v], [1, -1], n, n)/(2*h);

% Create diagonal matrices
px1 = ph1(x);
px2 = ph2(x);
py1 = ps1(y);
py2 = ps2(y);

Phi1 = spdiags(px1(:), 0, n, n);
Phi2 = spdiags(px2(:), 0, n, n);
Psi1 = spdiags(py1(:), 0, n, n);
Psi2 = spdiags(py2(:), 0, n, n);

M0 = I; N0 = A; M1 = Phi1 * B2; M2 = Phi2; N1 = Psi1; N2 = B2' * Psi2; 

M = {N0, M0, M1, M2};
N = {M0, N0, N1, N2};

% rhs set up
xtemp = linspace(0,1,n);
x = repmat(xtemp, 1, n);
y = reshape(repmat(xtemp, length(xtemp), 1), 1, length(xtemp)^2);
b = @(x, y) sin(pi*x).*cos(pi*y);
rhs = b(x, y)';
rhs1 = reshape(rhs, n, n);
rhs2 = I;

% preconditioner
P1 = speye(n); 
P = P1*P1';

% kronecker product form of matrix eq
AA  = kron(A', I) + kron(I', A) + kron(N1', M1) + kron(N2', M2);

%AB = kron(A, I) + kron(I, B);
