% Setup for MultiRB. 
% 
%       IXN0 + M1XN1 + M2XN2 + M3XN3 = C, 
% 
%   N = [N0, N1, N2, N3] and M = [I, M1, M2, M3].

n = 5; 
h = 1/n;

% A and B
eps = 1;
A = eps*(diag(2*ones(n, 1)) + diag (-1*ones(n-1, 1), 1) + diag (-1*ones(n-1, 1), -1))/h^2;
I = speye(n);

% for N and M
x = linspace(0,1,n);
y = fliplr(x);
ph1 = @(x) zeros(size(x));
ph2 = @(x) zeros(size(x));
ps1 = @(y) zeros(size(x));
ps2 = @(y) zeros(size(x));
% ph1 = @(x) 1-x.^2;
% ph2 = @(x) -2*x;
% ps1 = @(y) 2*y;
% ps2 = @(y) 1-y.^2;

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

M = {M0, N0, M1, M2};
N = {N0, M0, N1, N2};

% rhs set up
rhs1 = cos(2*n*pi*x)';
rhs2 = sin(2*n*pi*y)';
rhs = rhs1*rhs2';
rhs_vec = rhs(:);


% preconditioner
P1 = speye(n); 
P = P1*P1';

% kronecker product form of matrix eq
AA  = kron(A', I) + kron(I', A) + kron(N1', M1) + kron(N2', M2);

%AB = kron(A, I) + kron(I, B);