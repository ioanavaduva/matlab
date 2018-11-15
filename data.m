%function b
xtemp = [0.1:0.1:1];
x = repmat(xtemp, 1, 10);
y = reshape(repmat(xtemp, length(xtemp), 1), 1, length(xtemp)^2);
b = @(x, y) sin(pi*x).*cos(pi*y);
b_use = b(x, y)';
%initial guess x0
x0 = zeros(length(b_use), 1);
%tolerance
TOL = 10^-8;
%alowed number of iterations
maxit = 500;
%Poisson matrix
A = poissonmatrix(10); %actual matrix A
D = diag(diag(A)); %diagonal matrix from A
L = tril(A, -1); %strictly lower triangular matrix from A
U = triu(A, 1); %strictly upper triangular matrix from A
%test matrix
C = [4 1 1 0; 1 4 1 1; 1 1 4 1; 0 1 1 4]; %test matrix C
c = [-2; 4; 10; 14] ; %vector of results 
z = zeros(4, 1); %initial guess
D1 = diag(diag(C)); %diagonal matrix from C
L1 = tril(C, -1); %strictly lower triangular matrix from C
U1 = triu(C, 1); %strictly upper triangular matrix from C
% not positive definite matrix
N = [1, 2; 2, 1];

