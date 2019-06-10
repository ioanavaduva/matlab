nh = 50;
h = 1/nh;
T = diag(2*ones(nh,1))+diag(-ones(nh-1,1),1)+diag(-ones(nh-1,1),-1);
I = speye(nh);
A = -T/(h^2);

B=A;

[LA,UA] = lu(A);
[LB,UB] = lu(B);

x = linspace(0,1,nh);
y = linspace(1, 0, nh);
b = @(x, y) sin(pi.*x).*cos(pi.*y);
rhs1 = b(x, y)';
rhs2 = rhs1;
%rhs1 = ones(nh, 1);
%rhs2 = rhs1;

m = 300;
tol = 1e-9;
[X1,X2,r] = kpik_sylv(A,LA,UA,B,LB,UB,rhs1,rhs2,m,tol);

fprintf('final true absolute residual norm: \n')
disp(norm(A*X1*X2'+X1*X2'*B'+rhs1*rhs2'))    %this matrix should never be formed for n large 

% test if we get the same as when using backslash Kronecker form - only do
% for small problems
X = X1*X2';
X = X(:);
AA = -(kron(A,I)+kron(I,A));
rhs = rhs1*rhs2';
rhs = rhs(:);
U = AA\rhs;
norm(U-X)