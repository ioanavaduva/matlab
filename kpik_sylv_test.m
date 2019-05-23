nh=100;
T=diag(2*ones(nh,1))+diag(-ones(nh-1,1),1)+diag(-ones(nh-1,1),-1);
I=speye(nh);
A=-(kron(T,I)+kron(I,T));

n=nh^2;

B=A;

[LA,UA]=lu(A);
[LB,UB]=lu(B);

xtemp = linspace(0,1,nh);
x = repmat(xtemp, 1, nh);
y = reshape(repmat(xtemp, length(xtemp), 1), 1, length(xtemp)^2);
b = @(x, y) sin(pi.*x).*cos(pi.*y);
rhs1 = b(x, y)';

rhs2=I(:);

m=100;
tol=1e-9;
[X1,X2,r]=kpik_sylv(A,LA,UA,B,LB,UB,rhs1,rhs2,m,tol);

fprintf('final true absolute residual norm: \n')
disp(norm(A*X1*X2'+X1*X2'*B'+rhs1*rhs2'))    %this matrix should never be formed for n large 