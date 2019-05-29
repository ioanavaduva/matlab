nh=50;
T=diag(2*ones(nh,1))+diag(-ones(nh-1,1),1)+diag(-ones(nh-1,1),-1);
I=speye(nh);
A=-T;

n=nh^2;

B=A;

[LA,UA]=lu(A);
[LB,UB]=lu(B);

x = linspace(0,1,nh);
y = linspace(1, 0, nh);
b = @(x, y) sin(pi.*x).*cos(pi.*y);
rhs1 = b(x, y)';
rhs2=rhs1;

m=100;
tol=1e-9;
[X1,X2,r]=kpik_sylv(A,LA,UA,B,LB,UB,rhs1,rhs2,m,tol);

fprintf('final true absolute residual norm: \n')
disp(norm(A*X1*X2'+X1*X2'*B'+rhs1*rhs2'))    %this matrix should never be formed for n large 