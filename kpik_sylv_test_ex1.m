n = 1000;

h = 1/n; ep = 0.0083;

A = ep*(spdiags([-ones(n, 1) 2*ones(n, 1) -ones(n, 1)],-1:1,n,n))/(h^2);
B = spdiags([-ones(n, 1) zeros(n, 1) ones(n, 1)],-1:1,n,n)/(2*h);

% xtemp = linspace(0, 1, n);
% w1 = @(x) 1 + ((x+1).^2)/4;
% phi = diag(w1(xtemp));

w1 = @(x) 1 - (2*x + 1).^2;
phi = diag(w1(xtemp));

w2 = @(x) 1 - x.^2;
psi = diag(w2(xtemp));

M = A + phi*B; 
N = A + B*psi;

[LM,UM] = lu(M);
[LN,UN] = lu(N);

rhs1 = ones(n, 1);
rhs2 = -rhs1;

tic;[X1,X2,r] = kpik_sylv(M,LM,UM,N,LN,UN,rhs1,rhs2,300,1e-8);toc;