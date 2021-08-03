% n=1000;
h = 1/n; ep = 1;

A = ep*(spdiags([-ones(n, 1) 2*ones(n, 1) -ones(n, 1)],-1:1,n,n))/h^2;


%% Polynomial rhs
% xtemp = linspace(0,1,n);
% x = repmat(xtemp, 1, n);
% y = reshape(repmat(xtemp, length(xtemp), 1), 1, length(xtemp)^2);
% b = @(x, y) -2.*(6.*x.^2 - 6.*x + 1).*(y-1).^2.*y.^2-2.*(x-1).^2.*x.^2.*(6.*y.^2-6.*y+1);
% rhs = b(x, y)';
% rhs11 = reshape(rhs, n, n); rhs1 = rhs11(:,1);
% rhs2 = rhs1;

%% Ones
rhs1 = ones(n, 1); 
rhs2 = ones(n, 1);

%% Cos(pi*x)
% x = linspace(0, 1, n)';
% rhs1 = cos(pi*x); rhs2 = rhs1;

%% Alternate 1&-1
% rhs1 = ((-1).^(0:n-1))'; rhs2 = rhs1;