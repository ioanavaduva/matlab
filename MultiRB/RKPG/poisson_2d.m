% n=1000;
h = 1/n; ep = 1;

A = ep*(spdiags([-ones(n, 1) 2*ones(n, 1) -ones(n, 1)],-1:1,n,n))/h^2;

% A_rksm = -(eps*(diag(2*ones(n, 1)) + diag(-1*ones(n-1, 1), 1) + diag(-1*ones(n-1, 1), -1))/h^2);
% I = speye(n);

%% Identity
% rhs1 = I(:,6); rhs2=rhs1;

%% Polynomial
% xtemp = linspace(0,1,n);
% x = repmat(xtemp, 1, n);
% y = reshape(repmat(xtemp, length(xtemp), 1), 1, length(xtemp)^2);
% b = @(x, y) -2.*(6.*x.^2 - 6.*x + 1).*(y-1).^2.*y.^2-2.*(x-1).^2.*x.^2.*(6.*y.^2-6.*y+1);
% rhs = b(x, y)';
% rhs11 = reshape(rhs, n, n); rhs1 = rhs11(:,1);
% rhs2 = rhs1;

%% Sin(pi*x)Sin(pi*y) 
% xtemp = linspace(0,1,n);
% rhs1 = sqrt(2)*pi*sin(pi.*xtemp)';
% rhs2 = rhs1;

%% Sin(pi*x)Cos(pi*y) 
% xtemp = linspace(0,1,n);
% x = repmat(xtemp, 1, n);
% y = reshape(repmat(xtemp, length(xtemp), 1), 1, length(xtemp)^2);
% b = @(x, y) sin(pi.*x).*cos(pi.*y); % without differentiating -- I think satisfy bc
% % b = @(x, y) pi^2*sin(pi.*x).*cos(pi.*y) - pi^2*sin(pi.*x).*cos(pi.*y); % this is just 0 
% rhs = b(x, y)';
% rhs11 = reshape(rhs, n, n); rhs1 = rhs11(:,1);
% rhs2 = rhs1;

%% Ones
rhs1 = ones(n, 1); 
rhs2 = ones(n, 1);

%% Random from the standard normal distribution
% rhs1 = randn(n, 1);
% rhs2 = rhs1;

%% Another random
% NNn = 0.0001;
% rhs1 = NNn*randn(n, 1) + ones(n, 1);
% rhs2 = rhs1;

%% Linspace
% rhs1 = linspace(1, n, n)'; rhs2 = rhs1;

%% Cos(pi*x)
% x = linspace(0, 1, n)';
% rhs1 = cos(pi*x); rhs2 = rhs1;

%% Another random
% rhs1 = sprand(n,1,0.23); rhs2 = rhs1;

%% Alternate 1&-1
% rhs1 = ((-1).^(0:n-1))'; rhs2 = rhs1;