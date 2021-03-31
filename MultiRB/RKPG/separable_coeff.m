% n = 1001;
g = 0:1/n:1; 
gs = 1/(2*n):1/n:1;

% a = @(x) 1+x; 
a = @(x) sin(x);
% a = @(x) 1 +  50*exp(-5*x.^2);
% syms x
% a = @(x) ((0<x) & (x<0.5)).*1 + ((0.5<x) & (x<1)).*1e+4;

rhs = ones(n-1, 1);

% rhs = ((-1).^(0:n-2))'; 

% xtemp = linspace(0,1,n-1);
% x = repmat(xtemp, 1, n-1);
% y = reshape(repmat(xtemp, length(xtemp), 1), 1, length(xtemp)^2);
% b = @(x, y) -2.*(6.*x.^2 - 6.*x + 1).*(y-1).^2.*y.^2-2.*(x-1).^2.*x.^2.*(6.*y.^2-6.*y+1);
% rs = b(x, y)';
% rhs11 = reshape(rs, n-1, n-1); rhs = rhs11(:,1);

% xtemp = linspace(0, 1, n-1);
% rhs = sqrt(2)*pi*sin(pi.*xtemp)';

% x = linspace(0, 1, n-1)';
% rhs = cos(pi*x);

ag = a(g(2:n));
ags = a(gs);

D = diag(ag);
D_inv = diag(sqrt(1./ag));

h = 1/n;
T_orig = gallery('tridiag', -ags(2:n-1), ags(1:n-1)+ ags(2:n), -ags(2:n-1))/h^2;

A = D_inv*T_orig*D_inv; % Scaled version of T 

rhs1 = D_inv*rhs; % Scaled version of rhs
rhs2 = rhs1;
% Fd = f*f';



