% data for Poisson 2D problem with n degrees of freedom
n = 1023;
h = 1/n;

xtemp = linspace(0,1,n);
x = repmat(xtemp, 1, n);
y = reshape(repmat(xtemp, length(xtemp), 1), 1, n^2);
b = @(x, y) sin(pi*x).*cos(pi*y);
rhs = b(x, y)';

% initial guess x0
x0 = zeros(length(rhs), 1);

T = (diag(2*ones(n, 1)) + diag (-1*ones(n-1, 1), 1) + diag (-1*ones(n-1, 1), -1))/h^2;
B = speye(n);

% Poisson matrix
A = kr_pois(n);

