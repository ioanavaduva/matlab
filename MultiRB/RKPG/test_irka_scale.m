n=100;
A = full(gallery('tridiag',n,-1,2,-1));
A2 = A/67;
% rhs1 = ones(n, 1);
xtemp = linspace(0,1,n);
x = repmat(xtemp, 1, n);
y = reshape(repmat(xtemp, length(xtemp), 1), 1, length(xtemp)^2);
b = @(x, y) -2.*(6.*x.^2 - 6.*x + 1).*(y-1).^2.*y.^2-2.*(x-1).^2.*x.^2.*(6.*y.^2-6.*y+1);
rhs = b(x, y)';
rhs11 = reshape(rhs, n, n); rhs1 = rhs11(:,1);
rhs2 = rhs1;

opts.tol=1e-4;
emin = eigs(A, 1, 'smallestabs', opts);
emax = eigs(A, 1,'largestabs',opts);

m = 10;
xi = emin + (emax-emin)*rand(1,m);
[shifts, its] = irka_shifts(A, rhs1, xi, 1e-4);
shifts2 = irka_shifts(A2, rhs1, xi, 1e-4);