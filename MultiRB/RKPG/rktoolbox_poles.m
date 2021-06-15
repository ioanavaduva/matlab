n = 1000;
h = 1/n; eps = 1;
A = eps*(diag(2*ones(n, 1)) + diag(-1*ones(n-1, 1), 1) + diag(-1*ones(n-1, 1), -1))/h^2;

opts.tol=1e-4;
emin = eigs(A, 1, 'smallestabs', opts);
emax = eigs(A, 1,'largestabs',opts);

bb = emax - emin + 1;
k = 4;      % number of poles is 2*k
[roots_denom, extrema] = get_rootsden(k, bb);
roots_denom = roots_denom + emin - 1;
roots_denom = sort(roots_denom, 'desc');