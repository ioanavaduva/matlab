% Set up
n = 100; % can change
A = full(gallery('tridiag',n,-1,2,-1));
f = ones(n,1); % can change
In = eye(n);

% Compute spectral interval
lam_max = max(eig(A));
lam_min = min(eig(A));

%%% Random poles from spectral interval
m = 4; % number of poles; can change
xi = lam_min + (lam_max-lam_min)*rand(1,m);
%%%

%%% Ioana's roots denominator
bb = lam_max - lam_min + 1; 
b = bb; 
k = 2; % number of poles from this is 2*k 
r = rkfun.gallery('sign', k, b); % Zolotarev 'sign' approximation (Prob 4) 
% of degree 2*k on the union of [1,b] and [-b,-1]
po = imag(poles(r));
poles_Zolo = po(po >= 0);
% Compute the minimum of Z in Problem 3 (by transforming c, the min of Prob 4)
K = ellipke(1-1/b^2);
[sn, cn, dn] = ellipj((0:k)*K/k, 1-1/b^2);
extrema = b*dn;
vals = 1-r(extrema);
c = mean( vals(1:2:end) );
e = eig( [ 2-4/c^2 1 ; -1 0 ] );
Zk = min(abs(e));
% Obtain the polynomials p and q of the rational function r = p/q 
[p,q,pq] = poly(r);
pp = [0, p];
% Transform Z Problem 4 into Problem 3 
denom = q.*(1-Zk) - pp.*(1+Zk);
roots_denom = roots(denom)';
%%%

[shifts, it] = irka_shifts(A, f, roots_denom, 1e-4);
[shifts2, its] = irka_shifts(A, f, xi, 1e-4);

