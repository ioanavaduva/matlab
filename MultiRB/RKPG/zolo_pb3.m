clear all;
addpath(genpath('../../rktoolbox'));

% Setup
n = 1000; % size of matrix A
h = 1/n; eps = 1;
A = eps*(diag(2*ones(n, 1)) + diag (-1*ones(n-1, 1), 1) + diag (-1*ones(n-1, 1), -1))/h^2;

emin = 1e-6; 
opts.tol=1e-4;
emax = eigs(A, 1,'LA',opts);

% Compute the solution to Z Problem 4 of rational degree k
bb = emax - emin + 1;

k = 2;      % rational degree
b = bb;     % sign function on [-10,-1]\cup [1,10]
r = rkfun.gallery('sign', k, b);  % Zolotarev's sign approximation on +/-[1, b] 
    % of type (2*d-1, 2*d) -- this means the degree of p is 2d-1 and the degree of q is 2d.

% Compute the minimum of Z Problem 3 (by transforming c, the min of Prob 4)
K = ellipke(1-1/b^2);
[sn, cn, dn] = ellipj((0:k)*K/k, 1-1/b^2);
extrema = b*dn;
vals = 1-r(extrema);
c = mean( vals(1:2:end) );
e = eig( [ 2-4/c^2 1 ; -1 0 ] ); % changed to have a -1 in one of the off-diagonals
Zk = min(abs(e));

% Obtain the polynomials p and q of the rational function r = p/q 
[p,q,pq] = poly(r); % p has degree 1 smaller than q, so add a 0 
% poly gives the coefficients of 3x+1 as vector [3, 1]
pp = [0, p]; % now p and q have the same degree

% Transform Z Problem 4 into Problem 3 using eq. (2) - rearranged in Theorem 2.1
% (Istace/Thiran paper)and have the denominator given by
denom = q.*(1-Zk) - pp.*(1+Zk);
roots_denom = roots(denom);



