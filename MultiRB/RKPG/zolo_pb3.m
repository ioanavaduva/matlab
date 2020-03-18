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

k = 4;      % rational degree
b = bb;     % sign function on [-10,-1]\cup [1,10]
r = rkfun.gallery('sign', k, b);

% Compute the minimum of Z Problem 4
K = ellipke(1-1/b^2);
[sn, cn, dn] = ellipj((0:k)*K/k, 1-1/b^2);
extrema = b*dn;
vals = 1-r(extrema);
c = mean( vals(1:2:end) );
e = eig( [ 2-4/c^2 1 ; -1 0 ] );
Zk = min(abs(e));

% Obtain the polynomials p and q of the rational function r = p/q 
[p,q,pq] = poly(r);
pp = @(z) polyval(p, z);
qq = @(z) polyval(q, z);

% Transform Z Problem 4 into Problem 3 using eq. (2) - rearranged in Theorem 2.1
% (Istace/Thiran paper)and have the denominator given by
denom = @(z) qq(z).*(1-Zk) - pp(z).*(1+Zk);

%%%% Questions: What should I make z in order to find the roots of denom?
%%%% Should it be the same as k, i.e. the degree of the rational polynomial?

z = linspace(-b, b, k);
poly(denom(z)) % these values are very large. what am I doing wrong??

