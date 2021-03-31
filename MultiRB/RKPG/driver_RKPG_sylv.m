clear all;
addpath(genpath('../../rktoolbox'));

% Setup
n = 100; % size of matrix A
h = 1/n; eps = 1;
T = eps*(diag(2*ones(n, 1)) + diag (-1*ones(n-1, 1), 1) + diag (-1*ones(n-1, 1), -1))/h^2;
I = eye(n);
A = kron(T, I) + kron(I, T);

rhsA = ones(n^2, 1);
rhsT = ones(n, 1);

tol = 1e-9;
maxit = 300;

opts.tol = 1e-4;
emin = eigs(A, 1, 'smallestabs', opts);
emax = eigs(A, 1,'largestabs',opts);

bb = emax - emin + 1;
k = 4;      % number of poles is 2*k
[roots_denom, extrema] = get_rootsden(k, bb);
roots_denom = roots_denom + emin - 1;

tic;
[X1, X2, final_err, vec_res, it, inner_it, avg_inner] = RKPG_sylv(A, T, rhsA, rhsT, roots_denom, tol,  maxit);
time = toc;
fprintf('\n Total execution time: %9.4e seconds \n', time)

fprintf('final_err   avg_inner  \n')
fprintf('\n  %9.4e       %d    \n \n', [final_err, avg_inner])