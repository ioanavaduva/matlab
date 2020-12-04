% setup block rhs with irka
clear all;
addpath(genpath('../../rktoolbox'));

n = 100; % size of matrix A 
h = 1/n; eps = 1;
A = eps*(diag(2*ones(n, 1)) + diag (-1*ones(n-1, 1), 1) + diag (-1*ones(n-1, 1), -1))/h^2;

% rhs1 = [ones(n, 1), A*ones(n, 1)];
rhs1 = ones(n, 2);
B = rhs1*rhs1';

opts.tol=1e-4;
emin = eigs(A, 1,'smallestabs',opts);
emax = eigs(A, 1,'largestabs',opts);

tol = 1e-4;

% compute roots denom from Zolotarev problem
bb = emax - emin + 1;
k = 4;      % number of poles is 2*k
poles = get_rootsden(k, bb);

m = 8; % number of poles; can change
xi = emin + (emax - emin)*rand(1,m);
tic;
B_hat = tangential_dir(A, B, xi);
[shifts, it] = irka_shifts_block(A, B, B_hat, xi, tol);
toc;
