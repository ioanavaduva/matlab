% To generate the IRKA poles use:
clear all;

n = 100;
h = 1/n; ep = 0.0333;

A = ep*(spdiags([-ones(n, 1) 2*ones(n, 1) -ones(n, 1)],-1:1,n,n))/(h^2);
B = spdiags([-ones(n, 1) zeros(n, 1) ones(n, 1)],-1:1,n,n)/(2*h);

M = A + B;
N = A + B';

rhs1 = ones(n, 1);
rhs2 = rhs1;

% nonuniform_mesh3;

% n=1001;
% separable_coeff;

opts.tol = 1e-4;
emin = eigs(M, 1, 'smallestabs', opts);
emax = eigs(M, 1,'largestabs',opts);

m = 16; % number of poles, can change
xi = emin + (emax-emin)*rand(1,m); % initial choice of poles

tic;
[shifts1, shifts2, its] = irka_shifts2(M, N, rhs1, rhs2, xi,  xi, 1e-2); 
toc;
shifts1 = sort(shifts1, 'desc');
shifts2 = sort(shifts2, 'desc');

[shifts, it] = irka_shifts(M, rhs1, xi, 1e-2);
