% Driver for RKPG

clear all;
addpath(genpath('../rktoolbox'));

% Setup
n = 15; % size of matrix A
h = 1/n; eps = 1;
A = eps*(diag(2*ones(n, 1)) + diag (-1*ones(n-1, 1), 1) + diag (-1*ones(n-1, 1), -1))/h^2;
rhs1 = ones(n,1);
rhs2 = rhs1;

tol = 1e-10;
maxit = 20;

% Get smallest and largest eigenvalues
emin = 1e-6; 
opts.tol=1e-4;
emax = eigs(A, 1,'LA',opts);

% 6 logspace poles
poles_log = logspace(log10(emin), log10(emax), 6)';

% 4 positive imaginary parts of Zolotarev poles
bb = emax - emin + 1;

k = 4;      % rational degree
b = bb;     % sign function on [-10,-1]\cup [1,10]
r = rkfun.gallery('sign', k, b);
po = imag(poles(r));
poles_Zolo = po(po >= 0);

% time & solve using RKPG
tic;
[X1, X2, final_err, vec_res, it, inner_it, avg_inner] = RKPG(A, rhs1, rhs2, poles_Zolo, tol,  maxit);
time = toc;

fprintf('\n Total execution time: %9.4e seconds \n', time)

fprintf('final_err   avg_inner  \n')
fprintf('\n  %9.4e       %4.2f    \n \n', [final_err, avg_inner])

