% To generate the IRKA poles use:
clear all;

n = 1000;
% poisson_2d;

nonuniform_mesh3;

% n=1001;
% separable_coeff;

opts.tol = 1e-4;
emin = eigs(A, 1, 'smallestabs', opts);
emax = eigs(A, 1,'largestabs',opts);

m = 16; % number of poles, can change
xi = emin + (emax-emin)*rand(1,m); % initial choice of poles

tic;
[shifts, its] = irka_shifts(A, rhs1, xi, 1e-2); % original irka poles
toc;
shifts = sort(shifts, 'desc');

% tic;
% [shifts_p, its_p] = irka_shifts_p(A, rhs1, xi, 1e-2); % eigenvalue decomposition irka poles
% toc;
% shifts_p = sort(shifts_p, 'desc');