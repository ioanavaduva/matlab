%% Driver for RKPG (with 2-sided projection)

clear all;
addpath(genpath('../../rktoolbox'));

%% 2D Poisson
% n = 1000; 
% poisson_2d;

% B = spdiags([-ones(n, 1) zeros(n, 1) ones(n, 1)],-1:1,n,n)/(2*h);
% % B = zeros(n, n);
% 
% M = A + B;
% N = A + B';
%% small poisson prob
% ns=100;
% hs = 1/ns; eps = 1;
% A_small = eps*(diag(2*ones(ns, 1)) + diag(-1*ones(ns-1, 1), 1) + diag(-1*ones(ns-1, 1), -1))/hs^2;
% 
% rhs1_small = ones(ns, 1);
%% 2D Variable Diffusion Coefficients
n = 1001;
separable_coeff;

%% Non-uniform Mesh
% n = 1000;
% nonuniform_mesh3;

%% Sparse, random, SPD matrix
% density = 2/n; % for example
% rc = 0.1; % Reciprocal condition number
% A = sprandsym(m, density, rc, 1);

%% Indefinite matrix B
% A = rand(n,n);
% B = A+A';
% if all(eig(B)>0)
%     return
% end

%% %% Exact solution --- only need to check bounds for the 3 bases & to
%% compare with 1 sided projection
% AA = kron(A, speye(n))+kron(speye(n), A);
% rhs = rhs1*rhs2';
% rhss = rhs(:);
% Xex = AA\rhss;
% Xex_mat = reshape(Xex, n, n);
%%%%

% [Q,D]=eig(A); d = diag(D);

% A = sparse(A);
tol = 1e-9;
maxit = 300;

%% %% Get smallest and largest eigenvalues
opts.tol=1e-4;
emin = eigs(A, 1, 'smallestabs', opts);
emax = eigs(A, 1,'largestabs',opts);


% emax_rksm = eigs(A_rksm, 1,'smallestabs',opts);
% emin_rksm = eigs(A_rksm, 1, 'largestabs', opts);

% emin_small = eigs(A_small, 1, 'smallestabs', opts);
% emax_small = eigs(A_small, 1,'largestabs',opts);

%% %% Different poles
%% 8 random poles in the spectral interval
% poles_rand = emin + rand(1,16)*(emax - emin)';

%% 8 linspace poles in the spectral interval
% poles_linspace = linspace(emin, emax, 8)';

%% logspace poles
tic;
poles_log = logspace(log10(emin), log10(emax), 16)';
toc;
poles_log = sort(poles_log, 'desc');

%% Get nodes (Sabino thesis)
s_nodes = 16;                           % Choose 2 nodes (could vary)
tic;
snew = get_nodes2(emin,emax,s_nodes);      % Use interval for A_1;
toc;
s_parameter=sort(snew, 'desc');

tic; % for geometric mesh only
param = sabino_approx(emin, emax, 16);
toc;
param = sort(param, 'desc');
%% Zolotarev
% bb = emax - emin + 1;
% k = 4;      % number of poles is 2*k
% tic;
% [roots_denom, ~] = get_rootsden(k, bb);
% toc;
% roots_denom = roots_denom + emin - 1;
% roots_denom = sort(roots_denom, 'desc');
% % roots = sort(extrema + emin - 1);

%% IRKA
m = 16; % number of poles; can change
xi = emin + (emax-emin)*rand(1,m);
tic;
[shifts, its] = irka_shifts(A, rhs1, xi, 1e-2);
toc;
shifts = sort(shifts, 'desc');

%% Adaptive (Druskin & Simoncini)
% tic;
% [Z, nrmrestot, pol] = rksm(A_rksm, I, I, rhs1, maxit, tol, emin_rksm, emax_rksm, 0, 1e-12);
% toc;
% [pols,m1,n1] = uniquetol(pol);
% [c1,d1] = sort(m1);
% pols_dist = sort(pols(d1), 'desc');
% pols = sort(pols, 'desc');
% fprintf('\n Number of rksm unique poles: %d \n', length(pols_dist))
% fprintf('\n Number of rksm poles (total): %d \n', length(pol))

%% %% Different combinations of poles
%% work with 6 poles only
% sm6 = roots_denom(3:8); % 6 smallest poles
% la6 = roots_denom(1:6); % 6 largest poles
% s3l3 = zeros(6,1); 
% s3l3(1:3) = roots_denom(1:3);
% s3l3(4:6) = roots_denom(6:8);
% rev_s3l3 = flip(s3l3);

%% work with 4 poles only
% sm4 = roots_denom(5:8); % 4 smallest poles
% la4 = roots_denom(1:4); % 4 largest poles
% l2s2 = zeros(4,1); 
% l2s2(1:2) = roots_denom(1:2);
% l2s2(3:4) = roots_denom(7:8);
% s2l2 = flip(l2s2);

%% reversed order poles
% rev_poles = flip(roots_denom);
% % random order of poles
% rand_poles = roots_denom(randperm(length(roots_denom)));
% % specified order of poles so that smallest first, largest second, etc
% vec_srt = [8, 1, 7, 2, 6, 3, 5, 4];
% sort_poles = roots_denom(vec_srt);

%% time & solve using RKPG
tic;
[X1, X2, final_err, vec_res, it, inner_it, avg_inner, Ap, rhsp] = RKPG(A, rhs1, rhs2, shifts, 1e-8, 300);

% !! for beckermann bound need to add extra 'upper_bound' to outputs
% !! to check errors need error_vec in outputs and Xex_mat (exact solution in matrix form) in inputs
% !! plot Ritz values need e_Ap in outputs
time = toc;



fprintf('\n Total execution time: %9.4e seconds \n', time)

fprintf('final_err   avg_inner  \n')
fprintf('\n  %9.4e       %d    \n \n', [final_err, avg_inner])

% tic;
% [X1_p, X2_p, final_err_p, vec_res_p, it_p, inner_it_p, avg_inner_p, V_p] = RKPG_p(A, rhs1, rhs2, Q, d, shifts, 1e-8, maxit);
% toc;
%% plot residual vs. iterations
% iter = linspace(1, it, it);
% semilogy(iter, vec_res, 'p'); hold on
% xlabel('Iterations');
% ylabel('Residual');

%% plot error on top of residuals
% plot(error_vec, 'o');hold on;

%% plot Beckermann bound on top of residuals
% semilogy(upper_vec, 'x'); hold off;