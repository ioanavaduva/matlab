%% Driver for RKPG (with 2-sided projection)

clear all;
addpath(genpath('../../rktoolbox'));

%% %% Setup
n = 1000; % size of matrix A 
h = 1/n; eps = 1;
A = eps*(diag(2*ones(n, 1)) + diag(-1*ones(n-1, 1), 1) + diag(-1*ones(n-1, 1), -1))/h^2;
% A_rksm = -(eps*(diag(2*ones(n, 1)) + diag(-1*ones(n-1, 1), 1) + diag(-1*ones(n-1, 1), -1))/h^2);
I = speye(n);

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

%% %% Choose RHS
%% Identity
% rhs1 = I(:,6); rhs2=rhs1;

%% Polynomial
% xtemp = linspace(0,1,n);
% x = repmat(xtemp, 1, n);
% y = reshape(repmat(xtemp, length(xtemp), 1), 1, length(xtemp)^2);
% b = @(x, y) -2.*(6.*x.^2 - 6.*x + 1).*(y-1).^2.*y.^2-2.*(x-1).^2.*x.^2.*(6.*y.^2-6.*y+1);
% rhs = b(x, y)';
% rhs11 = reshape(rhs, n, n); rhs1 = rhs11(:,1);
% rhs2 = rhs1;

%% Sin(pi*x)Sin(pi*y) 
% xtemp = linspace(0,1,n);
% % x = repmat(xtemp, 1, n);
% % y = reshape(repmat(xtemp, length(xtemp), 1), 1, length(xtemp)^2);
% % b = @(x, y) sin(pi.*x).*sin(pi.*y); % without differentiating -- gives NaNs
% % b = @(x, y) pi^2*sin(pi.*x).*sin(pi.*y) + pi^2*sin(pi.*x).*sin(pi.*y);
% rhs1 = sqrt(2)*pi*sin(pi.*xtemp)';
% % rhs = b(x, y)';
% % rhs11 = reshape(rhs, n, n); rhs1 = rhs11(:,1);
% rhs2 = rhs1;

%% Sin(pi*x)Cos(pi*y) 
% xtemp = linspace(0,1,n);
% x = repmat(xtemp, 1, n);
% y = reshape(repmat(xtemp, length(xtemp), 1), 1, length(xtemp)^2);
% b = @(x, y) sin(pi.*x).*cos(pi.*y); % without differentiating -- I think satisfy bc
% % b = @(x, y) pi^2*sin(pi.*x).*cos(pi.*y) - pi^2*sin(pi.*x).*cos(pi.*y); % this is just 0 
% rhs = b(x, y)';
% rhs11 = reshape(rhs, n, n); rhs1 = rhs11(:,1);
% rhs2 = rhs1;

%% Ones
rhs1 = ones(n, 1); 
rhs2 = ones(n, 1);

%% Random from the standard normal distribution
% rhs1 = randn(n, 1);
% rhs2 = rhs1;

%% Another random
% NNn = 0.0001;
% rhs1 = NNn*randn(n, 1) + ones(n, 1);
% rhs2 = rhs1;

%% Linspace
% rhs1 = linspace(1, n, n)'; rhs2 = rhs1;

%% Cos(pi*x)
% x = linspace(0, 1, n)';
% rhs1 = cos(pi*x); rhs2 = rhs1;

%% Another random
% rhs1 = sprand(n,1,0.23); rhs2 = rhs1;

%% Alternate 1&-1
% rhs1 = ((-1).^(0:n-1))'; rhs2 = rhs1;

%% %% Exact solution --- only need to check bounds for the 3 bases & to
%% compare with 1 sided projection
% AA = kron(A, speye(n))+kron(speye(n), A);
% rhs = rhs1*rhs2';
% rhss = rhs(:);
% Xex = AA\rhss;
% Xex_mat = reshape(Xex, n, n);
%%%%

tol = 1e-9;
maxit = 300;

%% %% Get smallest and largest eigenvalues
opts.tol=1e-4;
emin = eigs(A, 1, 'smallestabs', opts);
emax = eigs(A, 1,'largestabs',opts);
% emax_rksm = eigs(A_rksm, 1,'smallestabs',opts);
% emin_rksm = eigs(A_rksm, 1, 'largestabs', opts);

%% %% Different poles
%% 8 random poles in the spectral interval
% poles_rand = emin + rand(1,8)*(emax - emin)';

%% 8 linspace poles in the spectral interval
% poles_linspace = linspace(emin, emax, 8)';

%% 6 logspace poles
% poles_log = logspace(log10(emin), log10(emax), 8)';

%% Get nodes (Sabino thesis)
% s_nodes = 8;                           % Choose 2 nodes (could vary)
% tic;
% snew = get_nodes2(emin,emax,s_nodes);      % Use interval for A_1;
% toc;
% s_parameter=sort(snew, 'desc');

%% Zolotarev
bb = emax - emin + 1;
k = 4;      % number of poles is 2*k
tic;
[roots_denom, extrema] = get_rootsden(k, bb);
toc;
roots_denom = roots_denom + emin - 1;
roots = sort(extrema + emin - 1);
%% IRKA
% n_small = 100;
% h_small = 1/n_small;
% A_small = eps*(diag(2*ones(n_small, 1)) + diag(-1*ones(n_small-1, 1), 1) + diag(-1*ones(n_small-1, 1), -1))/h_small^2;
% A_small_scaled = eps*(diag(2*ones(n_small, 1)) + diag(-1*ones(n_small-1, 1), 1) + diag(-1*ones(n_small-1, 1), -1))/h^2;
% rhs1_small = ones(n_small, 1);
% % % % rhs1_small = ((-1).^(0:n_small-1))';
% % % e_sm = eigs(A_small_scaled, 1, 'smallestabs', opts);
% % % e_lg =  eigs(A_small_scaled, 1, 'largestabs', opts);
% % % lam_min = eigs(A_small, 1, 'smallestabs', opts);
% % % lam_max =  eigs(A_small, 1, 'largestabs', opts);
% m = 10; % number of poles; can change
% % % % xi = e_sm + (e_lg-e_sm)*rand(1,m);
% xi = emin + (emax-emin)*rand(1,m);
% % % % tic;
% [shifts3,its3] = irka_shifts(A_small,rhs1_small, xi, 1e-4);
% % % shifts3 = [4.000e+5; shifts3];
% [shifts2,its2] = irka_shifts(A_small_scaled,rhs1_small, xi, 1e-4);
% % % % toc;
% % % % tic;
% % % % [shifts, its] = irka_shifts(A, rhs1, xi, 1e-2);
% % % % toc;
% % tic;
% [shifts, its] = irka_shifts(A, rhs1, xi, 1e-4);
% % % toc;
% % % % tic;
% % % % [shifts3, its3] = irka_shifts(A, rhs1, xi, 1e-1);
% % % % toc;
irka_small_prob_poles;
%% Adaptive (Druskin & Simoncini)
% [Z, nrmrestot, pol] = rksm(A_rksm, I, I, rhs1, maxit, tol, emin_rksm, emax_rksm, 0, 1e-12);
% [pols,m1,n1] = uniquetol(pol);
% [c1,d1] = sort(m1);
% pols_dist = pols(d1);
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
[X1, X2, final_err, vec_res, it, inner_it, avg_inner, e_Ap] = RKPG(A, rhs1, rhs2, shifts, tol,  maxit);
%!! for beckermann bound need to add extra 'upper_bound' to outputs
%!! to check errors need error_vec in outputs and Xex_mat (exact solution in matrix form) in inputs
%!! plot Ritz values need e_Ap in outputs
time = toc;

fprintf('\n Total execution time: %9.4e seconds \n', time)

fprintf('final_err   avg_inner  \n')
fprintf('\n  %9.4e       %d    \n \n', [final_err, avg_inner])

%% plot residual vs. iterations
iter = linspace(1, it, it);
semilogy(iter, vec_res, 'p'); hold on
xlabel('Iterations');
ylabel('Residual');

%% plot error on top of residuals
% plot(error_vec, 'o');hold on;

%% plot Beckermann bound on top of residuals
% semilogy(upper_vec, 'x'); hold off;