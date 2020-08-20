
% Driver for 1- sided RKPG
% 
clear all;
addpath(genpath('../../rktoolbox'));

% Setup
n = 500; % size of matrix A
h = 1/n; eps = 1;
A = eps*(diag(2*ones(n, 1)) + diag (-1*ones(n-1, 1), 1) + diag (-1*ones(n-1, 1), -1))/h^2;

% B = A+A';
% 
rhs1 = ones(n, 1);
rhs2 = ones(n, 1);

% rhs1 = randn(n, 1);
% rhs2 = rhs1;

% NNn = 0.0001;
% rhs1 = NNn*randn(n, 1) + ones(n, 1);
% rhs2 = rhs1;

% rhs1 = linspace(1, n, n)'; rhs2 = rhs1;

% x = linspace(0, 1, n)';
% rhs1 = cos(pi*x); rhs2 = rhs1;

% rhs1 = sprand(n,1,0.23); rhs2 = rhs1;
% 
% rhs1 = ((-1).^(0:n-1))'; rhs2 = rhs1;

% %%%% Exact solution --- only need to check bounds for the 3 bases
% AA = kron(A, speye(n))+kron(speye(n), A);
% rhs = rhs1*rhs2';
% rhss = rhs(:);
% Xex = AA\rhss;
% Xex_mat = reshape(Xex, n, n);
% %%%%

tol = 1e-9;
maxit = 300;

% Get smallest and largest eigenvalues
emin = 1e-6; 
opts.tol=1e-4;
% emin = eigs(B, 1, 'SM', opts);
emax = eigs(A, 1,'LA',opts);

% 8 random poles in the spectral interval
poles_rand = [70194.8071105534,61045.6219154057,122475.200195771,127224.208005786,29897.7893219326,78357.5139720186,71289.4347736214,103403.761391083]';

% 8 linspace poles in the spectral interval
poles_linspace = linspace(emin, emax, 8)';

% 6 logspace poles
poles_log = logspace(log10(emin), log10(emax), 6)';

% 4 positive imaginary parts of Zolotarev poles
bb = emax - emin + 1;

k = 4;      % number of poles is 2*k
b = bb;     % sign function on [-10,-1]\cup [1,10]

roots_denom = get_rootsden(k, b);

% work with 6 poles only
sm6 = roots_denom(3:8);
% time & solve using RKPG
tic;
[X1, X2, final_err, vec_res, it, inner_it, avg_inner] = oneside_RKPG(A, rhs1, rhs2, roots_denom, tol,  maxit);
%!!for beckermann bound need to add extra 'upper_bound' to outputs
time = toc;

fprintf('\n Total execution time: %9.4e seconds \n', time)

fprintf('final_err   avg_inner  \n')
fprintf('\n  %9.4e       %d    \n \n', [final_err, avg_inner])

% plot residual v iterations
iter = linspace(1, it, it);
plot(iter, vec_res, 'o');hold on
xlabel('Iterations');
ylabel('Residual');

% plot Beckermann bound on top of residuals
% semilogy(upper_bound, 'x'); hold off;