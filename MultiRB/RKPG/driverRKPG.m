
% Driver for RKPG (with 2-sided projection)

clear all;
addpath(genpath('../../rktoolbox'));

% Setup
n = 500; % size of matrix A 
h = 1/n; eps = 1;
A = eps*(diag(2*ones(n, 1)) + diag (-1*ones(n-1, 1), 1) + diag (-1*ones(n-1, 1), -1))/h^2;

%%% get indefinite matrix B
% A = rand(n,n);
% B = A+A';
% if all(eig(B)>0)
%     return
% end
%%%

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
 
% rhs1 = ((-1).^(0:n-1))'; rhs2 = rhs1;

%%%% Exact solution --- only need to check bounds for the 3 bases & to
%%%% compare with 1 sided projection
AA = kron(A, speye(n))+kron(speye(n), A);
rhs = rhs1*rhs2';
rhss = rhs(:);
Xex = AA\rhss;
Xex_mat = reshape(Xex, n, n);
%%%%

tol = 1e-9;
maxit = 300;

% Get smallest and largest eigenvalues
emin = 1e-6; 
opts.tol=1e-4;
% emin = eigs(B, 1, 'SM', opts);
emax = eigs(A, 1,'LA',opts);

% 8 random poles in the spectral interval
poles_rand = emin + rand(1,8)*(emax - emin)';

% 8 linspace poles in the spectral interval
poles_linspace = linspace(emin, emax, 8)';

% 6 logspace poles
poles_log = logspace(log10(emin), log10(emax), 6)';

% 4 positive imaginary parts of Zolotarev poles
bb = emax - emin + 1;

k = 4;      % rational degree
b = bb;     % sign function on [-10,-1]\cup [1,10]
r = rkfun.gallery('sign', k, b);
% poles(r)
po = imag(poles(r));
poles_Zolo = po(po >= 0);

% % Compute the minimum of Z Problem 3 (by transforming c, the min of Prob 4)
K = ellipke(1-1/b^2);
[sn, cn, dn] = ellipj((0:k)*K/k, 1-1/b^2);
extrema = b*dn;
vals = 1-r(extrema);
c = mean( vals(1:2:end) );
e = eig( [ 2-4/c^2 1 ; -1 0 ] );
Zk = min(abs(e));

% Obtain the polynomials p and q of the rational function r = p/q 
[p,q,pq] = poly(r);
pp = [0, p];

% Transform Z Problem 4 into Problem 3 using eq. (2) - rearranged in Theorem 2.1
% (Istace/Thiran paper)and have the denominator given by
denom = q.*(1-Zk) - pp.*(1+Zk);
roots_denom = roots(denom);

% %%% work with 6 poles only
% sm6 = roots_denom(3:8); % 6 smallest poles
% la6 = roots_denom(1:6); % 6 largest poles
% s3l3 = zeros(6,1); 
% s3l3(1:3) = roots_denom(1:3);
% s3l3(4:6) = roots_denom(6:8);
% rev_s3l3 = flip(s3l3);

% %%% work with 4 poles only
% sm4 = roots_denom(5:8); % 4 smallest poles
% la4 = roots_denom(1:4); % 4 largest poles
% l2s2 = zeros(4,1); 
% l2s2(1:2) = roots_denom(1:2);
% l2s2(3:4) = roots_denom(7:8);
% s2l2 = flip(l2s2);

%%%% Other orderings of the poles
% % reversed order poles
% rev_poles = flip(roots_denom);
% % random order of poles
% rand_poles = roots_denom(randperm(length(roots_denom)));
% % specified order of poles so that smallest first, largest second, etc
% vec_srt = [8, 1, 7, 2, 6, 3, 5, 4];
% sort_poles = roots_denom(vec_srt);
%%%%

% time & solve using RKPG
tic;
[X1, X2, final_err, vec_res, it, inner_it, avg_inner, error_vec] = RKPG(A, rhs1, rhs2, roots_denom, tol,  maxit, Xex_mat);
%!! for beckermann bound need to add extra 'upper_bound' to outputs
%!! to check errors need error_vec in outputs and Xex_mat (exact solution in matrix form) in inputs
%!! plot Ritz values need e_Ap in outputs
time = toc;

fprintf('\n Total execution time: %9.4e seconds \n', time)

fprintf('final_err   avg_inner  \n')
fprintf('\n  %9.4e       %d    \n \n', [final_err, avg_inner])

% plot residual v iterations
iter = linspace(1, it, it);
semilogy(iter, vec_res, 'v');hold on
xlabel('Iterations');
ylabel('Residual');

% plot Beckermann bound on top of residuals
% semilogy(upper_bound, 'x'); hold off;