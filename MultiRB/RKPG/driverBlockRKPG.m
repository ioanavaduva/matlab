% Driver for RKPG (with 2-sided projection)

clear all;
addpath(genpath('../../rktoolbox'));

% Setup
n = 500; % size of matrix A 
h = 1/n; eps = 1;
A = eps*(diag(2*ones(n, 1)) + diag (-1*ones(n-1, 1), 1) + diag (-1*ones(n-1, 1), -1))/h^2;

% choices for RHS
% rhs1=ones(n, 1); 
% rhs1 = [ones(n, 1), A*ones(n, 1)];

% rhs1 = [ones(n,1), rand(n,1)];
% rhs1 = rand(n, 4); rhs2=rhs1;
% rhs1 = [ones(n, 1),((-1).^(0:n-1))', [ones(n/2, 1); zeros(n/2, 1)]];
% rhs1 = [[ones(n/2, 1); zeros(n/2, 1)], [zeros(n/2, 1); ones(n/2, 1)]];

% alternate 1&-1 columns 
% rhs1 = ((-1).^(1:n))';
% rhs1 = [[((-1).^(1:n/2))'; zeros(n/2, 1)], [zeros(n/2, 1); ((-1).^(1:n/2))']];

% columns of identity combined
% I = speye(n);
% rhs1 = I(:, 5);
% rhs1 = [I(:,1), I(:,1)+I(:,2)];
% rhs1 = [I(:,1), I(:,1)+I(:,2), I(:,1)+I(:,2)+I(:,3)];

% column of identity (lin indep)
% rhs1 = I(:,1);
% rhs1 = I(:, 1:2);
% rhs1 = I(:, 1:3);

% build up rand rhs
% x1 = rand(n,1);
% x2 = rand(n,1);
% x3 = rand(n,1);
% rhs1 = x1;
% rhs1 = [x1,x2];
% rhs1 = [x1,x2,x3];

% eigenvectors
% [V, ~] = eigs(A, n);
% % rhs1 = V(:, 2);
% % rhs1 = V(:, 2:3);
% rhs1 = V(:, 78:80);

% random orthonormal vectors
% Q = rand(n, 3);
% Q = orth(Q);
% rhs1 = Q(:, 1);
% rhs1 = Q(:, 1:2);
% rhs1 = Q;

rhs2 = rhs1;

tol = 1e-9;
maxit = 300;

% Get smallest and largest eigenvalues
emin = 1e-6; 
opts.tol=1e-4;
emax = eigs(A, 1,'LA',opts);

% compute roots denom from Zolotarev problem
bb = emax - emin + 1;
k = 4;      % number of poles is 2*k
roots_denom = get_rootsden(k, bb);

tic;
[X1, X2, vec_res, it, final_err, upper_vec] = RKPGblock2(A, rhs1, rhs2, roots_denom, tol,  maxit);
time = toc;

fprintf('\n Total execution time: %9.4e seconds \n', time)

fprintf('final_err \n')
fprintf('\n  %9.4e \n \n', final_err)

% plot residual v iterations
iter = linspace(1, it, it);
semilogy(iter, vec_res, 'p');hold on
xlabel('Iterations');
ylabel('Residual');
