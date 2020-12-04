%% Driver for RKPG (with 2-sided projection)
clear all;
addpath(genpath('../../rktoolbox'));

%% Setup
n = 1000; % size of matrix A 
h = 1/n; eps = 1;
A = eps*(diag(2*ones(n, 1)) + diag (-1*ones(n-1, 1), 1) + diag (-1*ones(n-1, 1), -1))/h^2;
tol = 1e-9;
maxit = 300;

%% RHS choices -- [1, A*1, ...]
% rhs1 = ones(n, 1); 
% rhs1 = [ones(n, 1), A*ones(n, 1)];
% rhs1 = [ones(n, 1), A*ones(n, 1), A^2*ones(n, 1), A^3*ones(n, 1)];

%% RHS choices -- [alt, A*alt, ...]
% rhs = ((-1).^(1:n))';
% rhs1 = [rhs, A*rhs];
% % rhs1 = [rhs, A*rhs, A^2*rhs];
% rhs1 = [rhs, A*rhs, A^2*rhs, A^3*rhs];

%% RHS choices -- [rand, A*rand, ...]
% rhs1 = rand(n, 1);
% r = rhs1;
% rhs1 = [r, A*r];
% rhs1 = [r, A*r, A^2*r];
% rhs1 = [r, A*r, A^2*r, A^3*r];

%% RHS choices -- [cos, A*cos, ...]
% x = linspace(0, 1, n)';
% rhs = cos(pi*x);
% % rhs1 = [rhs, A*rhs];
% % % rhs1 = [rhs, A*rhs, A^2*rhs];
% rhs1 = [rhs, A*rhs, A^2*rhs, A^3*rhs];

%% RHS choices -- other combinations of columns
% rhs1 = rand(n, 4); rhs2=rhs1;
% rhs1 = [ones(n, 1),((-1).^(0:n-1))', [ones(n/2, 1); zeros(n/2, 1)]];
% rhs1 = [[ones(n/2, 1); zeros(n/2, 1)], [zeros(n/2, 1); ones(n/2, 1)]];

%% RHS choices -- alternate 1&-1 columns 
% rhs1 = ((-1).^(1:n))';
% rhs1 = [[((-1).^(1:n/2))'; zeros(n/2, 1)], [zeros(n/2, 1); ((-1).^(1:n/2))']];

%% RHS choices -- columns of identity combined
% I = eye(n);
% rhs1 = I(:, 5);
% rhs1 = [I(:,1), I(:,1)+I(:,2)];
% rhs1 = [I(:,1), I(:,1)+I(:,2), I(:,1)+I(:,2)+I(:,3)];

%% RHS choices -- column of identity (lin indep)
% rhs1 = I(:,1);
% rhs1 = I(:, 1:2);
% rhs1 = I(:, 1:3);

%% RHS choices -- build up rand rhs
% x1 = rand(n,1);
% x2 = rand(n,1);
% x3 = rand(n,1);
% rhs1 = x1;
% rhs1 = [x1,x2];
% rhs1 = [x1,x2,x3];

%% RHS choices -- eigenvectors
% [V, ~] = eigs(A, n);
% % rhs1 = V(:, 2);
% % rhs1 = V(:, 2:3);
% rhs1 = V(:, 78:80);

%% RHS choices -- random orthonormal vectors
% Q = rand(n, 3);
% Q = orth(Q);
% rhs1 = Q(:, 1);
% rhs1 = Q(:, 1:2);
% rhs1 = Q;

%% RHS choices -- purely random columns
rhs1 = rand(n, 3);

%% RHS choices -- two identical columns
% rhs1 = [ones(n, 1), 3*ones(n, 1), rand(n, 1)];

%% rhs factors are the same
rhs2 = rhs1;

%% Get smallest and largest eigenvalues
opts.tol = 1e-4;
emin = eigs(A, 1,'smallestabs',opts); 
emax = eigs(A, 1,'largestabs',opts);

%% Zolotarev shifts
bb = emax - emin + 1;
k = 4;      % number of poles is 2*k
[roots_denom, ~] = get_rootsden(k, bb);
roots_denom = roots_denom + emin - 1;
%% IRKA shifts
m = 12; % number of poles; can change
xi = emin + (emax - emin)*rand(1,m);
B = rhs1*rhs2';
B_hat = tangential_dir(A, B, xi);
[shifts, its] = irka_shifts_block(A, B, B_hat, xi, tol);

%% Solver
tic;
[X1, X2, vec_res, it, final_err] = RKPGblock2(A, rhs1, rhs2, shifts, tol,  maxit);
time = toc;

fprintf('\n Total execution time: %9.4e seconds \n', time)

fprintf('final_err \n')
fprintf('\n  %9.4e \n \n', final_err)

%% plot residual v iterations
iter = linspace(1, it, it);
semilogy(iter, vec_res, 'x');hold on
xlabel('Iterations');
ylabel('Residual');
