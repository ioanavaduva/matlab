clear;

addpath(genpath('../../../rktoolbox'))

n = 1000;

h = 1/n; ep = 0.0083;

A = ep*(spdiags([-ones(n, 1) 2*ones(n, 1) -ones(n, 1)],-1:1,n,n))/(h^2);
B = spdiags([-ones(n, 1) zeros(n, 1) ones(n, 1)],-1:1,n,n)/(2*h);
% B = zeros(n, n);
% separable_coeff;
% nonuniform_mesh3;

% X = randn(n); D = diag(rand(n, 1)); M = X*D*inv(X);
% Y = randn(n); D = diag(rand(n, 1)); N = Y*D*inv(Y);

M = A + B;
N = A + B';

rhs1 = ones(n, 1);
% rhs2 = randn(n, 1);
rhs2 = rhs1;

opts.tol = 1e-4;
emin = eigs(M, 1, 'smallestabs', opts);
emax = eigs(M, 1, 'largestabs',opts);

% s_nodes = 20;        
% tic;
% snew = get_nodes2(emin,emax,s_nodes);      
% Zolo_time = toc;
% s_parameter = sort(snew, 'desc');
% fprintf('\n Zolotarev poles time: %9.4e seconds \n', Zolo_time)

% bb = emax - emin + 1;
% k = 8;      % number of poles is 2*k
% tic;
% [roots_denom, ~] = get_rootsden(k, bb);
% rk_time = toc;
% roots_denom = roots_denom + emin - 1;
% roots_denom = sort(roots_denom, 'desc');
% fprintf('\n RK Toolbox poles time: %9.4e seconds \n', rk_time)

m = 16; % number of poles, can change
xi = emin + (emax-emin)*rand(1,m); 
tic;
[shifts, its] = irka_shifts(M, rhs1, xi, 1e-2);
original_shifts_time = toc;
fprintf('\n Original poles time: %9.4e seconds \n', original_shifts_time)
 
% tic;
% [shifts1, shifts2, its2] = irka_shifts2(M, N, rhs1, rhs2, xi,  xi, 1e-2); 
% shifts2_time = toc;
% fprintf('\n Shifts2 time: %9.4e seconds \n', shifts2_time)
% shifts1 = sort(shifts1, 'desc');
% shifts2 = sort(shifts2, 'desc');

%% Solve 
tic;
[X1, X2, final_err, vec_res, it] = RKPG_S(M, N, rhs1, rhs2, shifts, 1e-8, 300);
time = toc;

fprintf('\n Total execution time: %9.4e seconds \n', time)

fprintf(' Final residual:  %9.4e   \n \n', final_err)

%% plot residual vs. iterations
% iter = linspace(1, it, it);
% semilogy(iter, vec_res, 'p'); hold on
% xlabel('Iterations');
% ylabel('Residual');

