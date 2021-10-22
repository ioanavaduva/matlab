clear all;

addpath(genpath('../../../rktoolbox'))

n = 1000;

h = 1/n; ep = 0.0167;

A = ep*(spdiags([-ones(n, 1) 2*ones(n, 1) -ones(n, 1)],-1:1,n,n))/(h^2);
B = spdiags([-ones(n, 1) zeros(n, 1) ones(n, 1)],-1:1,n,n)/(2*h);

% B = zeros(n, n);
% separable_coeff;
% nonuniform_mesh3;

% X = randn(n); D = diag(rand(n, 1)); M = X*D*inv(X);
% Y = randn(n); D = diag(rand(n, 1)); N = Y*D*inv(Y);

% convection coefficients
xtemp = linspace(0, 1, n);
% w1 = @(x) 1 + ((x+1).^2)/4;
% phi = diag(w1(xtemp));

w1 = @(x) 1 - (2*x + 1).^2;
phi = diag(w1(xtemp));

w2 = @(x) 1 - x.^2;
psi = diag(w2(xtemp));

% w1 = @(x) -2*(2*x + 1);
% phi = diag(w1(xtemp));
% 
% w2 = @(x) x;
% psi = diag(w2(xtemp));

% M = A + B;
% N = A + B';

% M = A + phi*B; 
% N = A;

M = A + phi*B;
N = A + B*psi;

rhs1 = ones(n, 1);
% rhs2 = randn(n, 1);
rhs2 = rhs1;

opts.tol = 1e-2;
opts.p = 30;
opts.maxit = 10000;
% emin = eigs(M, 1, 'smallestabs', opts);
% emax = eigs(M, 1, 'largestabs',opts);
% 
% emin_N = eigs(N, 1, 'smallestabs', opts);
% emax_N = eigs(N, 1, 'largestabs', opts);

for j = 1:n
    evi(j) = M(1,1) + 2*sqrt(M(1, 2) * M(2, 1)) * cos(j*pi/(n+1));
end
evi = sort(evi, 'desc');
emin = min(evi);
emax = max(evi);

s_nodes = 12;        
tic;
snew = get_nodes2(emin,emax,s_nodes);    
% snew2 = get_nodes2(emin_N, emax_N, s_nodes);
Zolo_time = toc;
s_parameter = sort(snew, 'desc');
% s_parameter2 = sort(snew2, 'desc');
fprintf('\n Sabino poles time: %9.4e seconds \n', Zolo_time)

% bb = emax - emin + 1;
% k = 4;      % number of poles is 2*k
% tic;
% [roots_denom, ~] = get_rootsden(k, bb);
% bb_N = emax_N - emin_N + 1;
% [roots_denom2, ~] = get_rootsden(k, bb_N);
% rk_time = toc;
% roots_denom = roots_denom + emin - 1;
% roots_denom2 = roots_denom2 + emin_N -1;
% roots_denom = sort(roots_denom, 'desc');
% roots_denom2 = sort(roots_denom2, 'desc');
% fprintf('\n RK Toolbox poles time: %9.4e seconds \n', rk_time)


% m = 20; % number of poles, can change
% xi = emin + (emax-emin)*rand(1,m); 
% tic;
% [shifts, its] = irka_shifts(M, rhs1, xi, 1e-2);
% shifts = sort(shifts, 'desc');
% original_shifts_time = toc;
% fprintf('\n Original IRKA poles time: %9.4e seconds \n', original_shifts_time)

% xi_N = emin_N + (emax_N - emin_N)*rand(1, m);
% tic;
% [shifts1, shifts2, its2] = irka_shifts2(M, N, rhs1, rhs2, xi,  xi_N, 1e-2); 
% shifts2_time = toc;
% fprintf('\n Shifts2 time: %9.4e seconds \n', shifts2_time) 
% shifts1 = sort(shifts1, 'desc');
% shifts2 = sort(shifts2, 'desc');

%% Solve 
tic;
[X1, X2, final_err, vec_res, it] = RKPG_S(M, N, rhs1, rhs2, s_parameter, 1e-8, 300);
% [X1, X2, final_err, vec_res, it] = RKPG_S2(M, N, rhs1, rhs2, roots_denom, roots_denom, 1e-8, 300);

time = toc;

fprintf('\n Total execution time: %9.4e seconds \n', time)

fprintf(' Final residual:  %9.4e   \n \n', final_err)



%% plot residual vs. iterations
iter = linspace(1, it, it);
semilogy(iter, vec_res, 'p'); hold on
xlabel('Iterations');
ylabel('Residual');

