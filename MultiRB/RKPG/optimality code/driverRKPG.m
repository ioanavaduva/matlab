%% Driver for RKPG (with 2-sided projection)

%% 2D Poisson
n = 1000; 
poisson_2d;

tol = 1e-9;
maxit = 300;

%% Diffusion coefficients
% n = 1001;
% separable_coeff;

%% Nonuniform (for graded mesh need extra code from zip; for geometric mesh code attached)
% n = 1000;
% nonuniform_mesh3; 


%% %% Get smallest and largest eigenvalues
opts.tol=1e-4;
emin = eigs(A, 1, 'smallestabs', opts);
emax = eigs(A, 1,'largestabs',opts);

%% %% Different poles

%% logspace poles
% tic;
% poles_log = logspace(log10(emin), log10(emax), 16)';
% toc;
% poles_log = sort(poles_log, 'desc');

%% Get nodes (Sabino thesis)
s_nodes = 16;                           % Choose 2 nodes (could vary)
tic;
snew = get_nodes2(emin,emax,s_nodes);      % Use interval for A_1;
t = toc;
fprintf('\n Sabino poles time: %9.4e seconds \n', t)
s_parameter = sort(snew, 'desc');

%% IRKA
m = 16; % number of poles; can change
xi = emin + (emax-emin)*rand(1,m);
tic;
[shifts, its] = irka_shifts(A, rhs1, xi, 1e-6); % IRKA poles with double orthogonalization & small tolerance
tt = toc;
fprintf('\n IRKA poles time: %9.4e seconds \n', tt)
shifts = sort(shifts, 'desc');

%% time & solve using RKPG
tic;
[X1, X2, final_err, vec_res, it, inner_it, avg_inner, Ap, rhsp] = RKPG(A, rhs1, rhs2, shifts, 1e-8, 15);
time = toc;

fprintf('\n Total execution time: %9.4e seconds \n', time)

fprintf('final_err   avg_inner  \n')
fprintf('\n  %9.4e       %d    \n \n', [final_err, avg_inner])

%% Compute H_hat and H to check optimality condition
ritz = eig(Ap);
ritz = sort(ritz, 'desc');
fprintf('\n Norm (poles-Ritz: %9.4e \n', norm(shifts-ritz))
for i = 1:16
H_hat(i) = rhsp'*(ritz(i)*eye(size(Ap)) + Ap)*rhsp;
H(i) = rhs1'*(ritz(i)*eye(size(A)) + A)*rhs1;
diffH(i) = norm(H_hat(i) - H(i)); % Optimality norm for each Ritz value
nmH(i) = norm(H(i));
nmH_hat(i) = norm(H_hat(i));
end
fprintf('\n Optimality norm: %9.4e \n',norm(H_hat - H))
