% Approx large IRKA poles by small problem poles

%% Large problem poles
n = 1001;
separable_coeff;

%% Poisson 
% n = 1000; h = 1/n; 
% A = (diag(2*ones(n, 1)) + diag(-1*ones(n-1, 1), 1) + diag(-1*ones(n-1, 1), -1))/h^2;

% rhs1 = ones(n, 1);

% rhs1 = ((-1).^(0:n-1))';

% xtemp = linspace(0,1,n);
% rhs1 = sqrt(2)*pi*sin(pi.*xtemp)';
% rhs2 = rhs1;

% xtemp = linspace(0,1,n);
% x = repmat(xtemp, 1, n);
% y = reshape(repmat(xtemp, length(xtemp), 1), 1, length(xtemp)^2);
% b = @(x, y) -2.*(6.*x.^2 - 6.*x + 1).*(y-1).^2.*y.^2-2.*(x-1).^2.*x.^2.*(6.*y.^2-6.*y+1);
% rhs = b(x, y)';
% rhs11 = reshape(rhs, n, n); rhs1 = rhs11(:,1);

% xtemp = linspace(0,1,n);
% rhs1 = sqrt(2)*pi*sin(pi.*xtemp)';
% rhs2 = rhs1;

% x = linspace(0, 1, n)';
% rhs1 = cos(pi*x); rhs2 = rhs1;

%% Prep
opts.tol=1e-4;
emin = eigs(A, 1, 'smallestabs', opts);
emax = eigs(A, 1,'largestabs',opts);

m = 16;
xi = emin + (emax-emin)*rand(1,m);
tic;
[shifts, its] = irka_shifts(A, rhs1, xi, 1e-4);
toc;

%% Small problem poles -- unscaled

n_small = floor(n/10); 
h_small = 1/n_small;

%% Diffusion coefficients small

g_small = 0:1/n_small:1; 
gs_small = 1/(2*n_small):1/n_small:1;

% a = @(x) 1+x; 
a = @(x) sin(x);
% a = @(x) 1 +  50*exp(-5*x.^2);
% syms x
% a = @(x) ((0<x) & (x<0.5)).*1 + ((0.5<x) & (x<1)).*1e+4;

ag_small = a(g_small(2:n_small));
ags_small = a(gs_small);

D_small = diag(ag_small);
D_inv_small = diag(sqrt(1./ag_small));

T_orig_small = gallery('tridiag', -ags_small(2:n_small-1), ags_small(1:n_small-1)+ ags_small(2:n_small), -ags_small(2:n_small-1))/h_small^2;

A_small = D_inv_small*T_orig_small*D_inv_small;


%% Poisson small
% A_small = (diag(2*ones(n_small, 1)) + diag(-1*ones(n_small-1, 1), 1) + diag(-1*ones(n_small-1, 1), -1))/h_small^2;

%% RHS small
rhs1_small = ones(n_small-1, 1); 
rhs1_small = D_inv_small*rhs1_small;

% rhs1_small = ((-1).^(0:n_small-2))';
% rhs1_small = D_inv_small*rhs1_small;

% xtemp_small = linspace(0,1,n_small-1);
% x = repmat(xtemp_small, 1, n_small-1);
% y = reshape(repmat(xtemp_small, length(xtemp_small), 1), 1, length(xtemp_small)^2);
% b = @(x, y) -2.*(6.*x.^2 - 6.*x + 1).*(y-1).^2.*y.^2-2.*(x-1).^2.*x.^2.*(6.*y.^2-6.*y+1);
% rhs = b(x, y)';
% rhs11_small = reshape(rhs, n_small-1, n_small-1); rhs1_small = rhs11_small(:,1);
% rhs1_small = D_inv_small*rhs1_small;

% xtemp_small = linspace(0,1,n_small-1);
% rhs1_small = sqrt(2)*pi*sin(pi.*xtemp_small)';
% rhs1_small = D_inv_small*rhs1_small;

% x_small = linspace(0, 1, n_small-1)';
% rhs1_small = cos(pi*x_small); 
% rhs1_small = D_inv_small*rhs1_small;

%% Spectral interval 

% emin_small = eigs(A_small, 1, 'smallestabs', opts);
% emax_small =  eigs(A_small, 1, 'largestabs', opts);

%% Small shifts
tic;
[shifts_small, its_small] = irka_shifts(A_small, rhs1_small, xi, 1e-4);
toc;

%% Small problem poles -- scaled
% 
% A_scaled = (diag(2*ones(n_small, 1)) + diag(-1*ones(n_small-1, 1), 1) + diag(-1*ones(n_small-1, 1), -1))/h^2;
% rhs1_small = ones(n_small, 1);
% 
% % emin_scaled = eigs(A_scaled, 1, 'smallestabs', opts);
% % emax_scaled =  eigs(A_scaled, 1, 'largestabs', opts);
% 
% % m_scaled = m-2; % for poles 2
% % m_scaled = m-5; % for poles 1
% % xi_scaled = emin + (emax-emin)*rand(1,m_scaled);
% tic;
% [shifts_scaled, its_scaled] = irka_shifts(A_scaled, rhs1_small, xi, 1e-4);
% toc;

%% Small problems poles combinations

% 1&2 small, 2&3*10 small scaled
poles = sort([shifts_small(length(shifts_small)/2+1:length(shifts_small)); shifts_small(length(shifts_small)/4+1:3*length(shifts_small)/4)*10], 'desc');

% 1&2 small, 3&4*10 small 
% poles = sort([shifts_small(length(shifts_small)/2+1:length(shifts_small)); shifts_small(1:length(shifts_small)/2)*10]);

% poles1 = sort([shifts_small(1); shifts_small(length(shifts_small)-3:length(shifts_small)); shifts_scaled], 'descend');
% poles2 = sort([shifts_small(1); shifts_small(length(shifts_small)-3:length(shifts_small)); shifts_scaled(4:length(shifts_scaled))], 'descend');