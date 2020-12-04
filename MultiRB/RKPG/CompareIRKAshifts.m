clear all;
addpath(genpath('../../rktoolbox'));

%% Set up
n = 100; % can change
h = 1/n;
A = full(gallery('tridiag',n,-1,2,-1))/h^2;
In = eye(n);
% eps = 1;
A_rksm = -((diag(2*ones(n, 1)) + diag(-1*ones(n-1, 1), 1) + diag(-1*ones(n-1, 1), -1))/h^2);

% density = 2/n; % for example
% rc = 0.1; % Reciprocal condition number
% A = sprandsym(n, density, rc, 1);

%% %% RHS
%% Ones
f = ones(n,1); 
%% Rand
% f = randn(n, 1);
%% Linspace
% x = linspace(0, 1, n)';
%% Cos(pi*x)
% f = cos(pi*x);
%% Alternate
% f = ((-1).^(0:n-1))';

%% Polynomial
% xtemp = linspace(0,1,n);
% x = repmat(xtemp, 1, n);
% y = reshape(repmat(xtemp, length(xtemp), 1), 1, length(xtemp)^2);
% b = @(x, y) -2.*(6.*x.^2 - 6.*x + 1).*(y-1).^2.*y.^2-2.*(x-1).^2.*x.^2.*(6.*y.^2-6.*y+1);
% rhs = b(x, y)';
% rhs11 = reshape(rhs, n, n); 
% f = rhs11(:,1);

%% Sin(pi*x)Sin(pi*y)
% xtemp = linspace(0,1,n);
% x = repmat(xtemp, 1, n);
% y = reshape(repmat(xtemp, length(xtemp), 1), 1, length(xtemp)^2);
% % b = @(x, y) sin(pi.*x).*sin(pi.*y); % without differentiating -- gives NaNs
% b = @(x, y) pi^2*sin(pi.*x).*sin(pi.*y) + pi^2*sin(pi.*x).*sin(pi.*y);
% f = sqrt(2)*pi*sin(pi.*xtemp)';
% % rhs = b(x, y)';
% % rhs11 = reshape(rhs, n, n); 
% % f = rhs11(:,1);

%% Sin(pi*x)Cos(pi*x)
% xtemp = linspace(0,1,n);
% x = repmat(xtemp, 1, n);
% y = reshape(repmat(xtemp, length(xtemp), 1), 1, length(xtemp)^2);
% b = @(x, y) sin(pi.*x).*cos(pi.*y); % without differentiating -- I think satisfy bc
% % b = @(x, y) pi^2*sin(pi.*x).*cos(pi.*y) - pi^2*sin(pi.*x).*cos(pi.*y); % this is just 0 
% rhs = b(x, y)';
% rhs11 = reshape(rhs, n, n); 
% f = rhs11(:,1);

%% Plot eigenvalues 
% Re = real(eig(A));
% Im = imag(eig(A));
% 
% semilogx(Re, Im, 'o'); hold on;

%% Compute spectral interval
lam_max = max(eig(A));
lam_min = min(eig(A));
opts.tol = 1e-4;
emax_rksm = eigs(A_rksm, 1,'largestabs',opts);
emin_rksm = eigs(A_rksm, 1, 'smallestabs', opts);


%% %% Pole Choices

%% Ioana's roots denominator/Zolotarev
% bb = lam_max - lam_min + 1;  
% k = 4; % number of poles from this is 2*k 
% roots_denom = get_rootsden(k, bb) + lam_min - 1;

% re = real(roots_denom);
% im = imag(roots_denom);
% plot(re, im, 's'); hold on;

%% IRKA
m = 16; % number of poles; can change
xi = lam_min + (lam_max-lam_min)*rand(1,m);
tic;
[shifts2, it2] = irka_shifts(A, f, xi, 1e-4);
toc;
% r = real(shifts);
% i = imag(shifts);
% plot(r, i, 'p'); hold on;

%% Adaptive rksm.m
% [Z, nrmrestot, pol] = rksm(A_rksm, In, In, f, 300, 1e-9, emin_rksm, emax_rksm, 0, 1e-12);
% % typically s_0^(1) approximates the part of the spectrum closest to zero, 
% % and s_0^(2) that farthest from zero
% [pols,m1,n1] = uniquetol(pol);
% [c1,d1] = sort(m1);
% pols = pols(d1);
% R = real(pols);
% I = imag(pols);
% plot(R, I, 'o'); hold on;
% fprintf('\n Number of rksm unique poles: %d \n', length(pols))
% fprintf('\n Number of rksm poles (total): %d \n', length(pol))