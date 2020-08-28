addpath(genpath('../../rktoolbox'));

% Set up
n = 500; % can change
h = 1/n;
A = full(gallery('tridiag',n,-1,2,-1))/h^2;

% density = 2/n; % for example
% rc = 0.1; % Reciprocal condition number
% A = sprandsym(n, density, rc, 1);

f = ones(n,1); % can change
% f = randn(n, 1);
% x = linspace(0, 1, n)';
% f = cos(pi*x);
% f = ((-1).^(0:n-1))';
In = eye(n);

Re = real(eig(A));
Im = imag(eig(A));

semilogx(Re, Im, 'o'); hold on;

% Compute spectral interval
lam_max = max(eig(A));
lam_min = min(eig(A));

% %%% Random poles from spectral interval
m = 12; % number of poles; can change
xi = lam_min + (lam_max-lam_min)*rand(1,m);
%%%

%%% Ioana's roots denominator
bb = lam_max - lam_min + 1; 
b = bb; 
k = 4; % number of poles from this is 2*k 
roots_denom = get_rootsden(k, b);
%%%
re = real(roots_denom);
im = imag(roots_denom);
plot(re, im, 'x'); hold on;

% [shifts, it] = irka_shifts(A, f, roots_denom, 1e-4);
[shifts, it] = irka_shifts(A, f, xi, 1e-4);
r = real(shifts);
i = imag(shifts);
plot(r, i, 'p');
