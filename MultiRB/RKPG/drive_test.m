% Set-up by Jen (see OneTwoSided.m)

n = 10; k = 4;
A = full(gallery('tridiag',n,-1,2,-1));
f = ones(n,1);

In = eye(n); Ik = eye(k);

% Choose pole set
lam_max = max(eig(A));
lam_min = min(eig(A));
xi = lam_min + (lam_max-lam_min)*rand(1,k);

% Solve full problem
F = f*f'; b = F(:);
B = kron(In,A) + kron(A,In);
x = B\b; X = reshape (x,n,n);
R = A*X + X*A - F;

% Create rational krylov basis using rktoolbox
% W = rat_krylov(A,f,xi);

%Compute rational Krylov basis using RKPG
tol = 1e-9;
maxit = 300;

% [X1, X2, final_err, vec_res, it, inner_it, avg_inner, error_vec] = RKPG(A, f, f, xi, tol, maxit, X); %this keyboards to give basis V

[X1, X2, final_err, vec_res, it, inner_it, avg_inner, error_vec] = oneside_RKPG(A, f, f, xi, tol, maxit, X);