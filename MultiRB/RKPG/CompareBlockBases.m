% Compare block bases

% Set up
n = 10; 
A = full(gallery('tridiag',n,-1,2,-1));
f = [ones(n,1), rand(n,1)];
In = eye(n);

% Compute spectral interval
lam_max = max(eig(A));
lam_min = min(eig(A));

%%% Random poles from spectral interval
m = 2; % number of poles
xi = lam_min + (lam_max-lam_min)*rand(1,m);
%%%

%Stefan's basis
V = rat_krylov(A,f,-xi);

% Ioana's basis
W = get_start_basis(f);
for k = 1:length(f(1, :)) % k refers to what column we're filling in
    for i = 1:length(xi)
        W = get_rk_block_basis(A, xi(i), W, k);
    end
end

fprintf('Bases norm after all poles = %g\n',norm(V - W));
fprintf('Angle between bases = %g\n', subspace(V,W));
