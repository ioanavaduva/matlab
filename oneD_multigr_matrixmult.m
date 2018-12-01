% Multigrid for Poisson problem with 7 unknowns

% create vector b
x1 = [0.1:0.1429:1];
b = @(x) sin(pi*x).*cos(pi*x);
b_1d = b(x1)';

% create matrix T of coefficients (fine grid matrix)
T = oned_pois(7);

% pre-smoothing with damped Jacobi (do 50 iterations only)
x1 = damped_jacobiM(2/3, zeros(length(b_1d), 1), T, b_1d, 10^-7, 30);

% compute the residual
res = b_1d - T*x1;

% generate restriction matrix
n = 2^3-1;
k = log2(n+1);

N = 2^(k-1)-1;

RE = zeros(N,n);

for i = 1:N
   RE(i,2*i-1:2*i+1) = [1 2 1]; 
end

RE = RE/4;

% generate interpolation matrix
II = 2*RE';

% transfer residual to coarse grid; v is coarse grid residual
v = RE*res;

% transfer matrix T to coarse grid
TC = RE * T * II;

% solve residual equation to find error
err = TC\v;

% transfer error to fine grid; r is fine grid error
erf = II*err;


% correct approximation (initial guess for damped Jacobi)
x1 = x1 + erf;

% post-smoothing Jacobi (50 iterations)
x1 = damped_jacobiM(2/3, x1, T, b_1d, 10^-7, 30);
