% DKS bound from Theorem 4.8
function bound_vec = DKSbound(n, it)

% Set up of Lyapunov equation 
h = 1/n; eps = 1;
A = eps*(diag(2*ones(n, 1)) + diag (-1*ones(n-1, 1), 1) + diag (-1*ones(n-1, 1), -1))/h^2;

% First, we need to transform the spectrum into [delta, 1], with
% 0 < delta < 1. We divide through by emax to get delta = emin/emax
emin = 1e-6; 
opts.tol=1e-4;
emax = eigs(A, 1,'LA',opts);

delta = emin/emax;

% Compute coefficients a, b, c, d and modulus k
a = -(3 + delta)/(1 - delta) - sqrt(((3 + delta)/(1 - delta))^2 - 1);
b = -(1 + 3*delta)/(1 - delta) - sqrt(((1 + 3*delta)/(1 - delta))^2 - 1);
c = 1/b;
d = 1/a;

k = sqrt(((d - a)*(c - b))/((d - b)*(c - a)));

% Compute principal and complementary elliptic integrals of modulus k (they
% are both computed as being 'of the first kind'
m = k^2; % from documentation for ellipke & ellipCK
princ = ellipke(m);
compl = ellipticCK(m);

% Compute the right-hand side of equation (4.26) -- this needs raised to
% the power p at the pth iteration in order to bound ||X-Xp||.
e = exp(-(pi*princ)/(2*compl)); 

% Plot the bound vector (initialized as zeros(number iterations, 1)
bound_vec = zeros(it, 1);
for n = 1:it
    er = e^n;
bound_vec(n) = er;
end
plot(bound_vec, 'x');
end

% % Obtain exact X using Kronecker product
% AA = kron(A, eye(n))+kron(eye(n), A);
% rhss = ones(n^2, 1);
% Xex = AA\rhss;
% Xex_mat = reshape(Xex, n, n);
% Xapp = X1*X2; % X1 & X2 are computed from running driverRKPG with same n
% norm(Xex_mat - Xapp);
% % To find the bound take final iteration number and find e^p?
% % For n=10, RKPG take 4 iterations, norm(Xex_mat - Xapp) = 2.8700e-16 and
% % e^4 = 0.1550