% Simoncini Druskin Polynomial bound for error from Proposition 3.1
function bound_vec = polynomial_bound(n, it)
% Set up of Lyapunov equation 
h = 1/n; eps = 1;
A = eps*(diag(2*ones(n, 1)) + diag (-1*ones(n-1, 1), 1) + diag (-1*ones(n-1, 1), -1))/h^2;

% Need smallest eigenvalue of A (set small)
opts.tol=1e-4;
emin = eigs(A, 1, 'SA', opts); %1e-6 for problems larger than 800

% Need the extrem eigenvalues of A+eminI
mat = A + emin * eye(n);
lambda_min = eigs(mat, 1, 'SA', opts); 
lambda_max = eigs(mat, 1,'LA',opts);

kappa = lambda_max/lambda_min;

% Compute rhs of eq 3.1
fr1 = (sqrt(kappa)+1)/(lambda_min * sqrt(kappa)); % first fraction 
fr2 = (sqrt(kappa)-1)/(sqrt(kappa)+1); % second fraction (to be raised to the iteration number

% Plot the bound vector (initialized as zeros(iterations, 1)
bound_vec = zeros(it, 1);
for i = 1:it
er = fr1 * fr2^i;
bound_vec(i) = er;
end
plot(bound_vec, 'x');
end