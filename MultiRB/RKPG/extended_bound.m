% Knizhnerman & Simoncini Extended Krylov bound from "Convergence analysis
% of the extended Krylov method" eq (4.1) & (3.17)

% Set up of Lyapunov equation
n = 1000; 
h = 1/n; eps = 1;
A = eps*(diag(2*ones(n, 1)) + diag (-1*ones(n-1, 1), 1) + diag (-1*ones(n-1, 1), -1))/h^2;

% Need smallest eigenvalue of A (set small)
opts.tol=1e-4;
emin = 1e-6; %eigs(A, 1, 'SA', opts); %1e-6 for problems larger than 800
emax = eigs(A, 1, 'LA', opts);

kappa = emax/emin;
rho = ((kappa^(1/4)-1)/(kappa^(1/4)+1))^2;

% Plot the bound vector (initialized as zeros(iterations, 1)
bound_vec = zeros(58, 1);
for n = 1:57
er = n*rho^n; %conjectured: er = rho^n/sqrt(n)
bound_vec(n+1) = er;
end
plot(bound_vec, 'x');