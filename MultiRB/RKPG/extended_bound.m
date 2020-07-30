% Knizhnerman & Simoncini Extended Krylov bound from "Convergence analysis
% of the extended Krylov method" eq (4.1) & (3.17)
function bound_vec = extended_bound(n, it)
% Set up of Lyapunov equation
h = 1/n; eps = 1;
A = eps*(diag(2*ones(n, 1)) + diag (-1*ones(n-1, 1), 1) + diag (-1*ones(n-1, 1), -1))/h^2;

% Need smallest eigenvalue of A (set small)
opts.tol=1e-4;
emin = 1e-6; %eigs(A, 1, 'SA', opts); %1e-6 for problems larger than 800
emax = eigs(A, 1, 'LA', opts);

kappa = emax/emin;
rho = ((kappa^(1/4)-1)/(kappa^(1/4)+1))^2;

% Plot the bound vector (initialized as zeros(iterations, 1)
bound_vec = zeros(it, 1);
for i = 1:it
er = i*rho^i; %conjectured: er = rho^n/sqrt(n)
bound_vec(i) = er;
end
plot(bound_vec, 'x');
end