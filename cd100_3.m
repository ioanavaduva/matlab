% Convection-diffusion problem with wind as in example 3 in D.Pallita & V.
% Simoncini -- Matrix -equation-based strategies for convection-diffusion
% equation, with zero Dirichlet boundary conditions

% 2D convection - diffusion problem with 15 unknowns
n = 100;
h = 1/n;

% Poisson Matrix (needs multiplied by an epsilon between 0 and 1).
A = kr_pois(n);

% Right hand side vector b coming from function b(x, y)= sin(pi*x)*cos(pi*y)
xtemp = [0:0.0101:1];
x = repmat(xtemp, 1, n);
y = reshape(repmat(xtemp, length(xtemp), 1), 1, length(xtemp)^2);
b = @(x, y) sin(pi*x).*cos(pi*y);
b100 = b(x, y)';

% Create B matrix for convection part of the equation
B = zeros(n^2, n^2);
B = sparse(diag(ones(n^2-1, 1), 1)+diag(-1*ones(n^2-1, 1), -1));

% Generate convection vector w=(w1, w2), where w1(x, y) = ph1(x)*ps1(y) &
% w2 (x, y) = ph2(x)*ps2(y)
ph1 = @(x) 1-x.^2;
ph2 = @(x) -2*x;
ps1 = @(y) 2*y;
ps2 = @(y) 1-y.^2;

% Create diagonal matrices 
Phi1 = sparse(diag(ph1(x)));
Phi2 = sparse(diag(ph2(x)));
Psi1 = sparse(diag(ps1(y)));
Psi2 = sparse(diag(ps2(y)));

% Get final matrix C-D for epsilon = 0.0333
eps = 0.0083;
CD = eps*A/(h^2) + Phi1*B*Psi1/(2*h) + Phi2*B'*Psi2/(2*h);


