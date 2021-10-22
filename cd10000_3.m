% Convection-diffusion problem with wind as in example 3 in D.Pallita & V.
% Simoncini -- Matrix -equation-based strategies for convection-diffusion
% equation, with zero Dirichlet boundary conditions

% 2D convection - diffusion problem with 15 unknocwns
n = 1000;
h = 1/n;

% Poisson Matrix (needs multiplied by an epsilon between 0 and 1).


% Right hand side vector b coming from function b(x, y)= sin(pi*x)*cos(pi*y)
xtemp = linspace(0,1,n);
x = repmat(xtemp, 1, n);
% y = reshape(repmat(xtemp, length(xtemp), 1), 1, length(xtemp)^2);
% b = @(x, y) sin(pi*x).*cos(pi*y);
% b100 = b(x, y)';
b1000 = ones(n^2, 1);

% Create B matrix for convection part of the equation
% v = ones(n^2, 1);
% B = spdiags([v -v], [1, -1], n^2, n^2);

% Generate convection vector w=(w1, w2), where w1(x, y) = ph1(x)*ps1(y) &
% w2 (x, y) = ph2(x)*ps2(y)
ph1 = @(x) 1-(2*x + 1).^2;
ps2 = @(y) 1-x.^2;

% Create diagonal matrices 
Phi1 = sparse(diag(ph1(x)));
Psi2 = sparse(diag(ps2(x)));

% Get final matrix C-D for epsilon = 0.0333
eps = 0.0167;
CD = eps*Acal/(h^2) + Phi1*B/(2*h) + Phi2*B'*Psi2/(2*h);