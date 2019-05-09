% 2D convection - diffusion problem for example 3 in paper

n = 50;
h = 1/n;

% Poisson Matrix (needs multiplied by an epsilon between 0 and 1).
A = kr_pois(n);

% Right hand side vector b coming from function b(x, y)= sin(pi*x)*cos(pi*y)
xtemp = linspace(0,1,n);
x = repmat(xtemp, 1, n);
y = reshape(repmat(xtemp, length(xtemp), 1), 1, length(xtemp)^2);
b = @(x, y) sin(pi*x).*cos(pi*y);
b100 = b(x, y)';

% Create B matrix for convection part of the equation
v = ones(n^2, 1);
B = spdiags([v -v], [1, -1], n^2, n^2);

% Generate convection vector w=(w1, w2), where w1(x, y) = ph1(x)*ps1(y) &
% w2 (x, y) = ph2(x)*ps2(y)
ph1 = @(x) 1-x.^2;
ph2 = @(x) -2*x;
ps1 = @(y) 2*y;
ps2 = @(y) 1-y.^2;

% Create diagonal matrices 
px1 = ph1(x);
px2 = ph2(x);
py1 = ps1(y);
py2 = ps2(y);
Phi1 = spdiags(px1(:), 0, n^2, n^2);
Phi2 = spdiags(px2(:), 0, n^2, n^2);
Psi1 = spdiags(py1(:), 0, n^2, n^2);
Psi2 = spdiags(py2(:), 0, n^2, n^2);

% Get final matrix C-D for epsilon = 0.0333
eps = 0.0333;
CD = eps*A/(h^2) + Phi1*B*Psi1/(2*h) + Phi2*B'*Psi2/(2*h);