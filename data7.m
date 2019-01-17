% data for Poisson 2D problem with 31 unknowns
% function b
xtemp = [0:0.1429:1];
x = repmat(xtemp, 1, 7);
y = reshape(repmat(xtemp, length(xtemp), 1), 1, length(xtemp)^2);
b = @(x, y) sin(pi*x).*cos(pi*y);
b7 = b(x, y)';
% initial guess x0
x0 = zeros(length(b7), 1);
% tolerance
TOL = 10^-8;
% alowed number of iterations
maxit = 500;
% Poisson matrix
A7 = kr_pois(7);

TOL = 10^-8;
maxit = 500;