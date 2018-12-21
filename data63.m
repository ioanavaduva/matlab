% data for Poisson 2D problem with 31 unknowns
% function b
xtemp = [0:0.0159:1];
x = repmat(xtemp, 1, 63);
y = reshape(repmat(xtemp, length(xtemp), 1), 1, length(xtemp)^2);
b = @(x, y) sin(pi*x).*cos(pi*y);
b63 = b(x, y)';
% initial guess x0
x0 = zeros(length(b63), 1);
% tolerance
TOL = 10^-8;
% alowed number of iterations
maxit = 500;
% Poisson matrix
A63 = kr_pois(63);

TOL = 10^-8;
maxit = 500;