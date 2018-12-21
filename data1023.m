% data for Poisson 2D problem with 15 unknowns
% function b
xtemp = [0:9.7752e-04:1];
x = repmat(xtemp, 1, 1023);
y = reshape(repmat(xtemp, length(xtemp), 1), 1, length(xtemp)^2);
b = @(x, y) sin(pi*x).*cos(pi*y);
b1023 = b(x, y)';
% initial guess x0
x0 = zeros(length(b1023), 1);
% tolerance
TOL = 10^-8;
% alowed number of iterations
maxit = 500;
% Poisson matrix
A1023 = kr_pois(1023);

TOL = 10^-8;
maxit = 500;