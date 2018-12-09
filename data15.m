% data for Poisson 2D problem with 15 unknowns
% function b
xtemp = [0:0.0667:1];
x = repmat(xtemp, 1, 15);
y = reshape(repmat(xtemp, length(xtemp), 1), 1, length(xtemp)^2);
b = @(x, y) sin(pi*x).*cos(pi*y);
b15 = b(x, y)';
% initial guess x0
x0 = zeros(length(b15), 1);
% tolerance
TOL = 10^-8;
% alowed number of iterations
maxit = 500;
% Poisson matrix
A15 = kr_pois(15);

TOL = 10^-8;
maxit = 500;