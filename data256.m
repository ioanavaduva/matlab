% data for Poisson problem with 128 unknowns
% function b
xtemp = [0:0.00392:1];
x = repmat(xtemp, 1, 256);
y = reshape(repmat(xtemp, length(xtemp), 1), 1, length(xtemp)^2);
b = @(x, y) sin(pi*x).*cos(pi*y);
b256 = b(x, y)';
% initial guess x0
x0 = zeros(length(b256), 1);
% tolerance
TOL = 10^-8;
% alowed number of iterations
maxit = 2000;
% Poisson matrix
A256 = kr_pois(256);

TOL = 10^-8;
