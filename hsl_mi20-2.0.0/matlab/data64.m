% data for Poisson problem with 64 unknowns
% function b
xtemp = [0:0.0157:1];
x = repmat(xtemp, 1, 64);
y = reshape(repmat(xtemp, length(xtemp), 1), 1, length(xtemp)^2);
b = @(x, y) sin(pi*x).*cos(pi*y);
b64 = b(x, y)';
% initial guess x0
x0 = zeros(length(b64), 1);
% tolerance
TOL = 10^-8;
% alowed number of iterations
maxit = 500;
% Poisson matrix
A64 = poissonmatrix(64);

TOL = 10^-8;
maxit = 500;