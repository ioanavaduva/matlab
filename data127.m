% data for Poisson 2D problem with 127 unknowns
% function b
xtemp = [0:0.0079:1];
x = repmat(xtemp, 1, 127);
y = reshape(repmat(xtemp, length(xtemp), 1), 1, length(xtemp)^2);
b = @(x, y) sin(pi*x).*cos(pi*y);
b127 = b(x, y)';
% initial guess x0
x0 = zeros(length(b127), 1);
% Poisson matrix
A127 = kr_pois(127);
