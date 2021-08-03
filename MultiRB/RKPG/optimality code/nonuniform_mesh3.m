addpath(genpath('./gradedMesh'));

% n = 1000;

% graded 
curve = 2;
weight = 0.1;
gmesh = meshfunc(curve, weight);
[x, h] = gmesh(linspace(0, 1, n+1));

% geometric
% [x,h] = geometric_mesh3(n,0.99);


% f_orig = ones(n,1);

f_orig = ((-1).^(0:n-1))'; 

% xtemp = linspace(0,1,n);
% x = repmat(xtemp, 1, n);
% y = reshape(repmat(xtemp, length(xtemp), 1), 1, length(xtemp)^2);
% b = @(x, y) -2.*(6.*x.^2 - 6.*x + 1).*(y-1).^2.*y.^2-2.*(x-1).^2.*x.^2.*(6.*y.^2-6.*y+1);
% rs = b(x, y)';
% rhs11 = reshape(rs, n, n); f_orig = rhs11(:,1);

% xtemp = linspace(0, 1, n);
% f_orig = sqrt(2)*pi*sin(pi.*xtemp)';

% x = linspace(0, 1, n)';
% f_orig = cos(pi*x);

h_l = h(1:end-1);
h_r = h(2:end);
 
% D = diag(2./(h_l + h_r));
D_inv = sparse(diag(sqrt(2./(h_l + h_r))));
 
T_orig = gallery('tridiag',-1./h(2:n),(h_l+h_r)./(h_l.*h_r),-1./h(2:n));
A = D_inv*T_orig*D_inv; % Scaled version of T
rhs1 = D_inv*f_orig; % Scaled version of f

rhs2 = rhs1;
