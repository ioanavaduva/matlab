addpath(genpath('./gradedMesh'));


%% Mesh function
n = 1001;
curve = 2; 
weight = 0.1;
gmesh = meshfunc(curve, weight);
[x, h] = gmesh(linspace(0, 1, n+1)); 

%% geometric mesh
% n = 500;
% r = 0.9;
% [x, h] = geometric_mesh2(n, r);

%% RHS

% f_orig = ones(n-1,1);

% f_orig = ((-1).^(0:n-2))'; 

% xtemp = linspace(0,1,n-1);
% x = repmat(xtemp, 1, n-1);
% y = reshape(repmat(xtemp, length(xtemp), 1), 1, length(xtemp)^2);
% b = @(x, y) -2.*(6.*x.^2 - 6.*x + 1).*(y-1).^2.*y.^2-2.*(x-1).^2.*x.^2.*(6.*y.^2-6.*y+1);
% rs = b(x, y)';
% rhs11 = reshape(rs, n-1, n-1); f_orig = rhs11(:,1);

% xtemp = linspace(0, 1, n-1);
% f_orig = sqrt(2)*pi*sin(pi.*xtemp)';

% x = linspace(0, 1, n-1)';
% f_orig = cos(pi*x);

%% Obtaining the matrices and rhs 

h_l = h(1:end-2);
h_r = h(3:end);

D = diag(2./(h_l + h_r));
D_inv = diag(sqrt(2./(h_l + h_r)));

T_orig = gallery('tridiag',-1./h(2:n-1),(h_l+h_r)./(h_l.*h_r),-1./h(2:n-1));
A = D_inv*T_orig*D_inv; % Scaled version of T
% rhs1 = D_inv*f_orig; % Scaled version of f
% 
% rhs2 = rhs1;

