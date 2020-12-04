% Approx large IRKA poles by small problem poles

%% Large problem poles

n = 1400; h = 1/n; 
A = (diag(2*ones(n, 1)) + diag(-1*ones(n-1, 1), 1) + diag(-1*ones(n-1, 1), -1))/h^2;

rhs1 = ones(n, 1);

% xtemp = linspace(0,1,n);
% rhs1 = sqrt(2)*pi*sin(pi.*xtemp)';
% rhs2 = rhs1;

% xtemp = linspace(0,1,n);
% x = repmat(xtemp, 1, n);
% y = reshape(repmat(xtemp, length(xtemp), 1), 1, length(xtemp)^2);
% b = @(x, y) -2.*(6.*x.^2 - 6.*x + 1).*(y-1).^2.*y.^2-2.*(x-1).^2.*x.^2.*(6.*y.^2-6.*y+1);
% rhs = b(x, y)';
% rhs11 = reshape(rhs, n, n); rhs1 = rhs11(:,1);

opts.tol=1e-4;
emin = eigs(A, 1, 'smallestabs', opts);
emax = eigs(A, 1,'largestabs',opts);

m = 16;
xi = emin + (emax-emin)*rand(1,m);
tic;
[shifts, its] = irka_shifts(A, rhs1, xi, 1e-4);
toc;

%% Small problem poles -- unscaled

n_small = n/10; h_small = 1/n_small;
A_small = (diag(2*ones(n_small, 1)) + diag(-1*ones(n_small-1, 1), 1) + diag(-1*ones(n_small-1, 1), -1))/h_small^2;
rhs1_small = ones(n_small, 1);

% emin_small = eigs(A_small, 1, 'smallestabs', opts);
% emax_small =  eigs(A_small, 1, 'largestabs', opts);

tic;
[shifts_small, its_small] = irka_shifts(A_small, rhs1_small, xi, 1e-4);
toc;

%% Small problem poles -- scaled
% 
% A_scaled = (diag(2*ones(n_small, 1)) + diag(-1*ones(n_small-1, 1), 1) + diag(-1*ones(n_small-1, 1), -1))/h^2;
% rhs1_small = ones(n_small, 1);
% 
% % emin_scaled = eigs(A_scaled, 1, 'smallestabs', opts);
% % emax_scaled =  eigs(A_scaled, 1, 'largestabs', opts);
% 
% % m_scaled = m-2; % for poles 2
% % m_scaled = m-5; % for poles 1
% % xi_scaled = emin + (emax-emin)*rand(1,m_scaled);
% tic;
% [shifts_scaled, its_scaled] = irka_shifts(A_scaled, rhs1_small, xi, 1e-4);
% toc;

%% Small problems poles combinations

% 1&2 small, 2&3/10 small scaled
% poles = sort([shifts_small(length(shifts_small)/2+1:length(shifts_small)); shifts_small(length(shifts_small)/4+1:3*length(shifts_small)/4)*10]);

% 1&2 small, 3&4*10 small 
poles = sort([shifts_small(length(shifts_small)/2+1:length(shifts_small)); shifts_small(1:length(shifts_small)/2)*10]);

% poles1 = sort([shifts_small(1); shifts_small(length(shifts_small)-3:length(shifts_small)); shifts_scaled], 'descend');
% poles2 = sort([shifts_small(1); shifts_small(length(shifts_small)-3:length(shifts_small)); shifts_scaled(4:length(shifts_scaled))], 'descend');