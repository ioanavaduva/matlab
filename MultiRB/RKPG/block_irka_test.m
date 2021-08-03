% block irka poles

%Setup
n = 1000;  
h = 1/n; ep = 1;
A = ep*(spdiags([-ones(n, 1) 2*ones(n, 1) -ones(n, 1)],-1:1,n,n))/h^2;
opts.tol = 1e-4;
emin = eigs(A, 1,'smallestabs',opts); 
emax = eigs(A, 1,'largestabs',opts);

%RHS
xtemp = linspace(0,1,n);
x = repmat(xtemp, 1, n);
y = reshape(repmat(xtemp, length(xtemp), 1), 1, length(xtemp)^2);
b = @(x, y) -2.*(6.*x.^2 - 6.*x + 1).*(y-1).^2.*y.^2-2.*(x-1).^2.*x.^2.*(6.*y.^2-6.*y+1);
rhs = b(x, y)'; 
rhs11 = reshape(rhs, n, n); 
col1 = rhs11(:,1); % polynomial column

col2 = ones(n, 1); %ones

% col2 = ((-1).^(0:n-1))'; %alternating 1&-1

rhs1 = [col2, col1];

%IRKA 
m = 16; % number of poles; can change
xi = emin + (emax - emin)*rand(1,m);
B = rhs1*rhs1';

tic;
B_hat = tangential_dir(A, B, xi); %keyboard
[shifts, its] = irka_shifts_block(A, B, B_hat, xi, 1e-2);
shifts = sort(shifts, 'desc');
toc;