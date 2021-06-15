n = 10;

A = full(gallery('tridiag',n,-1,2,-1));

% Compute eigenvalues of A
theta = ((n:-1:1)*pi/(n+1))';
d1 = 2+2*cos(theta);


% Test
y = randn(n,1);
x = A\y;

x_dst = dst(d1.\(idst(y)));

fprintf('||x-x_dst|| = %g\n',norm(x-x_dst))

%% Note: to get eigenvector matrix explicitly use X = dst(eye(n));

evec = dst(eye(n));
[Q, D] = eig(A);d = diag(D);
