% Compares results from one-sided and two-sided Rational Krylov 
% projections of Lyapunov problem - requires RKToolbox
%
% J. Pestana 29/6/20

% Set up problem
n = 10; k = 4;
A = full(gallery('tridiag',n,-1,2,-1));
f = ones(n,1);

In = eye(n); Ik = eye(k);

% Solve full problem
F = f*f'; b = F(:);
B = kron(In,A) + kron(A,In);
x = B\b; X = reshape (x,n,n);
R = A*X + X*A - F;
fprintf('Norm of residual for full problem = %g\n',norm(R));

% Choose pole set
lam_max = max(eig(A));
lam_min = min(eig(A));
xi = lam_min + (lam_max-lam_min)*rand(1,k);

% Create rational krylov basis using rktoolbox
W = rat_krylov(A,f,xi);
W = W(:,1:k);

% Projected matrix and vector
Ak = V'*A*V;
fk = V'*f;

% Identity matrices of correct sizes
In = eye(n); Ik = eye(k);

% Solve one-sided problem
F1 = fk*f'; b1 = F1(:);
B1 = kron(In,Ak) + kron(A,Ik);
y1 = B1\b1; Y1 = reshape(y1,k,n);
X1 = V*Y1; %norm(A*X1 + X1*A - F)
R1 = V'*(A*X1 + X1*A - F);
fprintf('Norm of residual for one-sided problem = %g\n',norm(R1));

% Solve two-sided problem
F2 = fk*fk'; b2 = F2(:);
B2 = kron(Ik,Ak) + kron(Ak,Ik);
y2 = B2\b2; Y2 = reshape(y2,k,k);
X2 = V*Y2*V'; %norm(A*X2 + X2*A - F)
R2 = V'*(A*X2 + X2*A - F)*V;
fprintf('Norm of residual for two-sided problem = %g\n',norm(R2));
