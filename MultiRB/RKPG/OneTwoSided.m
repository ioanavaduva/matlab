% Compares results from one-sided and two-sided Rational Krylov 
% projections of Lyapunov problem - requires RKToolbox
%
% J. Pestana 29/6/20

% Set up problem
n = 10; k = 4;
A = full(gallery('tridiag',n,-1,2,-1));
f = ones(n,1);

In = eye(n); 

% Solve full problem
F = f*f'; b = F(:);
B = kron(In,A) + kron(A,In);
x = B\b; X = reshape (x,n,n);
R = A*X + X*A - F;
fprintf('Norm of residual for full problem = %g\n',norm(R));

% Choose pole set
lam_max = max(eig(A));
lam_min = min(eig(A));
% xi = lam_min + (lam_max-lam_min)*rand(1,k);

%%%
% roots denom poles
b = bb;     % sign function on [-10,-1]\cup [1,10]
r = rkfun.gallery('sign', k, b);
% poles(r)
po = imag(poles(r));
poles_Zolo = po(po >= 0);

% Compute the minimum of Z Problem 3 (by transforming c, the min of Prob 4)
K = ellipke(1-1/b^2);
[sn, cn, dn] = ellipj((0:k)*K/k, 1-1/b^2);
extrema = b*dn;
vals = 1-r(extrema);
c = mean( vals(1:2:end) );
e = eig( [ 2-4/c^2 1 ; -1 0 ] );
Zk = min(abs(e));

% Obtain the polynomials p and q of the rational function r = p/q 
[p,q,pq] = poly(r);
pp = [0, p];

% Transform Z Problem 4 into Problem 3 using eq. (2) - rearranged in Theorem 2.1
% (Istace/Thiran paper)and have the denominator given by
denom = q.*(1-Zk) - pp.*(1+Zk);
xi = roots(denom)';
%%% 

% Create rational krylov basis using rktoolbox
V = rat_krylov(A,f,xi);
%V = V(:,k+1);

% Projected matrix and vector
Ak = V'*A*V;
fk = V'*f;

% Identity matrices of correct sizes
Ik = eye(length(V(1,:)));

% Solve one-sided problem
F1 = fk*f'; b1 = F1(:);
B1 = kron(In,Ak) + kron(A,Ik);
y1 = B1\b1; Y1 = reshape(y1,length(V(1,:)),n);
X1 = V*Y1; 
R1 = V'*(A*X1 + X1*A - F);
fprintf('Norm of residual for one-sided problem = %g\n',norm(A*X1 + X1*A - F));

% Solve two-sided problem
F2 = fk*fk'; b2 = F2(:);
B2 = kron(Ik,Ak) + kron(Ak,Ik);
y2 = B2\b2; Y2 = reshape(y2,length(V(1,:)),length(V(1,:)));
X2 = V*Y2*V';
R2 = V'*(A*X2 + X2*A - F)*V;
fprintf('Norm of residual for two-sided problem = %g\n',norm(A*X2 + X2*A - F));
