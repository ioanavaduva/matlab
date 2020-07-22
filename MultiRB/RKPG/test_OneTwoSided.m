% Set up problem
n = 10; k = 4;
A = full(gallery('tridiag',n,-1,2,-1));
f = ones(n,1);

In = eye(n); Ik = eye(k);

% Choose pole set
lam_max = max(eig(A));
lam_min = min(eig(A));
% xi = lam_min + (lam_max-lam_min)*rand(1,k);
xi=[1.37449523950303,3.67993398690207,3.43280769287445,2.15605542644147]';
% basis 
V = f/norm(f);
for i = 1:4
    V = get_rk_basis(A, xi(i), V); 
end
% V = V(:,1:4);

W = rat_krylov(A, f, xi); 
W = W(:, 1:4); % V & W are clearly not the same

%solve original problem
F = f*f'; b = F(:);
B = kron(In,A) + kron(A,In); 
x = B\b; X = reshape (x,n,n);
R = A*X + X*A - F;

% Solve with each of the bases: 2-sided my basis
Ap1 = V'*A*V;
fp1 = V'*f;
F1 = fp1*fp1';

%solve with lyap
Y1_lyap = lyap(-Ap1, F1); 
%solve directly
b1 = F1(:);
B1 = kron(Ik,Ap1) + kron(Ap1,Ik);
y1_dir = B1\b1; Y1_dir = reshape(y1_dir,k,k); %both direct and lyap solvers give the same inner solution for two sided
Xk1 = V*Y1_lyap*V';

Ak = V'*A*V;
fk = V'*f;
Fk = fk*f';

Yk_lyap = lyap(-Ak, -A, Fk);
bk = Fk(:);
Bk = kron(In, Ak)+kron(A, Ik);
yk_dir = Bk\bk; 
Yk_dir = reshape(yk_dir, k, n);
Xk = V*Yk_lyap;

% %2-sided Stefan basis
% Ap2 = W'*A*W;
% fp2= W'*f;
% F2 = fp2*fp2';
% 
% %solve with lyap
% Y2_lyap = lyap(-Ap2, F2); 
% %solve directly
% b2 = F2(:);
% B2 = kron(Ik,Ap2) + kron(Ap2,Ik);
% y2_dir = B2\b2; Y2_dir = reshape(y2_dir,k,k);
