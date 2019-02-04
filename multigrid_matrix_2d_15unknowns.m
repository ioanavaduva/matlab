% multigrid v-cycle matrix M -- one level, 3 Jacobi iterations for 15
% unknowns 2D Poisson problem

%original matrix
A = kr_pois(15); 

%restriction (RE) and prolongation (II)
n = 15; 
k = log2(n+1);
N = 2^(k-1)-1;
RE = zeros(N,n);
for i = 1:N
    RE(i,2*i-1:2*i+1) = [1 2 1];
end
RE = RE/4;
II = 2*RE';
P = kron(sparse(II), sparse(II)); % prolongation
T = kron(sparse(RE), sparse(RE)); % restriction
AA = T*A*P;

D = diag(diag(A)); % diagonal matrix from A
Q = -(A - D); % negative remaining matrix of A

H = inv(D)*Q;
R = (H^2)*inv(D) + H*inv(D) + inv(D);

M = (H^3)*R + R + (H^3)*P*AA*T*(H')^3;