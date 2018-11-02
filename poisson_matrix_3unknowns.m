%poisson matrix with 3 unknowns
A=4*eye(3,3)+diag(v,1)+diag(v,-1);
B=-1*eye(3,3);
C= [A, B, zeros(3,3); B, A, B; zeros(3,3), B, A];