%  Code to run kpik.m for Poisson with 63 unknowns
cd_ioana3;
A = CD;
nh = n^2;
E = spdiags(rand(nh,1),0,nh,nh);
LE = chol(E,'lower');

m = 100;
tol = 1e-9;
tolY = 1e-12;

[Z,r]=kpik(A,E,LE,b100,m,tol,tolY);

fprintf('final true absolute residual norm: \n')
disp(norm(A*Z*Z'*E+E*Z*Z'*A'+B*B'))    %this matrix should never be formed for n large 