nh=10;
h = 1/nh;
eps = 1;
T = -eps*(diag(2*ones(nh, 1)) + diag (-1*ones(nh-1, 1), 1) + diag (-1*ones(nh-1, 1), -1))/h^2;
I=speye(nh);
A=-(kron(T,I)+kron(I,T));
%n=nh^2;
B=ones(nh, 1);
m=100;
tol=1e-10;
tolY=1e-12;
opts.tol=1e-2;
s1=eigs(-T,1,'largestreal',opts)
s2=eigs(-T,1,'smallestreal',opts)
params.tol=1e-8;
params.smax=s1;
params.smin=s2;
params.m=100;
params.ch=1;
params.period=1;
Z=RKSM_Lyap_real(T,B,params);

fprintf('final true absolute residual norm: \n')
disp(norm(T*Z*Z'+Z*Z'*T'+B*B'))    %this matrix should never be formed for n large 