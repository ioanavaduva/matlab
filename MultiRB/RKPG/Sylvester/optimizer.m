n = 3;

h = 1/n; ep = 0.0333;

A = ep*(spdiags([-ones(n, 1) 2*ones(n, 1) -ones(n, 1)],-1:1,n,n))/(h^2);
B = spdiags([-ones(n, 1) zeros(n, 1) ones(n, 1)],-1:1,n,n)/(2*h);

M = A + B;
N = A + B';

b = ones(n, 1);

opts.tol = 1e-2;
emin = min(eig(full(M)));
emax = max(eig(full(M)));



f = @(z) (z-conj(s))/(z+s);

x0 = 1;
Aa = []; % No constraints
ba = [];
Aeq = [];
beq = [];
lb = [];
ub = [];
nonlcon = [];
x1 = fminimax(f,x0,Aa,ba,Aeq,beq,lb,ub,nonlcon);