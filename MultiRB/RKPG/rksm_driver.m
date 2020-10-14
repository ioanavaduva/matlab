% rksm.m driver

% Input:
%
% A, E  coeff. matrices. A<0,  E is spd
% EL   Cholesky lower factor of E
% B     rhs factor
% m       max space dimension allowed
% tol     stopping tolerance, with stopping criterion
%          ||LE\A X LE  + LE' X A'/LE'-LE\BB'/LE'||
%          ----------------------------------------  < tol
%      ||LE\BB'/LE'|| + ||E^{-1}|| ||A|| ||LE'X LE ||
%         computed in a cheap manner
% s1,smax estimates for real spectral interval
%         associated with field of values of (A,E)
% ch      ch=1  complex poles  ch=0 real poles
% tolY    truncation tolerance for final solution, e.g., tolY=1e-12
clear all;

n = 100;
h = 1/n; eps = 1;
A = -(eps*(diag(2*ones(n, 1)) + diag(-1*ones(n-1, 1), 1) + diag(-1*ones(n-1, 1), -1))/h^2);
I = speye(n);
EL = chol(I);
B = ((-1).^(0:n-1))';%ones(n, 1);
m = 300; tol = 1e-9;
emin = 1e-6;  
opts.tol=1e-4;
emax = eigs(A, 1,'LA',opts);
ch = 0;
tolY = 1e-12;

[Z, nrmrestot, poles] = rksm(A, I, I, B, m, tol, emin, emax, 0, 1e-12);
% [V, T, s] = RKSM_adapt(A, E, E, E, B, m, emin, emax, ch);
% 
% fprintf('final true absolute residual norm: \n')
% disp(norm(A*Z*Z' + Z*Z'*A' + B*B'))  
% 
% % compare with 'backslash' solution
% AA = -(kron(A, I) + kron(I, A));
% rhs = kron(B, B);
% X = AA\rhs;
% X_resh = reshape(X, n, n);
% X_approx = Z*Z';
% fprintf('norm difference of exact and approx solutions: \n')
% disp(norm(X_resh - X_approx))
% 
% Re = real(eig(-A));
% Im = imag(eig(-A));
% 
% semilogx(Re, Im, 'o'); hold on;
% 
% re = real(poles);
% im = imag(poles);
% plot(re, im, 's'); hold on;