% rktoolbox Zolotarev third problem vs logspace poles
clear all

addpath(genpath('../../rktoolbox'));
addpath(genpath('../../MultiRB'));
% setup for logspace shifts
setup_Poisson_rank1rhs;

opts.tol=1e-4;
emin2 = 1e-6; 
emax2=eigs(M{2},M{1},1,'LA',opts);

aa=emin2;
bb=emax2;



%setup for Zolotarev shifts
k = 12;      % rational degree
b = emax2-emin2+1;     % sign function on [-10,-1]\cup [1,10]
r = rkfun.gallery('sign', k/2, b);

s_parameter = logspace(log10(1),log10(b),6);

% Extrema for [-1,-1/b]\cup [1/b,1]:
% K = ellipke(1-1/b^2);
% [sn, cn, dn] = ellipj((0:k)*K/k, 1-1/b^2);
% extrema = b*dn;   % Transplant to [-b,-1]\cup [1,b]
% 
% vals = 1-r(extrema);
% c = mean( vals(1:2:end) );
% e = eig( [ 2-4/c^2 1 ; 1 0 ] );
% Zk = min(abs(e))
% 
% % Mobius transformation of r(x):
% R = @(x) (1 + (1+Zk)/(1-Zk)*r(x))./(1 - (1+Zk)/(1-Zk)*r(x));
% x = linspace(-b, b, 5000);
% plot(x, R(x), 'linewidth', 2), ylim([-1e4,1e4])
% xlabel('x')
% title('solution to Zolotarev''s third problem'), hold on
% 
% semilogy(s_parameter, '-o')

