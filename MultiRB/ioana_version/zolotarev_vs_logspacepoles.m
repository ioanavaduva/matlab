% rktoolbox Zolotarev third problem vs logspace poles
clear all

addpath(genpath('../../rktoolbox'));
addpath(genpath('../../MultiRB'));

% setup for logspace shifts
setup_Poisson_rank1rhs;

e = sort(eig(M{2}));
emin2 = e(1); % min eigenvalue
emax2 = e(length(e)); %max eigenvalue

bb = emax2 - emin2 + 1;

s_parameter = logspace(log10(1), log10(bb), 6);

%setup for Zolotarev shifts
k = 12;      % rational degree
b = bb;     % sign function on [-10,-1]\cup [1,10]
r = rkfun.gallery('sign', k/2, b);

po = imag(poles(r));

poles_positive = po(po>=0 );

plot(poles_positive, 'o'); hold on;
plot(s_parameter, 'x'); hold off
set(gca, 'YScale', 'log');

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

