setup_Poisson_rank1rhs;

opts.tol=1e-4;
emin = 1e-6; 
emax=eigs(M{2},M{1},1,'LA',opts);

s_parameter = logspace(log10(emin)+0.1,log10(emax)-0.1,6);

u_inner = @ (z,s) abs((sqrt((z - emax)./(z - emin))*sqrt((s - emax)./(s - emin)) - 1)./((sqrt((z - emax)./(z - emin))*sqrt((conj(s) - emax)./(conj(s) - emin)) + 1)));

z = linspace(emin, emax, 100);

u_outer = ones(size(z));
for i = 1:6
    u_outer = u_outer.*u_inner(z,s_parameter(i));
end

plot(z, u_outer, 'linewidth', 2), ylim([-1e2,0.7*1e4]);