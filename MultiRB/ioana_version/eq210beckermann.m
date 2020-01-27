setup_Poisson_rank1rhs;

opts.tol=1e-4;
emin = 1e-6; 
emax=eigs(M{2},M{1},1,'LA',opts);

s_parameter = logspace(log10(emin),log10(emax),6);

u = cell(6,1);

for i = 1:6
   u{i} = @(z) abs((sqrt((z - emax)./(z - emin))*(sqrt(s_parameter(i) - emax)./(s_parameter(i) - emin)) - 1)./((sqrt((z - emax)./(z - emin))*(sqrt(s_parameter(i) - emax)./(s_parameter(i) - emin)) + 1)));
end


z = linspace(emin, emax, 100);
U =  prod(u);
plot(z, U(z), 'linewidth', 2), ylim([-1e4,1e4]);