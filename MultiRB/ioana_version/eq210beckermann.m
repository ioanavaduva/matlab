setup_Poisson_rank1rhs;

opts.tol=1e-4;
emin = 1e-6; 
emax=eigs(M{2},M{1},1,'LA',opts);

s_parameter = logspace(log10(emin)+1,log10(emax)-1,6);

% u = cell(6,1);

u_inner = @ (z,s) abs((sqrt((z - emax)./(z - emin))*sqrt((s - emax)./(s - emin)) - 1)./((sqrt((z - emax)./(z - emin))*sqrt((conj(s) - emax)./(conj(s) - emin)) + 1)));

z = linspace(emin, emax, 100);

u_outer = ones(size(z));
for i = 1:6
    u_outer = u_outer.*u_inner(z,s_parameter(i));
end
    
% for i = 1:6
%    u{i} = @(z) abs((sqrt((z - emax)./(z - emin))*(sqrt(s_parameter(i) - emax)./(s_parameter(i) - emin)) - 1)./((sqrt((z - emax)./(z - emin))*(sqrt(s_parameter(i) - emax)./(s_parameter(i) - emin)) + 1)));
% end

% U =  prod(u);
% plot(z, U(z), 'linewidth', 2), ylim([-1e4,1e4]);

plot(z, u_outer, 'linewidth', 2), ylim([-1e4,1e4]);