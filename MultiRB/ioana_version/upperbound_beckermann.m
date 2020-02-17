setup_Poisson_rank1rhs; % set-up of Poisson/convection-diffusion problem

% opts.tol=1e-4;

e = sort(eig(M{2}));
eig_min = e(1); % min eigenvalue
eig_max = e(length(e)); %max eigenvalue

% emin = 1e-6; % min eigenvalue
% emax=eigs(M{2},M{1},1,'LA',opts); % max eigenvalue

s = logspace(log10(eig_min),log10(eig_max), 6); % might need to change to take in the Zolotarev poles 

val =  -fminbnd(@(z) u_out (z, s(1)), -eig_max, -eig_min);

const = 4 + 4*sqrt(2*cond(M{2}));

fval = u_out(val, s(1));

upper_bound = const*fval;