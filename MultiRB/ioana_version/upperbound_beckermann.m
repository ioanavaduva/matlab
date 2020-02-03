setup_Poisson_rank1rhs; % set-up of Poisson/convection-diffusion problem

opts.tol=1e-4;
emin = 1e-6; % min eigenvalue
emax=eigs(M{2},M{1},1,'LA',opts); % max eigenvalue

s = logspace(log10(emin)+0.1,log10(emax)-0.1, 6); % might need to change to take in the Zolotarev poles 

gamma = - fminbnd(@u_out, -emax, -emin);

upper_bound = (4 + 4*sqrt(2*cond(M{2})))*gamma;

