% n (number of unknowns) must be of the form n = 2^p-1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% data for Poisson 1D - 7 unknowns

% create vector b
    x = [0.1:0.1429:1];
    b = @(x) sin(pi*x).*cos(pi*x);
    b_1d7 = b(x)';
    
% create matrix T of coefficients (fine grid matrix)
    T7 = oned_pois(7);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% data for Poisson 1D - 15 unknowns

% create vector b
    x = [0.1:0.0625:1];
    b = @(x) sin(pi*x).*cos(pi*x);
    b_1d15 = b(x)';
    
% create matrix T of coefficients (fine grid matrix)
    T15 = oned_pois(15);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% data for Poisson 1D - 31 unknowns

% create vector b
    x = [0.1:0.0294:1];
    b = @(x) sin(pi*x).*cos(pi*x);
    b_1d31 = b(x)';
    
% create matrix T of coefficients (fine grid matrix)
    T31 = oned_pois(31);
    