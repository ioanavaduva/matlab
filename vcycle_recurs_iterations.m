function [x,it] = vcycle_recurs_iterations(n, b, T, w, maxit, TOL,levels)

x = zeros(n,1);
if norm(b) < TOL, it = 0; return; end

for it = 1:maxit
    [x] = vcycle1d_original(n, b, T, w, TOL,levels,x);
    if norm(b-T*x) < TOL, return; end
    
end

end