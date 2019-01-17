data1Dpois; 
x_exact = T1023\b_1d1023;
for maxlev = 2:2:8
    disp("number of levels")
    disp(maxlev)
    [x, iter] = Vmg1d(1023, b_1d1023, T1023, 2/3, 100, 1e-4, maxlev);
    n_it = iter
    n = norm(x_exact - x, 2)/norm(x_exact, 2)
    r_it = b_1d1023 - T1023*x;
    r_b = norm(r_it)/norm(b_1d1023)
end
    
    