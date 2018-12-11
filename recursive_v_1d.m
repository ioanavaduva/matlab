
data1Dpois;
x1 = zeros(1023, 1);

for i = 1:2
    res = b_1d1023 - T1023*x1;
    if norm(res) > 1e-4
        x1 = vcycle1d(1023, x1, b_1d1023, T1023, 2/3, 4);
    else
        break;
    end
end
