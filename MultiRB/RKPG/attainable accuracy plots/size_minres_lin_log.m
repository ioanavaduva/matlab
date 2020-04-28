n = [50; 100; 200; 400; 600; 800; 1000; 1200; 1600; 2000; 2400];
rk_res = [1.607115e-13; 6.3571e-13; 2.2636e-12; 9.7273e-12; 4.5327e-11; 7.8354e-11; 1.0150e-10; 1.2137e-10; 3.5882e-10; 4.0422e-10; 8.4849e-10];
ek_res = [1.9556e-13; 9.8885e-13; 4.1082e-12; 1.8025e-11; 6.0872e-11; 1.0500e-10; 1.5122e-10; 1.9931e-10; 4.7599e-10; 6.3803e-10; 1.1256e-09];
nn = [50; 100; 200; 400; 600; 800; 1000; 1200; 1600];
pk_res = [1.1856e-13; 1.5054e-12; 8.7500e-12; 6.4914e-11; 1.6859e-10; 3.3783e-10; 5.5451e-10; 8.6700e-10; 1.7052e-09];

%subplot(121); %linear 
plot(n, rk_res, 'x'); hold on;
p = polyfit(n, rk_res, 3);
f = polyval(p, n);
plot(n, f); hold on;

plot(n, ek_res, 'o'); hold on;
p2 = polyfit(n, ek_res, 3);
f2 = polyval(p2, n);
plot(n, f2); hold on;

plot(nn, pk_res', 'v'); hold on;
p3 = polyfit(nn, pk_res, 3);
f3 = polyval(p3, nn);
plot(nn, f3); hold on;


% subplot(122); %logarithmic
% semilogy(n, rk_res, 'x'); hold on;
% semilogy(n, ek_res, 'o'); hold on;
% semilogy(nn, pk_res, '+'); hold on;