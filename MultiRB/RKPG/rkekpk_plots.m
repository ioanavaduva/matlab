load('rk_rootsden_mat.mat')
subplot(321);
semilogy(resid_mat(:, 6), 'x');hold on;
[m,i] = min(resid_mat(:,6));
semilogy(i,m, 'o'); hold on; title('Rational Krylov, size 800');
subplot(322);
semilogy(resid_mat(:, 10), 'x');hold on;
[m,i] = min(resid_mat(:,10));
semilogy(i,m, 'o'); hold on; title('Rational Krylov, size 1600')
clear
load('extendedk__mat.mat')
subplot(323);
semilogy(resid_mat(:, 6), 'x');hold on;
[m,i] = min(resid_mat(:,6));
semilogy(i,m, 'o'); hold on; title('Extended Krylov, size 800');
subplot(324);
semilogy(resid_mat(:, 10), 'x');hold on;
[m,i] = min(resid_mat(:,10));
semilogy(i,m, 'o'); hold on; title('Extended Krylov, size 1600');
clear
load('polyk_mat.mat')
subplot(325); semilogy(resid_mat(:, 6), 'x');hold on;
[m,i] = min(resid_mat(:,6));
semilogy(i,m, 'o'); hold on; title('Polynomial Krylov, size 800');
subplot(326);
semilogy(resid_mat(:, 10), 'x');hold on;
[m,i] = min(resid_mat(:,10));
semilogy(i,m, 'o'); hold on; title('Polynomial Krylov, size 1600');

for i = 1:1:6
subplot(3, 2, i)
xlabel('Iterations');
ylabel('Residual');
hold on
end