% To obtain the Ritz values plots run the following (need to have e_Ap from
% RKPG outputs

 for i=1:it
semilogy(i*ones(length(e_Ap{i})), e_Ap{i}, 'o'); hold on;
end
ee = eig(A);
plot(zeros(length(A)), ee, 'v')