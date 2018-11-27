n = 2^3-1;
k = log2(n+1);

N = 2^(k-1)-1;

II = zeros(N,n);

for i = 1:N
   II(i,2*i-1:2*i+1) = [1 2 1]; 
end