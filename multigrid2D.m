% 2D interpolation - bits and pieces

AA = kr_pois(3);

xtemp = [0.1:0.3333:1];
x = repmat(xtemp, 1, 3);
y = reshape(repmat(xtemp, length(xtemp), 1), 1, length(xtemp)^2);
b = @(x, y) sin(pi*x).*cos(pi*y);
b_use = b(x, y)';

x0 = zeros(9, 1);

res = b_use - AA*x0;

% interpolation operation
v(1,:) = res(1)/2;
for i = 1:3
    v(2*i, :) = res(i, :);
    v(2*i+1, :) = (res(i, :)+ res(i+1, :))/2;
end

v(:, 1) = v(:, 1)/2;
for j = 1:3
    v(:, 2*j) = v(:, j);
    v(:, 2*j+1) = (v(:, j)+ v(:, j+1))/2;
end

n = 2^3-1;
k = log2(n+1);

N = 2^(k-1)-1;

RE = zeros(N,n);

for i = 1:N
   RE(i,2*i-1:2*i+1) = [1 2 1]; 
end
RE = RE/4;

II = 2*RE';

% interpolation matrix is 
II2d = kron(II, II);

% restriction matrix is
RE2d = kron(RE, RE);

% coarse grid matrix A is 
A_coarse = RE2d*A*II2d;

% restriction operation
for i = 1:n/2-1
    for j = 1:n/2-1
        v (i, j) = (res(2*i-1, 2*j-1) + res(2*i-1, 2*j+1) + res(2*i+1, 2*j-1) + res(2*i+1, 2*j+1) + 2*(res(2*i, 2*j-1)+res(2*i, 2*j+1) + res(2*i-1, 2*j) + res(2*i+1, 2*j)) + 4*res(2*i, 2*j))\16;
    end
end



