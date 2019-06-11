function [y] = kronprod(A, B, z)
    Z = vec2mat(z, length(A));
    Y = B*Z*A' + A*Z*B';
    y = Y(:);
end

