% standard Krylov subspace

function K = ks(A, b, n) 
    k{1} = b;
    K = [k{1}];
    for i = 1:n-1
        k{i+1} = A^i * b;
        K = [K, k{i+1}];
    end
end

