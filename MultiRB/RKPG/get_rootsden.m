function [roots_denom, extrema] = get_rootsden(k, b)
    r = rkfun.gallery('sign', k, b);
    % poles(r)
    po = imag(poles(r));
    poles_Zolo = po(po >= 0);

    % % Compute the minimum of Z Problem 3 (by transforming c, the min of Prob 4)
    K = ellipke(1-1/b^2); %keyboard
    [sn, cn, dn] = ellipj((0:k)*K/k, 1-1/b^2); %keyboard
    extrema = b*dn; %keyboard
    vals = 1-r(extrema); %keyboard
    c = mean( vals(1:2:end) );
    e = eig( [ 2-4/c^2 1 ; -1 0 ] );
    Zk = min(abs(e));

    % Obtain the polynomials p and q of the rational function r = p/q 
    [p,q,pq] = poly(r);
    pp = [0, p];

    % Transform Z Problem 4 into Problem 3 using eq. (2) - rearranged in Theorem 2.1
    % (Istace/Thiran paper)and have the denominator given by
    denom = q.*(1-Zk) - pp.*(1+Zk);
    roots_denom = roots(denom);
end