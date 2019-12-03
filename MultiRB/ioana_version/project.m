function [Nm] = project(N, Nm, Pvnew2, V, W, nterm, n, vnew2, iw)

    iwnew = size(vnew2,2);
    for ind = 1:nterm
        wrk2 = N{ind}*Pvnew2;
        Pwrk2 = wrk2;
        newk2 = V(1:n,1:(iw-iwnew))'*Pwrk2; 
        Nm{ind}(1:iw-iwnew,iw-iwnew+1:iw) = newk2;  
        Nm{ind}(iw-iwnew+1:iw,1:iw-iwnew) = newk2';
        Nm{ind}(iw-iwnew+1:iw,iw-iwnew+1:iw) = Pwrk2'*W(:,iw-iwnew+1:iw);  
    end
end

