function [V,T,s]=RKSM_adapt(A,E,EL,EU,b,m,s1,s2,ch) % RKSM from paper

V(:,1)= EU(((A-s1*E)\(EL*b)));
V(:,1)=V(:,1)/norm(V(:,1));
H=zeros(m+2,m+1);
[n]=size(A,1); s=s1; snew=s1; first=1; i=0;
if imag(snew)==0, cmplxflag=0; else cmplxflag=1; end
% additional steps
while i < m
    i=i+1;
    % Select new pole
    if (cmplxflag)   % complex conj pole
        snew=conj(snew);
        cmplxflag=0;
    else
        if first,
            snew=s2; first=0;
        else % New pole
            T=V'*(EL\(A*(EU\V)));  %Note: some redundant ops here
            eH=eig(T);
            if (ch)       % Complex poles. Compute convex hull
                eHpoints = sort([s1; s2.';-eH]);
                if (any(abs(imag(eHpoints)))> 1e-6),
                    ij=convhull(real(eHpoints),imag(eHpoints));
                    eHpoints=eHpoints(ij);
                end
            else   % Real poles s from real set.
                eHpoints=sort([s2,s]);
            end
            snew = newpole(eHpoints,eH,s.');
        end  %first
        % If pole is complex, include its conjugate
        if (imag(snew) ~=0), cmplxflag=1; else cmplxflag=0; end
    end % if (cmplxflag)
    if (ch), s=[s,snew]; else, s=[s,real(snew)]; end
    
    % Arnoldi-type iteration
    wrk = EU*((A-snew*E)\(EL*V(:,i)));
    for it=1:2
        for j=1:i
            t = V(:,j)'*wrk; H(j,i)=H(j,i)+t;
            wrk=wrk-t*V(:,j);
        end;
    end;
    t = norm(wrk);
    H(i+1,i)=t;
    if (t ~= 0.)
        t=1.0/t; V(:,i+1)=wrk*t;
    end;
end;  % while i < m

% Make basis real
mmr=size(V,2);
[V,rr]=qr([real(V),imag(V)],0); V=V(:,1:mmr); T=V'*(EL\(A*(EU\V)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r=ratfun(x,eH,s)

for j=1:length(x)
    r(j)=prod( (x(j)-s)./(x(j)-eH) );
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function snew=newpole(eHpoints,eH,s)

for j=1:length(eHpoints)-1
    sval=linspace(eHpoints(j),eHpoints(j+1),500);
    [sf,jx] = max (abs(ratfun(sval,eH,s)));
    snew(j)=sval(jx);
end
[sn,jx]=max(abs(ratfun(snew,eH,s)));
snew=snew(jx);
return
