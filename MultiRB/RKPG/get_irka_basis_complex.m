function V = get_irka_basis_complex(A, rhs, poles)
% the basis will be given by span{(A+s1I)^{-1}b, ...,(A+skI)^{-1}b}
I = speye(size(A));
V = zeros(length(A(:, 1)), length(poles));

poles(imag(poles)<0)=[];

ind = 0;
v0 = (A + poles(1)*I)\rhs; %compute first column to add to basis
if isreal(poles(1))
    V(:, 1) = v0/norm(v0);
    ind = ind+1;
else
    wr = real(v0); wi = imag(v0);
    wr = wr/norm(wr); %keyboard
    
    wi = wi - wr*(wr'*wi); %keyboard
    wi = wi - wr*(wr'*wi);%keyboard
    wi = wi/norm(wi); %keyboard
    V(:,1:2) = [wr,wi];
    ind = ind+2;
end


for i=2:length(poles)
    w = V(:, ind); %keyboard
    w = (A + poles(i)*I)\w;
    
    if isreal(poles(i))
        w = w - V*(V'*w); %keyboard
        w = w - V*(V'*w);%keyboard
        w = w/norm(w); %keyboard
        V(:, ind+1) = w(:);
        ind = ind+1;
    else
        wr = real(w); wi = imag(w);
        
        wr = wr - V*(V'*wr); %keyboard
        wr = wr - V*(V'*wr);%keyboard
        wr = wr/norm(wr); %keyboard
        V(:,ind+1) = wr;
        
        wi = wi - V*(V'*wi); %keyboard
        wi = wi - V*(V'*wi);%keyboard
        wi = wi/norm(wi); %keyboard
        V(:,ind+2) = wi;
        ind = ind+2;
    end
end
end