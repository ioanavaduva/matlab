function V = get_start_basis(rhs)
    V = rhs(:,1)/norm(rhs(:, 1));
    for i = 2:length(rhs(1,:))
        v = rhs(:,i) - V*(V'*rhs(:,i));
        v = rhs(:,i) - V*(V'*rhs(:,i));
        v = v/norm(v);
        V = [V, v(:)];
    end
end

