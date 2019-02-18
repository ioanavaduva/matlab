n = size(A,1);

dom_flag = 0;
dom_row = 0;

for i = 1:n
if abs(A(i, i)) >= (sum(abs(A(i, 1:i-1))) + sum(abs(A(i, i+1:n))))
    continue 
else
    dom_flag = 1;
    dom_row = i;
    break
end
end

if dom_flag == 1
    disp("Not strictly dominant")
end
