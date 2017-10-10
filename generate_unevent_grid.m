function [ interval, h ] = generate_unevent_grid( N, q )
b1 = @(q, n) (q-1)./(q.^n-1);
h_fun = @(b1, q, N) b1*(q.^(N-1));
% interval = 0;
% interval(2) = b1(q, N);
h = b1(q, N);
% h(2) =  h_fun(b1(q, N), q, 2);
for i = 2:N
    h(i) = h_fun(b1(q, N), q, i);
%     interval(i) = interval(i-1) + h(i);   
end
interval = cumsum(h);
end

