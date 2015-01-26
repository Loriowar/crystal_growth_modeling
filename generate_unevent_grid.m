function [ interval, h ] = generate_unevent_grid( N, q )
b1 = @(q, n) (q-1)./(q.^n-1);
h_fun = @(b1, q, N) b1*(q.^(N-1));
h = b1(q, N);
for i = 2:N
    h(i) = h_fun(b1(q, N), q, i);
end
interval = cumsum(h);
end

