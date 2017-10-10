function Y = prog(A, B, C, F, kappa, mu)

N = length(A) + 2;
Alphas = zeros(1, N-1);
Betas = zeros(1, N-1);
Alphas(1) = kappa(1);
Betas(1) = mu(1);
Y = zeros(1, N);

for i = 1:N-2
    denom = C(i) - Alphas(i) * A(i);
    Alphas(i+1) = B(i) / denom;
    Betas(i+1) = (A(i) * Betas(i) + F(i)) / denom;
end

Y(N) = (mu(2) + kappa(2) * Betas(N-1)) / (1 - Alphas(N-1) * kappa(2));

for i = N-1:-1:1
    Y(i) = Alphas(i) * Y(i+1) + Betas(i);
end

end