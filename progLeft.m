function Y = progLeft(A, B, C, F, kappa, mu)

N = length(A) + 2;
Alphas = zeros(1, N-1);
Betas = zeros(1, N-1);
Alphas(end) = kappa(2);
Betas(end) = mu(2);
Y = zeros(1, N);

for i = N-2:-1:1
    denom = C(i) - Alphas(i+1) * B(i);
    Alphas(i) = A(i) / denom;
    Betas(i) = (B(i) * Betas(i+1) + F(i)) / denom;
end

Y(1) = (mu(1) + kappa(1) * Betas(2)) / (1 - Alphas(2) * kappa(1));

for i = 1:N-1
    Y(i+1) = Alphas(i) * Y(i) + Betas(i);
end

end