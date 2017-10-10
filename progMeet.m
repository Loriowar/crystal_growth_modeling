function Y = progMeet(A, B, C, F, kappa, mu)

N = length(A) + 2;
Nhalf = ceil(N / 2);
Alphas = zeros(1, N-1);
Betas = zeros(1, N-1);
Alphas(end) = kappa(2);
Alphas(1) = kappa(1);
Betas(end) = mu(2);
Betas(1) = mu(1);
Y = zeros(1, N);

for i = N-2:-1:Nhalf
    denom = C(i) - Alphas(i+1) * B(i);
    Alphas(i) = A(i) / denom;
    Betas(i) = (B(i) * Betas(i+1) + F(i)) / denom;
end
AlHalf2 = Alphas(Nhalf);
BeHalf2 = Betas(Nhalf);
for i = 1:Nhalf-1
    denom = C(i) - Alphas(i) * A(i);
    Alphas(i+1) = B(i) / denom;
    Betas(i+1) = (A(i) * Betas(i) + F(i)) / denom;
end
AlHalf1 = Alphas(Nhalf);
BeHalf1 = Betas(Nhalf);

Y(Nhalf) = (AlHalf1 * BeHalf2 + BeHalf1) / (1 - AlHalf1 * AlHalf2);

for i = Nhalf:N-1
    Y(i+1) = Alphas(i) * Y(i) + Betas(i);
end

for i = Nhalf-1:-1:1
    Y(i) = Alphas(i) * Y(i+1) + Betas(i);
end

end