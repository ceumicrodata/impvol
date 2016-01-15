function P_nj = recover_sectoral_prices(D, P_n)

global beta theta xi B

[N, ~, J] = size(D);
P_nj = zeros(N, J);

for j = 1:J
    P_nj(:, j) = xi * B(j) * ...
                 (D(:, :, j) * P_n.^(theta * (beta(j) - 1))).^(- 1 / theta);
end % j
end