function zeta = compute_zeta(d, va, psi, pwt, p_sectoral, parameters)

theta = parameters.theta;
B = parameters.B;
xi = parameters.xi;
beta = parameters.beta;
kappa = parameters.kappa;
gammas = parameters.gammas;

[n_countries, n_sectors, n_years] = size(va);

zeta = zeros(n_countries, n_countries, n_sectors - 1, n_years);

for k = 1:n_countries
    for n = 1:n_countries
        for j = 1:(n_sectors - 1)
            for t = 1:n_years
                zeta(k, n, j, t) = ...
                    (B(j, t) * xi)^theta * ...
                    d(k, n, j, t) * ...
                    kappa(k, n, j, t)^(-theta) * ...
                    (va(n, j, t) / psi(n, j, t))^(theta * beta(j)) * ...
                    prod(p_sectoral(n, :, t).^(gammas(:, j)'))^theta;
            end % t
        end % j
    end % n
end % k

zeta = log(zeta);

end