function zeta = computeZeta(d, va, psi, pwt, pBase, parameters)

theta = parameters.theta;
B = parameters.B;
xi = parameters.xi;
beta = parameters.beta;
kappa = parameters.kappa;

[nCountries, nSectors, nYears] = size(va);

zeta = zeros(nCountries, nCountries, nSectors - 1, nYears);

for k = 1:nCountries
    for n = 1:nCountries
        for j = 1:(nSectors - 1)
            for t = 1:nYears
                zeta(k, n, j, t) = ...
                    (B(j) * xi)^theta * ...
                    d(k, n, j, t) * ...
                    kappa(k, n, j, t)^(-theta) * ...
                    (va(n, j, t) / psi(n, j, t))^(theta * beta(j)) * ...
                    (pwt(t, n) * pBase(t))^(theta * (1 - beta(j)));
            end % t
        end % j
    end % n
end % k

zeta = log(zeta);

end