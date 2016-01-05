function Z1972 = computeZ1972(d, va, psi, pwt, parameters)

B = parameters.B;
xi = parameters.xi;
theta = parameters.theta;
iBase = parameters.iBase;
iServices = parameters.iServices;
beta = parameters.beta;
kappa = parameters.kappa;
pBase1972 = parameters.pBase1972;
noBaseIndex = parameters.noBaseIndex;

[nCountries, nSectors, ~] = size(va);

t0 = 1;

Z1972 = zeros(nCountries, nSectors);

for n = 1:nCountries
%     n = noBaseIndex(nn);
    for j = 1:nSectors
        Z1972(n, j) = ...
            B(j)^theta * ...
            xi^theta * ...
            d(iBase, n, j, t0) * ...
            pBase1972^(theta * (1 - beta(j))) * ...
            pwt(t0, n)^(theta * (1 - beta(j))) * ...
            kappa(iBase, n, j, t0)^(- theta) * ...
            ((va(n, j, 1) ./ psi(n, j, 1))').^(theta * beta(j));   
    end % for j
end % for n

Z1972(noBaseIndex, iServices) = 0;

end