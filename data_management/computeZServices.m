function ZServices = computeZServices(va, psi, pwt, pBase, pSectoral, ...
                                          parameters)

B = parameters.B;
xi = parameters.xi;
theta = parameters.theta;
beta = parameters.beta; 
iServices = parameters.iServices;

[nCountries, ~, nYears] = size(va);

ZServices = zeros(nCountries, nYears);

for n = 1:nCountries
    for t = 1:nYears
        ZServices(n, t) = ...
            xi^theta * ...
            B(iServices)^theta * ...
            (va(n, iServices, t) / psi(n, iServices, t))^...
                                            (theta * beta(iServices)) * ...
            (pwt(t, n) * pBase(t))^(theta * (1 - beta(iServices))) * ...
            pSectoral(n, iServices, t)^(-theta);
    end % t
end % n

end