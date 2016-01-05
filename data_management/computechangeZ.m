function changeZ = computechangeZ(dToBaseChange, vaChange, psiChange,...
                                  pBaseChange, pSectoralBaseChange,...
                                  pwtChange, kappaToBaseChange, parameters)

theta = parameters.theta;
beta = parameters.beta; 
noBaseIndex = parameters.noBaseIndex;
iServices = parameters.iServices;

[nCountries, nSectors, nYears_minus_1] = size(dToBaseChange);
nYears = nYears_minus_1 + 1;

changeZ = zeros(nCountries, nSectors, nYears - 1);

for n = 1:nCountries
    for j = 1:nSectors
        for t = 1:(nYears - 1)
            changeZ(n, j, t) = ...
                dToBaseChange(n, j, t) - ...
                theta * kappaToBaseChange(n, j, t) + ...
                theta * beta(j) * (vaChange(n, j, t) - psiChange(n, j, t)) + ...
                theta * (1 - beta(j)) * (pBaseChange(t) + pwtChange(t, n)) - ...
                theta * pSectoralBaseChange(t, j);       
        end % for t
    end % for j
end % for n

changeZ(noBaseIndex, iServices, :) = 0;

end