function out = computePhi(Z, va, psi, pwt, pBase, d, parameters)

B = parameters.B;
theta = parameters.theta;
beta = parameters.beta; 
% iUs = parameters.iUs;
% kappa = parameters.kappa;

[nCountries, nSectors, nYears] = size(Z);

out = zeros(nCountries, nSectors - 1, nYears);

for n = 1:nCountries
    for j = 1:(nSectors - 1)
        for t = 1:nYears
            out(n, j, t) = ...
                B(j)^(- theta) * ...
                Z(n, j, t) * ...
                (psi(n, j, t) / va(n, j, t))^(theta * beta(j)) * ...
                (pwt(t, n) * pBase(t))^(theta * (beta(j) - 1)) * ...
                (1 / d(n, n, j, t)); 
        end % for t
    end % for j
end % for n

