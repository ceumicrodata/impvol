function kappa = computetradecost(d, parameters) 

theta = parameters.theta;
numericalZero = parameters.numericalZero;
iServices = parameters.iServices;

% Compute trade costs from expenditure shares

[nCountries, ~, nSectors, nYears] = size(d);


% Trade costs are calculated based on the symmetric formula of (26)
kappa = zeros(nCountries, nCountries, nSectors, nYears);
for j = 1:(nSectors - 1)
    for t = 1:nYears
        kappa(:, :, j, t) = ...
            ((d(:, :, j, t) .* d(:, :, j, t)') ./ ...
             (diag(d(:, :, j, t)) * diag(d(:, :, j, t))')).^(1 / (2 * theta));
         
         % normalize
%          kappa(:, :, j, t) = ...
%              kappa(:, :, j, t) ./ ...
%              repmat(max(kappa(:, :, j, t), [], 2), [1, nCountries]);
%          
%          for n = 1:nCountries
%              kappa(n, n, j, t) = 1;
%          end % n
    end % t
end % j

kappa(:, :, iServices, :) = repmat(eye(nCountries), [1, 1, 1, nYears]);

kappa(:, :, 1:(nSectors - 1), :) = ...
    max(kappa(:, :, 1:(nSectors - 1), :), numericalZero);
kappa = min(kappa, 1);

end