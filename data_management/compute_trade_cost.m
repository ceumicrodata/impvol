function kappa = compute_trade_cost(d, parameters) 

theta = parameters.theta;
numerical_zero = parameters.numerical_zero;
i_services = parameters.i_services;

% Compute trade costs from expenditure shares

[n_countries, ~, n_sectors, n_years] = size(d);


% Trade costs are calculated based on the symmetric formula of (26)
kappa = zeros(n_countries, n_countries, n_sectors, n_years);
for j = 1:(n_sectors - 1)
    for t = 1:n_years
        kappa(:, :, j, t) = ...
            ((d(:, :, j, t) .* d(:, :, j, t)') ./ ...
             (diag(d(:, :, j, t)) * diag(d(:, :, j, t))')).^(1 / (2 * theta));
         
         % normalize
%          kappa(:, :, j, t) = ...
%              kappa(:, :, j, t) ./ ...
%              repmat(max(kappa(:, :, j, t), [], 2), [1, n_countries]);
%          
%          for n = 1:n_countries
%              kappa(n, n, j, t) = 1;
%          end % n
    end % t
end % j

kappa(:, :, i_services, :) = repmat(eye(n_countries), [1, 1, 1, n_years]);

kappa(:, :, 1:(n_sectors - 1), :) = ...
    max(kappa(:, :, 1:(n_sectors - 1), :), numerical_zero);
kappa = min(kappa, 1);

end