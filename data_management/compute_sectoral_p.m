function [p_sectoral, p_base, p_sectoral_base] = compute_sectoral_p(p_sectoral_data, pwt, d, parameters)
    n_countries = parameters.n_countries;
    n_years = parameters.n_years;
    i_services = parameters.i_services;
    final_expenditure_share = parameters.final_expenditure_share;
    has_prices = parameters.has_prices;
    i_base = parameters.i_base;
    rho = parameters.rho;
    
    % index of base country in sectoral price data
    ii = sum(has_prices(1:i_base));
    p_sectoral_base = p_sectoral_data((ii - 1) * n_years + 1 : ii * n_years, :);

    %% Process imported matrices
    % Normalize sectoral price index in base country 
    p_sectoral_base = p_sectoral_base ./ repmat(p_sectoral_base(1, :), [n_years 1]);
    
    % Compute aggregate price index in base country
    p_base = p_sectoral_base(:, i_services) .* squeeze(final_expenditure_share(i_base, i_services, :)).^(1/(rho - 1));
    
    %% Calculate sectoral prices
    kappa = parameters.kappa;
    theta = parameters.theta;

    p_sectoral = ...
        exp(squeeze(mean(1/theta * log(bsxfun(@rdivide, d, d(i_base, :, :, :))) - ...
                    log(bsxfun(@rdivide, kappa, kappa(i_base, :, :, :))), 2)) + ...
            permute(repmat(log(p_sectoral_base), [1, 1, n_countries]), [3, 2, 1]));

    p_sectoral(:, i_services, :) = ...
        compute_p_services(pwt, p_base, parameters);
end

