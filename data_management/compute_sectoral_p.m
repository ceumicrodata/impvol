function [p_sectoral, p_base, p_sectoral_base] = compute_sectoral_p(p_sectoral_data, pwt, d, c, parameters)
    
    n_countries = parameters.n_countries;
    n_years = parameters.n_years;
    i_base = parameters.i_base;
    i_services = parameters.i_services;
    final_exp_share = parameters.alpha_t;
    alpha_replace = parameters.alpha_r;
    
    % index of base country in sectoral price data
    ii = sum(c.has_prices(1:c.i_base));
    p_sectoral_base = p_sectoral_data((ii - 1) * n_years + 1 : ii * n_years, :);

    % aa = cumsum(has_prices);
    % 
    % p_sectoral_data2 = zeros(n_countries, n_sectors, n_years);
    % for n = 1:n_countries
    %     if has_prices(n) == 1
    %         p_sectoral_data2(n, :, :) =  p_sectoral_data((aa(n) - 1) * n_years + 1 : aa(n) * n_years, : )';
    % %         p_sectoral_data2(n, :, :) = bsxfun(@rdivide, p_sectoral_data2(n, :, :), p_sectoral_data2(n, :, 1));  
    %     end % if
    % end % for n


    % p_sectoral_data2 = bsxfun(@rdivide, p_sectoral_data2, p_sectoral_data2(25, :, 1));

    %% Process imported matrices
    % Normalize sectoral price index in base country 
    p_sectoral_base = p_sectoral_base ./ repmat(p_sectoral_base(1, :), [n_years 1]);
    
    % Compute aggregate price index in base country
    p_base = prod((p_sectoral_base' ./ alpha_replace) .^ alpha_replace, 1)';
    
%     % Compute aggregate price index in base country - is this OK?
%     p_base = prod((p_sectoral_base' ./ squeeze(final_exp_share(i_base,:,:))) .^ squeeze(final_exp_share(i_base,:,:)), 1)';
    
    %% Calculate sectoral prices
    kappa = parameters.kappa;
    theta = parameters.theta;

    p_sectoral = ...
        exp(squeeze(mean(1/theta * log(bsxfun(@rdivide, d, d(i_base, :, :, :))) - ...
                    log(bsxfun(@rdivide, kappa, kappa(i_base, :, :, :))), 2)) + ...
            permute(repmat(log(p_sectoral_base), [1, 1, n_countries]), [3, 2, 1]));

    p_sectoral(:, i_services, :) = ...
        compute_p_services(pwt, p_sectoral, p_base, parameters);
%     p_sectoral(:, i_services, :) = compute_p_services(pwt, parameters);

    % why have we calculated service prices if we fix them to unity?
    % p_sectoral(:,i_services, :) = 1;
    
    
    % nu(services) = 1 for all n,t,
    % other sectoral prices normalized to 1 initially

    % 
    % %% check sectoral prices
    % close all
    % 
    % for j = 1:24
    %     figure()
    % 
    %     n = 2;
    %     n1 = sum(has_prices(1:n));
    %    
    % 
    %     a = p_sectoral_data(1 + (n1 - 1) * 36 : n1 * 36, j);
    %     a = a / a(1);
    % 
    %     aa = squeeze(p_sectoral(n, j, :));
    %     aa = aa / aa(1);
    % 
    %     plot([a, aa])
    % 
    % end
end

