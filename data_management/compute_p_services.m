function out = compute_p_services(pwt, p_sectoral, p_base, parameters)

alpha = parameters.alpha;
i_services = parameters.i_services;

[n_countries, n_sectors, n_years] = size(p_sectoral);
   

out = zeros(n_countries, n_years);

for t = 1:n_years
    for n = 1:n_countries
        out(n, t) = ...
            (pwt(t, n) * p_base(t))^(1 / alpha(i_services, t)) * ...
            prod(alpha(:, t).^(- alpha(:, t)))^(- 1 / alpha(i_services, t)) * ...
            prod(p_sectoral(n, 1:(n_sectors - 1), t) .^ ...
                 (alpha(1:(n_sectors - 1), t)'))^(- 1 / alpha(i_services, t));
    end % n 
end % t 

end