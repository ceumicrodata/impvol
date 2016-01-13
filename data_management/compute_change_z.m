function change_z = compute_change_z(d_to_base_change, va_change, psi_change,...
                                  p_base_change, p_sectoral_base_change,...
                                  pwt_change, kappa_to_base_change, parameters)

theta = parameters.theta;
beta = parameters.beta; 
no_base_index = parameters.no_base_index;
i_services = parameters.i_services;

[n_countries, n_sectors, n_years_minus_1] = size(d_to_base_change);
n_years = n_years_minus_1 + 1;

change_z = zeros(n_countries, n_sectors, n_years - 1);

for n = 1:n_countries
    for j = 1:n_sectors
        for t = 1:(n_years - 1)
            change_z(n, j, t) = ...
                d_to_base_change(n, j, t) - ...
                theta * kappa_to_base_change(n, j, t) + ...
                theta * beta(j) * (va_change(n, j, t) - psi_change(n, j, t)) + ...
                theta * (1 - beta(j)) * (p_base_change(t) + pwt_change(t, n)) - ...
                theta * p_sectoral_base_change(t, j);       
        end % for t
    end % for j
end % for n

change_z(no_base_index, i_services, :) = 0;

end