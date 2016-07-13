function z_services = compute_z_services(va, psi, p_sectoral, parameters)

B = parameters.B;
xi = parameters.xi;
theta = parameters.theta;
beta = parameters.beta; 
i_services = parameters.i_services;
gammas = parameters.gammas;

[n_countries, ~, n_years] = size(va);

z_services = zeros(n_countries, n_years);

for n = 1:n_countries
    for t = 1:n_years
        z_services(n, t) = ...
            xi^theta * ...
            B(i_services, t)^theta * ...
            (va(n, i_services, t) / psi(n, i_services, t))^...
                                            (theta * beta(i_services)) * ...
            prod(p_sectoral(n, :, t).^(gammas(:, i_services, t)'))^theta * ...
            p_sectoral(n, i_services, t)^(-theta);
    end % t
end % n

end