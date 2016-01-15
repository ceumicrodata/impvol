function zServices = compute_z_services(va, psi, pwt, p_base, p_sectoral, ...
                                          parameters)

B = parameters.B;
xi = parameters.xi;
theta = parameters.theta;
beta = parameters.beta; 
i_services = parameters.i_services;

[n_countries, ~, n_years] = size(va);

zServices = zeros(n_countries, n_years);

for n = 1:n_countries
    for t = 1:n_years
        zServices(n, t) = ...
            xi^theta * ...
            B(i_services)^theta * ...
            (va(n, i_services, t) / psi(n, i_services, t))^...
                                            (theta * beta(i_services)) * ...
            (pwt(t, n) * p_base(t))^(theta * (1 - beta(i_services))) * ...
            p_sectoral(n, i_services, t)^(-theta);
    end % t
end % n

end