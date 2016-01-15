function z1972 = compute_z_1972(d, va, psi, pwt, parameters)

B = parameters.B;
xi = parameters.xi;
theta = parameters.theta;
i_base = parameters.i_base;
i_services = parameters.i_services;
beta = parameters.beta;
kappa = parameters.kappa;
p_base_1972 = parameters.p_base_1972;
no_base_index = parameters.no_base_index;

[n_countries, n_sectors, ~] = size(va);

t0 = 1;

z1972 = zeros(n_countries, n_sectors);

for n = 1:n_countries
%     n = no_base_index(nn);
    for j = 1:n_sectors
        z1972(n, j) = ...
            B(j)^theta * ...
            xi^theta * ...
            d(i_base, n, j, t0) * ...
            p_base_1972^(theta * (1 - beta(j))) * ...
            pwt(t0, n)^(theta * (1 - beta(j))) * ...
            kappa(i_base, n, j, t0)^(- theta) * ...
            ((va(n, j, 1) ./ psi(n, j, 1))').^(theta * beta(j));   
    end % for j
end % for n

z1972(no_base_index, i_services) = 0;

end