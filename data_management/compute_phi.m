function out = compute_phi(z, va, psi, pwt, p_base, d, parameters)

B = parameters.B;
theta = parameters.theta;
beta = parameters.beta; 
% iUs = parameters.iUs;
% kappa = parameters.kappa;

[n_countries, n_sectors, n_years] = size(z);

out = zeros(n_countries, n_sectors - 1, n_years);

for n = 1:n_countries
    for j = 1:(n_sectors - 1)
        for t = 1:n_years
            out(n, j, t) = ...
                B(j)^(- theta) * ...
                z(n, j, t) * ...
                (psi(n, j, t) / va(n, j, t))^(theta * beta(j)) * ...
                (pwt(t, n) * p_base(t))^(theta * (beta(j) - 1)) * ...
                (1 / d(n, n, j, t)); 
        end % for t
    end % for j
end % for n

