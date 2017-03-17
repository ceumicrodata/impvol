function out = compute_p_services(pwt, p_base, parameters)

n_countries = parameters.n_countries;
n_sectors = parameters.n_sectors;
n_years = parameters.n_years;
final_expenditure_share = parameters.final_expenditure_share;
i_services = parameters.i_services;
i_base = parameters.i_base;
rho = parameters.rho;

out = zeros(n_countries,n_years);
for t = 1:n_years
    for n = 1:n_countries
        out(n,t) = pwt(t,n) .* p_base(t,1) .*(final_expenditure_share(n,i_services,t)) .^ (1/(1 - rho));
    end
end