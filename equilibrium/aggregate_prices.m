function out = aggregate_prices(p_sectoral, parameters)
	assert_all(size(p_sectoral), [parameters.n_countries, parameters.n_sectors]);
	assert(parameters.rho > 0);
	assert_all(parameters.nu > 0);

	nu = parameters.nu;
	rho = parameters.rho;

	out = sum(nu .* p_sectoral .^ (1-rho),2)) .^ (1/(1-rho);

	assert_all(size(out), [parameters.n_countries, 1])
end
