function out = aggregate_prices(p_sectoral, rho, nu)
	assert_all(size(p_sectoral)==size(nu));
	assert(rho > 0);
	assert_all(nu > 0);

	out = sum(nu .* p_sectoral .^ (1-rho),2) .^ (1/(1-rho));
end
