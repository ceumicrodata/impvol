function out = get_expenditure_shares(p_sectoral, rho, nu_input)
	assert_all(size(p_sectoral)==size(nu_input));
	assert_all(p_sectoral>0);
	assert_all(nu_input>0);
	assert(rho>0);
	[N, J] = size(nu_input);

	%out = nu_input .* p_sectoral .^ (1-rho);
	%out = out ./ repmat(sum(out,2), [1, J]);

	%% DEBUG: equal shares
	out = repmat(1/J, N, J);

	assert_all(size(out)==[N, J]);
	assert_all(out>0);
	assert_all(out<1);
end 