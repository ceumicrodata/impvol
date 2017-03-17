function nu = compute_nus(final_expenditure_share, p_sectoral, parameters, c)
    n_countries = parameters.n_countries;
    n_sectors = parameters.n_sectors;
    n_years = parameters.n_years;
    numerical_zero = parameters.numerical_zero;
    rho = c.rh;
    
    assert_all(size(final_expenditure_share) == [n_countries, n_sectors, n_years]);
    assert_all(final_expenditure_share(:) > 0);
    assert_all(abs(sum(final_expenditure_share, 2) - 1) < numerical_zero);
    assert_all(size(p_sectoral) == [n_countries, n_sectors, n_years]);
    assert_all(p_sectoral(:) > 0);
    
    nu = final_exp_share .* p_sectoral.^(rho - 1);
    
    assert(size(nu) == [n_countries, n_sectors, n_years]);
end

