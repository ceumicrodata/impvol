function final_expenditure_share = compute_final_expenditure_share(D, va, gammas, beta, parameters, c)
    n_countries = parameters.n_countries;
    n_sectors = parameters.n_sectors;
    n_years = parameters.n_years;
    numerical_zero = parameters.numerical_zero;
    weights = c.filter_weights;
    
    assert_all(size(D) == [n_countries*n_sectors, n_countries*n_sectors, n_years]);
    assert_all(size(va) == [n_countries*n_sectors, n_years]);
    assert_all(size(gammas) == [n_sectors, n_sectors, n_years]);
    assert_all(size(beta) == [n_sectors, n_years]);
    assert_all(abs(squeeze(sum(gammas, 1)) + beta - 1) < numerical_zero);
    
    final_exp = zeros(n_countries*n_sectors, n_years);    
    % Compute final expenditure shares
    for t = 1:n_years
        final_exp(:,t) = ( (inv(D(:,:,t)) - kron(eye(n_countries), squeeze(gammas(:,:,t)))) / kron(eye(n_countries), diag(beta(:,t))) ) * va(:,t);
        
        %final_expenditure_share(:,:,t) = reshape(final_exp(:,t), [n_sectors,n_countries])' ./ repmat(sum(reshape(final_exp(:,t), [n_sectors,n_countries])', 2), [1, n_sectors]);
    end
    
    % Detrend the final expenditure
    [final_exp_tr,~] = detrend_series(final_exp, weights);

    % Replace the final expenditure with 0 whenever below 0
    final_exp_tr(final_exp_tr < numerical_zero) = numerical_zero;

    % Calculate shares
    final_exp_tr_matr = wide(final_exp_tr, parameters);

    final_expenditure_share = ones(n_countries, n_sectors, n_years);
    for t = 1:n_years
        final_expenditure_share(:,:,t) = squeeze(final_exp_tr_matr(:,:,t)) ./ repmat(sum(final_exp_tr_matr(:,:,t),2), [1, n_sectors, 1]);
    end
    
    assert_all(size(final_expenditure_share) == [n_countries, n_sectors, n_years]);
    assert_all(final_expenditure_share(:) > 0);
    assert_all(abs(sum(final_expenditure_share, 2) - 1) < numerical_zero);
end

