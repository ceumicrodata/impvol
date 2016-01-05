function [out, inner_iteration, lambda] = get_prices(D, Ptheta_n_previous, t)

global lambda_p c

% initialize price vector
Ptheta_n_start = Ptheta_n_previous;

% technical values
inner_tol = c.inner_tol;
inner_maxiter = c.inner_maxiter;
inner_convergence = 0;

for lambda = lambda_p
    Ptheta_n = Ptheta_n_start;
    inner_iteration = 0;
    inner_dif = c.inner_dif;
    % loop contraction (?) until convergence or maximum number of iterations
    while inner_dif > inner_tol
        inner_iteration = inner_iteration + 1;
        if inner_iteration > inner_maxiter
            break
        end % if
        Ptheta_new = price_contraction(Ptheta_n, D, t);
        inner_dif = max(abs((Ptheta_new - Ptheta_n) ./ (1 + Ptheta_n)));
        Ptheta_n = lambda * Ptheta_new + (1 - lambda) * Ptheta_n;
    end % while
    if inner_dif < inner_tol
        inner_convergence = 1;
        break
    end %if
end % for

if inner_convergence == 0
        fprintf('Maximum number of iterations (%d) exceeded in price loop. \n', inner_maxiter)
        fprintf('The difference is %e. \n', inner_dif)
        error('No convergence in price loop.')
end % if convergence == 0

out = Ptheta_n; 
end