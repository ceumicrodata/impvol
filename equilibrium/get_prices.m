function [P_nj_new, inner_iteration] = get_prices(P_nj, D, t)

assert_all(P_nj>0);
[N, J] = size(P_nj);
assert_all(size(D)==[N, N, J]);

global c theta

lambda_p = c.dampening_price_loop;

% technical values
inner_tol = c.inner_tol;
inner_maxiter = c.inner_maxiter;
inner_convergence = 1;

inner_iteration = 0;
inner_dif = c.dif;


while inner_dif > inner_tol

    inner_iteration = inner_iteration + 1;
    if inner_iteration > inner_maxiter
        inner_convergence = 0;
        break
    end % if
    
    Pmtheta_nj = P_nj.^(-theta);
    
    Pmtheta_nj_update = real(price_eq_conditions(Pmtheta_nj, D, t));
    
    P_nj_update = Pmtheta_nj_update.^(- 1 / theta);
    
   
    step = P_nj_update - P_nj;

    inner_dif = norm(step(:)) / (1 + norm(P_nj(:)));

    P_nj = P_nj + lambda_p * step;    
    
end % while

if inner_convergence == 0
        fprintf('Maximum number of iterations (%d) exceeded in price loop. \n', inner_maxiter)
        fprintf('The difference is %e. \n', inner_dif)
%         error('No convergence in price loop.')
end % if convergence == 0

P_nj_new = P_nj;
% P_nj_new = Pmtheta_nj.^(- 1 / theta);
end