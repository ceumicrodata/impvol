function [P_nj_new, inner_iteration] = get_prices(P_nj, D, t)

global c theta

% technical values
inner_tol = c.inner_tol;
inner_maxiter = c.inner_maxiter;
inner_convergence = 1;

Pmtheta_nj = P_nj.^(-theta);
inner_iteration = 0;
inner_dif = c.dif;

while inner_dif > inner_tol
    inner_iteration = inner_iteration + 1;
    if inner_iteration > inner_maxiter
        inner_convergence = 0;
        break
    end % if
    
    Pmtheta_nj_update = price_eq_conditions(Pmtheta_nj, D, t);
      
    step = Pmtheta_nj_update - Pmtheta_nj;
   
    inner_dif = max(abs(step(:)));
     
    Pmtheta_nj_new = Pmtheta_nj + 1 * step;
    
    Pmtheta_nj = Pmtheta_nj_new;
    
end % while

if inner_convergence == 0
        fprintf('Maximum number of iterations (%d) exceeded in price loop. \n', inner_maxiter)
        fprintf('The difference is %e. \n', inner_dif)
        error('No convergence in price loop.')
end % if convergence == 0
 
P_nj_new = Pmtheta_nj.^(- 1 / theta);
end