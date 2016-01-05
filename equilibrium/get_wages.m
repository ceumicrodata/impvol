function [w_njt, w_nt, P_nt, P_njt, d_mnjt] = get_wages(L_njt, L_nt, Z_njt, outer_iteration)
% This function calculates equilibrium sector-specific wages for any given set of
% sector-specific labor allocation.
% The corresponding equilibrium aggregate wages and aggregate prices are also returned along the wages.

global verbose lambda_w c

[N, J, T] = size(L_njt);

P_nt = zeros(N, T); % Aggregate prices
P_njt = zeros(N, J, T); % Sector specific prices
w_nt = zeros(N, T); % Aggregate wages
w_njt = zeros(N, J, T); % Sector specific wages
d_mnjt = zeros(N, N, J, T);

for t = 1:T 
    
    if verbose >= 3
        fprintf('-------  WAGE LOOP, period %d, labor iteration %d  --------\n', t, outer_iteration)
    end
    
    % Get values from appropriate time periods
    L_nj = L_njt(:, :, t);
    L_n = L_nt(:, t);
    Z_nj = Z_njt(:, :, t);
    
    % Initial guess on sector-specific wages
    % Use previous period values in later periods
    % Note: Subsequent dampening parameters will start with values of the
    % previous dampening parameters
    if t > 1
        w_nj_start = w_njt(:, :, t - 1);
        P_n = P_nt(:, t - 1);
    else 
        w_nj_start = 1 * ones(N, J);
        P_n = 80 * ones(N, 1);
    end % if t == T0
    
    w_nj = w_nj_start;
    
    % technical values for the loop
    middle_tol = c.middle_tol; 
    middle_maxiter = c.middle_maxiter; 
    middle_convergence = 0;
    
    damp_step = 0;
    for lambda = lambda_w
        damp_step = damp_step + 1;
        middle_dif = c.middle_dif;
        middle_iteration = 0; % set current iteration to zero
        
        % update wage until convergence or maximum number of iterations
        while middle_dif > middle_tol       
            middle_iteration = middle_iteration + 1;
            
            if middle_iteration > middle_maxiter
                break
            end
            
            if middle_iteration > middle_maxiter
                fprintf('Maximum number of iterations (%d) exceeded in wage loop.\n', middle_maxiter)
                fprintf('The difference is %e \n', middle_dif)
                error('No convergence in wage loop.')
            end
            
            % calculate new sector specific wages (and associated prices)
            % based on current values
            [w_nj_new, P_n, price_iterations, price_lambda, P_nj, d] = ...
                wage_update(w_nj, L_nj, L_n, Z_nj, P_n, t);
            
            middle_dif = max(abs((w_nj_new(:) - w_nj(:)) ./ (1 + w_nj(:))));
            
            % update current values
            w_nj = lambda * w_nj_new + (1 - lambda) * w_nj;
            
            if (verbose >= 3) && (mod(middle_iteration, c.middlePrintEvery) == 0)
                fprintf('    Wage iteration %d, difference is %e\n', middle_iteration, middle_dif)
                if (verbose == 4)
                    fprintf('        Price loop converged in %d iterations and with %1.2f as the dampening parameter.\n', price_iterations, price_lambda)
                    fprintf('        Mean aggregate price is %e \n', mean(P_n))
                    fprintf('\n')
                end % if (verbose == 4)
            end % (verbose == 3) && (mod(iteration, 10) == 0)
        end % while dif > tol
        
        if middle_dif < middle_tol
            middle_convergence = 1;
            break
        end %if
        
    end % for lambda = lambda_w
    
    
    
    if middle_convergence == 0
        fprintf('Maximum number of iterations (%d) exceeded in wage loop. \n', middle_maxiter)
        fprintf('The difference is %e. \n', middle_dif)
        error('No convergence in wage loop.')
    end % if convergence == 0
    
    % calculate aggregate wages
    wL_n = sum(w_nj .* L_nj, 2);
    w_n = wL_n ./ L_n;  
    
    % Save results from period t to output
    P_nt(:, t) = P_n;
    P_njt(:, :, t) = P_nj;
    w_nt(:, t) = w_n;
    w_njt(:, :, t) = w_nj;
    d_mnjt(:, :, :, t) = d;
    
    
    if verbose == 2
        fprintf('    WAGE loop, period %d, labor iteration %d: Convergence in %d iterations. Dampening: %1.2f \n',...
            t, outer_iteration,...
            (damp_step - 1) * middle_maxiter + middle_iteration, lambda)
    elseif verbose >= 3
        fprintf('    Convergence in %d iterations and with %1.1f as the dampening parameter.\n', middle_iteration, lambda)
        %fprintf('The minimum aggregate wage is %f \n', min(w_n))
        %fprintf('The maximum aggregate wage is %f \n', max(w_n))
        fprintf('----  END OF WAGE LOOP, period %d, labor iteration %d  ----\n', t, outer_iteration)
        fprintf('\n')
    end % if verbose == 2

end % for t = 1:T

end
