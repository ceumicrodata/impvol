function [w_njt, w_nt, P_nt, P_njt] = get_wages(equilibrium_input, L_njt, L_nt, z_njt, outer_iteration, w_njt, P_njt)
% This function calculates equilibrium sector-specific wages for any given set of
% sector-specific labor allocation.
% The corresponding equilibrium aggregate wages and aggregate prices are also returned along the wages.

global verbose c alpha beta gammas S

% lambda_w = c.dampening_wage_loop;

[N, J, T] = size(L_njt);

% P_nt = zeros(N, T); % Aggregate prices
% P_njt = zeros(N, J, T); % Sector specific prices
% w_nt = zeros(N, T); % Aggregate wages
% w_njt = zeros(N, J, T); % Sector specific wages
% d_mnjt = zeros(N, N, J, T);


for t = 1:T 
    
    if verbose >= 3
        fprintf('-------  WAGE LOOP, period %d, labor iteration %d  --------\n', t, outer_iteration)
    end
    
    % Get values from appropriate time periods
    L_nj = L_njt(:, :, t);
    L_n = L_nt(:, t);
    z_nj = z_njt(:, :, t);
    
    % Initial guess on sector-specific wages
    % Use previous period values in later periods
    % Note: Subsequent dampening parameters will start with values of the
    % previous dampening parameters
    if (t > 1) && (outer_iteration == 1) 
        w_nj_start = w_njt(:, :, t - 1);
        P_nj = P_njt(:, :, t - 1);
    else 
        w_nj_start = w_njt(:, :, t);
        P_nj = P_njt(:, :, t);
    end % if
%     if t > 1
%         w_nj_start = w_njt(:, :, t - 1);
%         P_nj = P_njt(:, :, t - 1);
%     else 
%         w_nj_start = 1e-1 * ones(N, J);
%         P_nj = 1e-1 * ones(N, J);
% %         w_nj_start = 10 * rand(N, J);
% %         P_nj = 10 * rand(N, J);
%     end % if t > 1
    
    w_nj = w_nj_start;
    
    % technical values for the loop
    middle_tol = c.middle_tol; 
    middle_maxiter = c.middle_maxiter; 
    middle_convergence = 0;
    
    middle_dif = c.dif;
    middle_iteration = 0; % set current iteration to zero
    
    load([c.model_folder, '/data_rgdp_and_volatility.mat'], 'va_total', 'p_base')
    
    va_t_total = sum(va_total(t, :));
%     va_t_us = va_total(t, 25);
    
    va_to_fit = va_t_total;
%     va_to_fit = va_t_us;
    p_to_fit = p_base(t);
    
    
    B_gamma = kron(eye(N), sparse(gammas(:, :, t)));
    B_beta = kron(eye(N), sparse(repmat(beta', [J, 1])));
    %% FIXME: D_alpha is now country and time varying
    D_alpha = diag(sparse(repmat(alpha(:, t), [N, 1])));
    S_full = kron(S(:, t), ones([J, 1]));
    beta_full = repmat(beta, [N, 1]);

    
    % update wage until convergence or maximum number of iterations
    while middle_dif > middle_tol
        middle_iteration = middle_iteration + 1;
        
%         if middle_iteration > middle_maxiter
%             break
%         end
        
        if middle_iteration > middle_maxiter
            fprintf('Maximum number of iterations (%d) exceeded in wage loop.\n', middle_maxiter)
            fprintf('The difference is %e \n', middle_dif)
            fprintf('The mean is %e. \n', mean(w_nj(:)))
%             error('No convergence in wage loop.')
            break
        end
        
        % calculate new sector specific wages (and associated prices)
        % based on current values
        [w_nj_new, P_nj, price_iterations, P_n] = wage_update(equilibrium_input, w_nj, L_nj, z_nj, P_nj, t, va_to_fit, p_to_fit, B_gamma, B_beta, D_alpha, S_full, beta_full);
                
        step = w_nj_new - w_nj;
        
        middle_dif = norm(step(:)) / (1 + norm(w_nj(:)));   
     
%         w_nj = w_nj + 0.2 * step;
        
        % update current values
        
%         w_nj = w_nj + (lambda_w/2 + lambda_w * rand()) * step;
        
        mm = mod(floor(middle_iteration / 100), 5);
        if mm == 0
            w_nj = w_nj + 0.15 * step;
        elseif mm == 1
            w_nj = w_nj + 0.1 * step;
        elseif mm == 2
            w_nj = w_nj + 0.05 * step;
        elseif mm == 3
            w_nj = w_nj + 0.025 * step;
        elseif mm == 4
            w_nj = w_nj + 0.01 * step;
        end %if        
%             
         
        if (verbose >= 3) && (mod(middle_iteration, c.middle_print_every) == 0)
            fprintf('    Wage iteration %d, difference is %e\n', middle_iteration, middle_dif)
            if (verbose == 4)
                fprintf('        Price loop converged in %d iterations.\n', price_iterations)
                fprintf('        Mean aggregate price is %e \n', mean(P_n))
                fprintf('\n')
            end % if (verbose == 4)
        end % (verbose == 3) && (mod(iteration, 10) == 0)
    end % while dif > tol
    
    if middle_dif < middle_tol
        middle_convergence = 1;
    end %if
  
    
%     if middle_convergence == 0
%         fprintf('Maximum number of iterations (%d) exceeded in wage loop. \n', middle_maxiter)
%         fprintf('The difference is %e. \n', middle_dif) 
% %         error('No convergence in wage loop.')
%     end % if convergence == 0
    
    % calculate aggregate wages
    wL_n = sum(w_nj .* L_nj, 2);
    w_n = wL_n ./ L_n;  
    
    % Save results from period t to output
    P_nt(:, t) = P_n;
    P_njt(:, :, t) = P_nj;
    w_nt(:, t) = w_n;
    w_njt(:, :, t) = w_nj;    
    
    if verbose == 2
        fprintf('    WAGE loop, period %d, labor iteration %d: Convergence in %d iterations.\n',...
            t, outer_iteration, middle_iteration)
    elseif verbose >= 3
        fprintf('    Convergence in %d iterations.\n', middle_iteration)
        %fprintf('The minimum aggregate wage is %f \n', min(w_n))
        %fprintf('The maximum aggregate wage is %f \n', max(w_n))
        fprintf('----  END OF WAGE LOOP, period %d, labor iteration %d  ----\n', t, outer_iteration)
        fprintf('\n')
    end % if verbose == 2

end % for t = 1:T

end