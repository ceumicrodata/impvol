function equilibrium_out = equilibrium_algorithm(equilibrium_input)

% declare the constants used by most functions as global variables
global alpha beta theta xi kappa B verbose lambda_w lambda_p c gammas S
% NOTE: see the definition of K in the documentation

theta = c.theta;

lambda_L = c.dampening_labor_loop;
lambda_w = c.dampening_wage_loop;
lambda_p = c.dampening_price_loop;

% load data:
% aggregate labor: L (25 countries for 36 time periods)
% shocks: z (24 sectors for 25 countries and 36 time periods)
% calibrated parameters: alpha, beta, theta, kappa, xi

% ln_rGDP = zeros(36, 25, 8);

% in_file = [in_folder, 'data_theta_', num2str(theta), '.mat'];
% load(in_file, 'baseline')

scenario = equilibrium_input.scenario;

B = equilibrium_input.B;
% K = equilibrium_input.K;
alpha  = equilibrium_input.alpha;
beta = equilibrium_input.beta;
kappa =equilibrium_input.kappa;
xi = equilibrium_input.xi;
gammas = equilibrium_input.gammas;
S = equilibrium_input.trade_balance;

L = equilibrium_input.L;
z = equilibrium_input.z;
psi = equilibrium_input.psi;

% verbose is a global to control the amount of written output
% 0: no output
% 1: iteration and convergence info on labor loop
% 2: 1 + convergence info on wage loop
% 3: 2 + iteration info on wage loop 
% 4: 3 + convergence info on price loop 
verbose = c.verbosity;

% save output
% if verbose > 0
    %diary('covergence_log.txt')
    %diary on
% end % if verbose > 0


% For readability throughout the code all variables have subindices consistent with the 
% main article (Jan 2014 version)
% n: countries
% j: sectors
% t: time periods
% This makes it easy to infer the content and the dimension of each variable immediately.
L_nt = L; % total equipped labor
z_njt = z; % compound shocks

[N, J, T] = size(z_njt);

% initial guess on sector specific labor allocation
% L_njt = zeros(N, J, T);
% for t = 1:T
%    for j = 1:J
%        L_njt(:, j, t) = alpha(j, t) * L_nt(:, t);
% %        L_njt(:, j, t) = 1/24 * L_nt(:, t);
%    end
% end

L_share_njt = psi;

L_nt_full = permute(repmat(L_nt, [1 1 J]), [1 3 2]);

L_njt = L_share_njt .* L_nt_full;

%%
%%%%% Outer loop: Search for sector specific labor allocation (L_njt)

% Store resource allocation in each iteration for convergence analysis
% L_njt_iterations = zeros(N, J, T, 1);
% L_njt_iterations(:, :, :, 1) = L_njt;

% technical values for the loop
outer_dif = c.dif; % a big number
outer_tol = c.outer_tol; % set the convergence tolerance
outer_maxiter = c.outer_maxiter; % limit the maximum number of iterations

outer_iteration = 0; % set current iteration to zero

if verbose > 0
    fprintf([scenario, '\n'])
    fprintf('\n')
    fprintf('#######################   LABOR LOOP   #######################\n')
    fprintf('\n')
end % if verbose > 0

while outer_dif > outer_tol
    outer_iteration = outer_iteration + 1;
    if outer_iteration > outer_maxiter
        fprintf('Maximum number of iterations (%d) exceeded in labor loop.\n', outer_maxiter)
        fprintf('The difference is %e \n', outer_dif)
        error('No convergence in labor loop.')
    end % if iteration > maxiter
        
    % get sectoral wages, aggregate prices and aggregate wages that correspond
    % to the current value of sectoral labor allocation
%     tic
    [w_njt, w_nt, P_nt, P_njt, d] = get_wages(L_njt, L_nt, z_njt, outer_iteration);
%     toc
    
%     [w_njt, w_nt, P_nt, P_njt] = fit_to_data(w_njt, w_nt, P_nt, P_njt);
    
    % calculate sectoral value added
    wL_njt = w_njt .* L_njt;
    
    % get wagebill shares (in logs)
    log_value_added_share = log(wL_njt ./ repmat(sum(wL_njt, 2), [1, J, 1]));
    
    % wage gap in percentage points (see labourrel.tex)
    sectoral_wage_gap = (repmat(mean(w_njt, 2), [1, J, 1])-w_njt).*L_nt_full./repmat(sum(wL_njt, 2), [1, J, 1]);

    % for each time period stack sectors (size: 25 x 1) on top of each other
    % these are going to be our dependent variables
    
%     L_njt_new = update_resource_allocation_bp(log_value_added_share, L_nt);
    
    L_share_njt_new = update_resource_allocation_bp(log_value_added_share, sectoral_wage_gap);
    

    %L_njt_iterations(:, :, :, outer_iteration + 1) = L_njt_new;
    %save('results/iterations.mat', 'L_njt_iterations', 'L_nt')
       
    step = L_share_njt_new - L_share_njt;
    
    % calculate difference from last iteration
    outer_dif = norm(step(:)) / (1 + norm(L_share_njt(:)));
    
    % update sectoral labor allocations
%     L_njt = lambda_L * L_njt_new + (1 - lambda_L) * L_njt;
    L_share_njt = L_share_njt + 1 * step;
    
    L_njt = L_share_njt .* L_nt_full;
    
    if (verbose > 0) && (mod(outer_iteration, c.outer_print_every) == 0)
        fprintf('\n')
        fprintf('>>>> LABOR loop: Iteration %d, rel. diff.: %e \n',...
            outer_iteration, outer_dif)
        fprintf('Mean sectoral wage gap: %e \n', mean(sectoral_wage_gap)) 
%         fprintf('\n')
    end % if (verbose > 0)
end % while

if verbose > 0
    fprintf('Convergence in %d iterations. \n', outer_iteration)
    fprintf('####################  END OF LABOR LOOP   ####################\n')
    fprintf('\n')
end % if verbose > 0

%diary off

equilibrium_out = struct('scenario', scenario,...
                        'P_nt', P_nt, 'P_njt', P_njt,...
                        'w_nt', w_nt, 'w_njt', w_njt,...
                        'L_nt', L_nt, 'L_njt', L_njt,...
                        'd', d);

end