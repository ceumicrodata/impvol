function equilibriumOut = equilibrium_algorithm(equilibriumInput)

% declare the constants used by most functions as global variables
global alpha beta theta xi kappa K B verbose lambda_w lambda_p c
% NOTE: see the definition of K in the documentation

theta = c.theta;

lambda_L = c.dampeningLaborLoop;
lambda_w = c.dampeningWageLoop;
lambda_p = c.dampeningPriceLoop;

% load data:
% aggregate labor: L (25 countries for 36 time periods)
% shocks: Z (24 sectors for 25 countries and 36 time periods)
% calibrated parameters: alpha, beta, theta, kappa, xi

% ln_rGDP = zeros(36, 25, 8);

% inFile = [inFolder, 'data_theta_', num2str(theta), '.mat'];
% load(inFile, 'baseline')

scenario = equilibriumInput.scenario;

B = equilibriumInput.B;
K = equilibriumInput.K;
alpha  = equilibriumInput.alpha;
beta = equilibriumInput.beta;
kappa =equilibriumInput.kappa;
xi = equilibriumInput.xi;

L = equilibriumInput.L;
Z = equilibriumInput.Z;
psi = equilibriumInput.psi;

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
Z_njt = Z; % compound shocks

[N, J, T] = size(Z_njt);

% initial guess on sector specific labor allocation
% L_njt = zeros(N, J, T);
% for t = 1:T
%    for j = 1:J
%        L_njt(:, j, t) = alpha(j, t) * L_nt(:, t);
% %        L_njt(:, j, t) = 1/24 * L_nt(:, t);
%    end
% end
L_njt = psi .* permute(repmat(L_nt, [1 1 J]), [1 3 2]);

%%%%% Outer loop: Search for sector specific labor allocation (L_njt)

% Store resource allocation in each iteration for convergence analysis
L_njt_iterations = zeros(N, J, T, 1);
L_njt_iterations(:, :, :, 1) = L_njt;

% technical values for the loop
outer_dif = c.outer_dif; % a big number
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
    tic
    [w_njt, w_nt, P_nt, P_njt, d] = get_wages(L_njt, L_nt, Z_njt, outer_iteration);
    toc
    
    [w_njt, w_nt, P_nt, P_njt] = fit_to_data(w_njt, w_nt, P_nt, P_njt);
    
    % calculate sectoral value added
    wL_njt = w_njt .* L_njt;
    
    % get wagebill shares (in logs)
    log_value_added_share = log(wL_njt ./ repmat(sum(wL_njt, 2), [1, J, 1]));
    
    % for each time period stack sectors (size: 25 x 1) on top of each other
    % these are going to be our dependent variables
    
    L_njt_new = update_resource_allocation_bp(log_value_added_share, L_nt);
    
    %L_njt_iterations(:, :, :, outer_iteration + 1) = L_njt_new;
    %save('results/iterations.mat', 'L_njt_iterations', 'L_nt')
    
    % calculate difference from last iteration
    outer_dif = max(abs((L_njt_new(:) - L_njt(:)) ./ (1 + L_njt(:))));
    
    % update sectoral labor allocations
    L_njt = lambda_L * L_njt_new + (1 - lambda_L) * L_njt;
    
    
    if (verbose > 0) && (mod(outer_iteration, c.outerPrintEvery) == 0) 
        fprintf('\n')
        fprintf('>>>> LABOR loop: Iteration %d, rel. diff.: %e \n',...
            outer_iteration, outer_dif)
%         fprintf('\n')
    end % if (verbose > 0)
end % while

if verbose > 0
    fprintf('Convergence in %d iterations. \n', outer_iteration)
    fprintf('####################  END OF LABOR LOOP   ####################\n')
    fprintf('\n')
end % if verbose > 0

%diary off

equilibriumOut = struct('scenario', scenario,...
                        'P_nt', P_nt, 'P_njt', P_njt,...
                        'w_nt', w_nt, 'w_njt', w_njt,...
                        'L_nt', L_nt, 'L_njt', L_njt,...
                        'd', d);

end