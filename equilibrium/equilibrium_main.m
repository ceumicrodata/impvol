function equilibrium_main

global c

% include_counterfactuals = c.include_counterfactuals;
load([c.model_folder, 'alg_inputs.mat']);
equilibrium_input = baseline;
%%
% Calculate baseline equilibrium
equilibrium = equilibrium_algorithm(equilibrium_input);
save([c.model_folder, 'equilibrium.mat'], 'equilibrium')

% Save nominal GDP and deflator in *.csv
country_names = importdata([c.data_folder_original, 'country_name.txt']);
fid = fopen([c.model_folder, 'ngdp.csv'], 'w');
fprintf(fid, '%s,',  'year', country_names{1:end-1});
fprintf(fid, '%s\n', country_names{end});
fclose(fid);
dlmwrite([c.model_folder, 'ngdp.csv'], [(1972:2007)', (equilibrium.L_nt .* equilibrium.w_nt)'], '-append', 'precision','%12.4f')


fid = fopen([c.model_folder, 'deflator.csv'], 'w');
fprintf(fid, '%s,',  'year', country_names{1:end-1});
fprintf(fid, '%s\n', country_names{end});
fclose(fid);
dlmwrite([c.model_folder, 'deflator.csv'], [(1972:2007)', equilibrium.P_nt'], '-append', 'precision','%10.4f')

