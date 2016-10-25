function load_model_specifications(table_name)

fid = fopen(['model_specifications/', table_name, '.csv']);
model_names = textscan(fid, '%s', 5, 'Delimiter', ',');
fclose(fid);

specs_table = csvread(['model_specifications/',table_name, '.csv'], 1, 1);

for i = 1:4
    spec.model = [table_name, '_', model_names{1}{i + 1}];
    %% Model type
    % io links (bool)    
    spec.iol = specs_table(1, i);

    % Balanced trade (bool) 
    spec.bt = specs_table(2, i);

    % Labor adjustment cost (0: infinity, int: rho)
    spec.lac = specs_table(3, i);



    %% Fundamental parameters
    % theta in {2, 4, 8}
    spec.th = specs_table(4, i);
    
    % eta (int?)
    spec.et = specs_table(5, i);
    
    
    
    %% Conterfactual parameters
    % trade cost (0: actual calibrated, 1: 1972, 2: free trade)
    spec.tc = specs_table(6, i);
    
    % productivity shock type (0: actual calibrated,
    %                          1: no sectoral shocks,
    %                          2: no sectoral and residual shocks)
    spec.sh = specs_table(7, i);
    
    save(['model_specifications/', spec.model, '.mat'], 'spec')
end




