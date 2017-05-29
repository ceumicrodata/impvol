%After running the full model, run this code to extract the decade averages
th = 4;
t1 = ['table_th_', num2str(th), '_lac_inf_bt'];
load_model_specifications(t1)

model = [t1, '_baseline'];
c = init_globals(model);
load([c.model_folder, 'equilibrium.mat'])
realGDP = equilibrium.L_nt .* equilibrium.w_nt ./ equilibrium.P_nt;
cd data_management
[~, cycle_baseline] = detrend_series(log(realGDP), c.filter_weights);
cd ..
cycle_baseline = mean(cycle_baseline, 1);

model = [t1, '_nosectoral'];
c = init_globals(model);
load([c.model_folder, 'equilibrium.mat'])
realGDP = equilibrium.L_nt .* equilibrium.w_nt ./ equilibrium.P_nt;
cd data_management
[~, cycle_nosectoral] = detrend_series(log(realGDP), c.filter_weights);
cd ..
cycle_nosectoral = mean(cycle_nosectoral, 1);

model = [t1, '_kappa1972'];
c = init_globals(model);
load([c.model_folder, 'equilibrium.mat'])
realGDP = equilibrium.L_nt .* equilibrium.w_nt ./ equilibrium.P_nt;
cd data_management
[~, cycle_kappa1972] = detrend_series(log(realGDP), c.filter_weights);
cd ..
cycle_kappa1972 = mean(cycle_kappa1972, 1);

model = [t1, '_kappa1972_nosectoral'];
c = init_globals(model);
load([c.model_folder, 'equilibrium.mat'])
realGDP = equilibrium.L_nt .* equilibrium.w_nt ./ equilibrium.P_nt;
cd data_management
[~, cycle_kappa1972_nosectoral] = detrend_series(log(realGDP), c.filter_weights);
cd ..
cycle_kappa1972_nosectoral = mean(cycle_kappa1972_nosectoral, 1);

volatilities = [var(cycle_baseline(1,2:8)); var(cycle_baseline(1,9:18)); var(cycle_baseline(1,19:28)); var(cycle_baseline(1,29:36))];
volatilities_nosectoral = [var(cycle_nosectoral(1,2:8)); var(cycle_nosectoral(1,9:18)); var(cycle_nosectoral(1,19:28)); var(cycle_nosectoral(1,29:36))];
volatilities_kappa1972 = [var(cycle_kappa1972(1,2:8)); var(cycle_kappa1972(1,9:18)); var(cycle_kappa1972(1,19:28)); var(cycle_kappa1972(1,29:36))];
volatilities_kappa1972_nosectoral = [var(cycle_kappa1972_nosectoral(1,2:8)); var(cycle_kappa1972_nosectoral(1,9:18)); var(cycle_kappa1972_nosectoral(1,19:28)); var(cycle_kappa1972_nosectoral(1,29:36))];

vol = [volatilities, volatilities_nosectoral, volatilities_kappa1972, volatilities_kappa1972_nosectoral];
decades = cellstr({'1973-1979'; '1980-1989'; '1990-1999'; '2000-2008'});

model = [t1, '_baseline'];
c = init_globals(model);
fid = fopen([c.model_folder, 'decade_vol.csv'], 'w');
fprintf(fid, '%25s, %25s, %25s, %25s, %25s\n', 'decade', 'baseline', 'nosectoral', 'kappa1972', 'kappa1972_nosectoral');

for i = 1:size(vol,1)
	fprintf(fid, '%25s, %12.10f, %12.10f, %12.10f, %12.10f\n', decades{i}, vol(i,1), vol(i,2), vol(i,3), vol(i,4));
end
fclose(fid);
