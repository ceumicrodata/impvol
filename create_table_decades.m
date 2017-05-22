%After running the full model, run this code to extract the decade averages
th = 4;
t1 = ['table_th_', num2str(th), '_lac_inf_bt'];
load_model_specifications(t1)

model = [t1, '_baseline'];
c = init_globals(model);
load([c.model_folder, 'equilibrium.mat'])
realGDP = equilibrium.L_nt .* equilibrium.w_nt ./ equilibrium.P_nt;
cd data_management
[~, cycle] = detrend_series(log(realGDP), c.filter_weights);
cd ..
cycle70 = cycle(:,2:8);
cycle80 = cycle(:,9:18);
cycle90 = cycle(:,19:28);
cycle00 = cycle(:,29:36);

model = [t1, '_nosectoral'];
c = init_globals(model);
load([c.model_folder, 'equilibrium.mat'])
realGDP = equilibrium.L_nt .* equilibrium.w_nt ./ equilibrium.P_nt;
cd data_management
[~, cycle] = detrend_series(log(realGDP), c.filter_weights);
cd ..
cycle70_nosectoral = cycle(:,2:8);
cycle80_nosectoral = cycle(:,9:18);
cycle90_nosectoral = cycle(:,19:28);
cycle00_nosectoral = cycle(:,29:36);

model = [t1, '_kappa1972'];
c = init_globals(model);
load([c.model_folder, 'equilibrium.mat'])
realGDP = equilibrium.L_nt .* equilibrium.w_nt ./ equilibrium.P_nt;
cd data_management
[~, cycle] = detrend_series(log(realGDP), c.filter_weights);
cd ..
cycle70_kappa1972 = cycle(:,2:8);
cycle80_kappa1972 = cycle(:,9:18);
cycle90_kappa1972 = cycle(:,19:28);
cycle00_kappa1972 = cycle(:,29:36);

model = [t1, '_kappa1972_nosectoral'];
c = init_globals(model);
load([c.model_folder, 'equilibrium.mat'])
realGDP = equilibrium.L_nt .* equilibrium.w_nt ./ equilibrium.P_nt;
cd data_management
[~, cycle] = detrend_series(log(realGDP), c.filter_weights);
cd ..
cycle70_kappa1972_nosectoral = cycle(:,2:8);
cycle80_kappa1972_nosectoral = cycle(:,9:18);
cycle90_kappa1972_nosectoral = cycle(:,19:28);
cycle00_kappa1972_nosectoral = cycle(:,29:36);

volatilities = [var(cycle70(:)); var(cycle80(:)); var(cycle90(:)); var(cycle00(:))];
volatilities_nosectoral = [var(cycle70_nosectoral(:)); var(cycle80_nosectoral(:)); var(cycle90_nosectoral(:)); var(cycle00_nosectoral(:))];
volatilities_kappa1972 = [var(cycle70_kappa1972(:)); var(cycle80_kappa1972(:)); var(cycle90_kappa1972(:)); var(cycle00_kappa1972(:))];
volatilities_kappa1972_nosectoral = [var(cycle70_kappa1972_nosectoral(:)); var(cycle80_kappa1972_nosectoral(:)); var(cycle90_kappa1972_nosectoral(:)); var(cycle00_kappa1972_nosectoral(:))];

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
