CODES		= run_model.m init_globals.m calibrate_shocks.m create_counterfactual_scenarios.m equilibrium_main.m
MATLAB		= matlab
# /Applications/MATLAB_R2014b.app/bin/matlab
MAP_BEGIN	= $(MATLAB) -nodisplay -r "try, run_model('
MAP_END		= '), catch, exit(1), end, exit(0);"

.PHONY: all table1

all: table1

table1: models/table1_baseline/volatilities.csv models/table1_kappa1972/volatilities.csv

models/%/volatilities.csv: model_specifications/%.m 
	$(MAP_BEGIN)$(notdir $<)$(MAP_END)