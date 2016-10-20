CODES		= run_model.m init_globals.m calibrate_shocks.m create_counterfactual_scenarios.m equilibrium_main.m
MATLAB		= matlab
# /Applications/MATLAB_R2014b.app/bin/matlab
MAP_BEGIN	= $(MATLAB) -nodisplay -r "try, run_model('
MAP_END		= '), catch, exit(1), end, exit(0);"
TABLE1		= table1_baseline table1_kappa1972 table1_nosectoral table1_kappa1972_nosectoral

.PHONY: all table1

all: table1

table1: $(foreach specification,$(TABLE1),models/$(specification)/volatilities.csv)

models/%/volatilities.csv: model_specifications/%.m 
	$(MAP_BEGIN)$(notdir $<)$(MAP_END)