CODES		= run_model.m init_globals.m calibrate_shocks.m create_counterfactual_scenarios.m equilibrium_main.m compute_volatilities.m
MATLAB		=  matlab
#/Applications/MATLAB_R2014b.app/bin/matlab
MATLAB_BEGIN	= $(MATLAB) -nodisplay -r "try, 
MATLAB_END		= , catch, exit(1), end, exit(0);"
COLUMNS 		= baseline nosectoral kappa1972 kappa1972_nosectoral

.PHONY: all

all: tables/table1.csv tables/table2a.csv tables/table3.csv

tables/%.csv: create_table.m $(foreach specification,$(COLUMNS),models/%_$(specification)/data_rgdp_and_volatility.mat)
	echo $(filter-out $<,$(patsubst models/%/volatilities.csv,%,$^)) > $*.txt
	$(MATLAB_BEGIN) create_table('$*.txt', '$@')$(MATLAB_END)

models/%/data_rgdp_and_volatility.mat: model_specifications/%.m $(CODES)
	$(MATLAB_BEGIN) run_model('$(notdir $(basename $<))')$(MATLAB_END)