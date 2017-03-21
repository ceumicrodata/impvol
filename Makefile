CODES		= run_model.m init_globals.m data_management/calibrate_shocks.m data_management/create_counterfactual_scenarios.m equilibrium/equilibrium_main.m compute_volatilities.m
MATLAB		=  matlab
#/Applications/MATLAB_R2014b.app/bin/matlab
MATLAB_BEGIN	= $(MATLAB) -nodisplay -r "try, 
MATLAB_END		= , catch, exit(1), end, exit(0);"
COLUMNS 		= baseline nosectoral kappa1972 kappa1972_nosectoral
TABLES  		= table1 table2a table2b table3 table4

.PHONY: all

all: $(foreach table,$(TABLES),tables/$(table).csv)

tables/%.csv: create_table.m $(foreach specification,$(COLUMNS),models/%_$(specification)/data_rgdp_and_volatility.mat)
	echo $(filter-out $<,$(patsubst models/%/volatilities.csv,%,$^)) > $*.txt
	$(MATLAB_BEGIN) create_table('$*.txt', '$@')$(MATLAB_END)

models/%/data_rgdp_and_volatility.mat: model_specifications/%.mat $(CODES)
	$(MATLAB_BEGIN) run_model('$(notdir $(basename $<))')$(MATLAB_END)

$(foreach specification,$(COLUMNS),model_specifications/%_$(specification).mat): load_model_specifications.m model_specifications/%.csv
	$(MATLAB_BEGIN) load_model_specifications('$*')$(MATLAB_END)

