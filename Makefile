CODES		= create_table.m run_model.m init_globals.m data_management/calibrate_shocks.m data_management/create_counterfactual_scenarios.m equilibrium/equilibrium_main.m equilibrium/equilibrium_algorithm.m equilibrium/get_wages.m equilibrium/wage_update.m compute_volatilities.m
MATLAB		=  matlab
MATLAB_BEGIN	= $(MATLAB) -nodisplay -r "try, 
MATLAB_END		= , catch, exit(1), end, exit(0);"
COLUMNS 		= baseline nosectoral kappa1972 kappa1972_nosectoral
TABLES  		= 1 2 3 4 5 6 

tables_ready_to_commit: $(foreach table,$(TABLES),tables/table_$(table).csv)
	git add $?
	git commit -m "autocommit: $?"
	git push
	touch $@

tables/table_1.csv: table_1.m  $(CODES)
	$(MATLAB_BEGIN) $(basename $<) $(MATLAB_END)	
	cat tables/table_th_4_lac_inf_bt.csv > $@

tables/table_2.csv: table_2.m $(CODES)
	$(MATLAB_BEGIN) $(basename $<) $(MATLAB_END)	
	cat tables/table_th_4_lac_inf.csv > $@

tables/table_3.csv: table_3.m $(CODES)
	$(MATLAB_BEGIN) $(basename $<) $(MATLAB_END)	
	cat tables/table_th_2_lac_inf_bt.csv > $@
	cat tables/table_th_8_lac_inf_bt.csv >> $@

tables/table_4.csv: table_4.m $(CODES)
	$(MATLAB_BEGIN) $(basename $<) $(MATLAB_END)	
	cat tables/table_th_4_lac_inf_bt_noio.csv > $@

tables/table_5.csv: table_5.m $(CODES)
	$(MATLAB_BEGIN) $(basename $<) $(MATLAB_END)	
	cat tables/table_th_4_lac_2000_bt.csv > $@
	cat tables/table_th_4_lac_1000_bt.csv >> $@
	cat tables/table_th_4_lac_500_bt.csv >> $@

tables/table_6.csv: table_6.m $(CODES)
	$(MATLAB_BEGIN) $(basename $<) $(MATLAB_END)	
	cat tables/table_th_4_lac_inf_bt_china_1.csv > $@
	cat tables/table_th_4_lac_inf_bt_china_2.csv >> $@