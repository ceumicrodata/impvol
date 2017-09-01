local column1 table_th_4_lac_inf_bt_baseline
local column2 table_th_4_lac_inf_bt_kappa1972
local column3 table_th_4_lac_inf_bt_nosectoral
local column4 table_th_4_lac_inf_bt_kappa1972_nosectoral

tempfile ngdp table
clear
gen str country=""
gen year=.
save `table', replace emptyok

forval i=1/4 {
	import delimited `column`i''/ngdp.csv, clear
	ren * ngdp*
	ren ngdpyear year
	reshape long ngdp, i(year) j(country) string
	save `ngdp', replace

	import delimited `column`i''/deflator.csv, clear
	ren * deflator*
	ren deflatoryear year
	reshape long deflator, i(year) j(country) string
	merge 1:1 country year using `ngdp'
	drop _m
	
	gen log_real_gdp_`i' = log(ngdp/deflator)
	keep country year log*
	
	merge 1:1 country year using `table'
	drop _m
	save `table', replace
}

gen decade = int(year/10)*10
egen ccode = group(country)
tsset ccode year

forval i=1/4 {
	gen growth`i' = D.log_real_gdp_`i'
}
collapse (sd) growth1 growth2 growth3 growth4, by(country decade)

forval i=1/4 {
	gen var`i' = growth`i'^2
}

gen total_change = 100*(var1-var2)/var2
gen nosectoral_change = 100*(var3-var4)/var2
gen difference = total_change-nosectoral_change

collapse (mean) total_change nosectoral_change difference, by(decade)

export delimited ../tables/table_7.csv, replace
