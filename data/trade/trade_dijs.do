* Here I calculate the surspluses and the dij's from the aggregated trade dataset

clear
set mem 1000m
*cd "/Users/federicorossi/Dropbox/impvol/data/federico/Stata Files/"
cd "C:\Users\rossif1\Dropbox\impvol\data\federico\Stata Files"

use "Trade/trade_aggr.dta"

/*
*********************************** Re-organize dataset ***********************************
****Here I consider info on imports and exports separately
* Consider import info first

drop if Trade_Flow_Code==2

rename Value_ALL import_AfromB_all
rename Value_Agriculture import_AfromB_agriculture
rename Value_Manufacturing import_AfromB_manufacturing

global isic_num 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36
foreach i of global isic_num{
rename Value_Manufacturing_`i' import_AfromB_manufacturing_`i' 
}

drop Trade_Flow_Code

sort countrya year countryb

save "Trade/trade_imp_info", replace

* Consider export info 

use "Trade/trade_aggr.dta", clear

drop if Trade_Flow_Code==1

//I reverse the "direction" of exports (and switch the names of the countries variables accordingly)

rename countryb country
rename countrya countryb
rename country countrya

rename Value_ALL export_BtoA_all
rename Value_Agriculture export_BtoA_agriculture
rename Value_Manufacturing export_BtoA_manufacturing

foreach i of global isic_num{
rename Value_Manufacturing_`i' export_BtoA_manufacturing_`i' 
}


drop Trade_Flow_Code

sort countrya year countryb

save "Trade/export_imp_info", replace

merge countrya year countryb using "Trade/trade_imp_info"

* Use import data when possible
gen import_agriculture=import_AfromB_agriculture
gen import_manufacturing=import_AfromB_manufacturing
foreach i of global isic_num{
gen import_manufacturing_`i'=import_AfromB_manufacturing_`i'
}


gen import_all=import_AfromB_all

* Integrate with export data when import data not available
replace import_agriculture=export_BtoA_agriculture if import_agriculture==.
replace import_manufacturing=export_BtoA_manufacturing if import_manufacturing==.
foreach i of global isic_num{
replace import_manufacturing_`i'=export_BtoA_manufacturing_`i' if import_manufacturing_`i'==.
}


replace import_all=export_BtoA_all if import_all==.

label variable import_agriculture "Imports in agriculture of countrya from countryb"
label variable import_manufacturing "Imports in manufacturing of countrya from countryb"
label variable import_all "Total imports of countrya from countryb"

drop export_BtoA* import_AfromB* _merge

sort countrya year countryb

save "Trade/trade_flows_pre_imputation", replace
*/

****************** Rename Value_ to import_ **************************************
rename Value_Agriculture import_agriculture
rename Value_Manufacturing* import_manufacturing*
rename Value_ALL import_all



*************************************Surpluses********************************************

egen totimport_agriculture = sum(import_agriculture), by (countrya year) 
egen totimport_manufacturing = sum(import_manufacturing), by (countrya year)
foreach i of global isic_num{
egen totimport_manufacturing_`i' = sum(import_manufacturing_`i'), by (countrya year)
}

 
egen totimport_all = sum(import_all), by (countrya year) 

label variable totimport_agriculture "Total imports in agriculture of countrya"
label variable totimport_manufacturing "Total imports in manufacturing of countrya"
label variable totimport_all "Total imports of countrya"

sort countrya year countryb

save "Trade/temp", replace

**Now consider exports separately

drop totimport*

egen totexport_agriculture = sum(import_agriculture), by (countryb year) 
egen totexport_manufacturing = sum(import_manufacturing), by (countryb year) 
foreach i of global isic_num{
egen totexport_manufacturing_`i' = sum(import_manufacturing_`i'), by (countryb year) 
}
egen totexport_all = sum(import_all), by (countryb year) 

keep countryb year totexport_agriculture totexport_manufacturing* totexport_all 

rename countryb countrya

duplicates drop countrya year totexport_agriculture totexport_manufacturing* totexport_all, force
sort countrya year 

merge countrya year using "Trade/temp"
tab _m
drop _m

sort countrya year countryb

label variable totexport_agriculture "Total exports in agriculture of countrya"
label variable totexport_manufacturing "Total exports in manufacturing of countrya"
label variable totexport_all "Total exports of countrya"


**** Now compute surpluses

gen surplus_agriculture = totexport_agriculture - totimport_agriculture
gen surplus_manufacturing = totexport_manufacturing - totimport_manufacturing
foreach i of global isic_num{
gen surplus_manufacturing_`i' = totexport_manufacturing_`i' - totimport_manufacturing_`i'
}
gen surplus_all = totexport_all - totimport_all

label variable surplus_agriculture "Total surplus (exports-imports) in agriculture of countrya"
label variable surplus_manufacturing "Total surplus (exports-imports) in manufacturing of countrya"
label variable surplus_all "Total surplus (exports-imports) of countrya"

sort countrya year countryb

save "Trade/trade_surpluses", replace

***********************************************************************************************
************ Merge with output ***************************************

use "Output/Output_usd_aggr", clear

sort country year

rename country countrya

merge countrya year using "Trade/trade_surpluses"

/*
***Convert units (from millions of US dollars to US dollars)
gen agriculture=agriculture_usd*1000000
gen manufacturing=manufacturing_usd*1000000
foreach i of global isic_num{
gen manufacturing_`i'=manufacturing_usd_`i'*1000000
}
gen services=services_usd*1000000
*/
**don't convert anymore since units changed in trade data; rename instead
rename agriculture_usd agriculture
rename manufacturing_usd* manufacturing*
rename services_usd services

*drop agriculture_usd manufacturing_usd* services_usd 
drop _m
drop if countryb==""


label variable agriculture "Output in agriculture (millions of US dollars)"
label variable manufacturing "Output in manufacturing (millions of US dollars)"
label variable services "Output in services (millions of US dollars)"

sort countrya year countryb

****************************** Compute dij's **************************************

gen dij_agriculture=(import_agriculture)/(agriculture - surplus_agriculture)
gen dij_manufacturing=(import_manufacturing)/(manufacturing - surplus_manufacturing)
foreach i of global isic_num{
gen dij_manufacturing_`i'=(import_manufacturing_`i')/(manufacturing_`i' - surplus_manufacturing_`i')
}
egen di_others_agriculture = sum (dij_agriculture), by (countrya year)
egen di_others_manufacturing = sum (dij_manufacturing), by (countrya year)
foreach i of global isic_num{
egen di_others_manufacturing_`i' = sum (dij_manufacturing_`i'), by (countrya year)
}
gen dii_agriculture= 1 - di_others_agriculture
gen dii_manufacturing= 1 - di_others_manufacturing
foreach i of global isic_num{
gen dii_manufacturing_`i'= 1 - di_others_manufacturing_`i'
}

label variable dij_agriculture "d_ij in agriculture for countrya with respect to countryb"
label variable dii_agriculture "d_ii in agriculture for countrya"
label variable dij_manufacturing "d_ij in manufacturing for countrya with respect to countryb"
label variable dii_manufacturing "d_ii in manufacturing for countrya"

drop di_others_agriculture di_others_manufacturing* 

drop totexport_agriculture totexport_manufacturing* totexport_all import_agriculture import_manufacturing* import_all totimport_agriculture totimport_manufacturing* totimport_all agriculture manufacturing* services 

sort countrya year countryb

save "Trade/trade_dij", replace
