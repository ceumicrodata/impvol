import delimited import_share.txt
edit
ren v1 year 
ren v2 importer
ren v3 exporter
edit if abs(year-1985)<5 & importer==3
egen average_share = rmean(v*)
scatter average_share year  if abs(year-1985)<5 & importer==3
lowess average_share year  if abs(year-1985)<10 & importer==3
lowess average_share year  if abs(year-1985)<12 & importer==3, bw(0.1)
lowess average_share year  if abs(year-1985)<12 & importer==3, bw(0.3)
lowess average_share year  if abs(year-1985)<12 & exporter==3, bw(0.3)
lowess average_share year  if abs(year-1985)<10 & exporter==3, bw(0.3)
tw (lowess average_share year  if abs(year-1985)<12 & importer==3, bw(0.3))
tw (lowess average_share year  if abs(year-1985)<20 & importer==3, bw(0.3))
tw (lowess average_share year  if abs(year-1985)<20 & importer==3, bw(0.1))
tw (lowess average_share year  if abs(year-1985)<20 & importer==3, bw(0.2))
drop average_share 
reshape long v, i(year importer exporter ) j(sector)
replace sector = sector -3
edit
tw (lowess v year  if abs(year-1985)<20 & importer==3, bw(0.3))
 lowess v year  if abs(year-1985)<20 & importer==3, bw(0.3)
 gen byte negative = v<0
  lowess negative  year  if abs(year-1985)<20 & importer==3, bw(0.3)
  tw ( lowess negative  year  if abs(year-1985)<20 & importer==3, bw(0.3))
  tw ( lowess negative  year  if abs(year-1985)<20 & importer==16, bw(0.3))
  tw ( lowess negative  year  if abs(year-1985)<20 & importer==12, bw(0.3))
  tw ( lowess negative  year  if abs(year-1985)<20 & importer==9, bw(0.3))
  tw ( lowess negative  year  if abs(year-1985)<20 & importer==10, bw(0.3))
  tw ( lowess negative  year  if abs(year-1985)<20 & importer==17, bw(0.3))
  tw ( lowess negative  year  if abs(year-1985)<30 & importer==17, bw(0.3))
  tw ( lowess negative  year  if abs(year-1985)<30 & importer==3, bw(0.1))
  table year importer, c(mean negative )
  di 1/0.043478
  table year exporter, c(mean negative )
  table year sector, c(mean negative )
  table year sector if importer ==3, c(mean negative )
  table year sector if importer ==16, c(mean negative )
  table year sector if importer ==3, c(mean negative )
  table year if importer ==3 & sector==12, c(mean negative )
  table year if importer ==3 & sector==12, c(mean v )
  table year if importer ==3 & inlist(sector,11,12,13), c(mean v )
  table year sector if importer ==3 & inlist(sector,11,12,13), c(mean v )
  egen total_import_share = sum(v), by(importer sector year)
  table year sector if importer ==3 & inlist(sector,11,12,13), c(mean v mean total_import_share  )
  table year sector if importer ==3 & inlist(sector,11,12,13), c(mean total_import_share  )
  replace total_import_share = 0.999 if total_import_share >0.999 & total_import_share <0
  replace total_import_share = 0.999 if total_import_share >0.999 | total_import_share <0
  table year sector if importer ==3 & inlist(sector,11,12,13), c(mean total_import_share  )
  edit if importer ==3 & inlist(sector,11,12,13)
  egen sum_v = sum(v), by(importer sector year)
  gen share_within_imports = v/sum_v 
  su share_within_imports , d
  drop total_import_share 
  su sum_v , d
  gen domestic_per_total_import = 1/sum_v -1
  su domestic_per_total_import , d
  gen domestic_per_total_import_2 = domestic_per_total_import  
  replace domestic_per_total_import_2 = 0.001 if domestic_per_total_import_2 < 0.001
  gen v2 = share_within_imports / (1+domestic_per_total_import_2 )
  scatter v2 v
  table year sector if importer ==3 & inlist(sector,11,12,13), c(mean v mean v2)
  table year sector if importer ==16, c(mean v)
  table year sector if importer ==16 & inlist(sector,5), c(mean v mean v2)
  save explore_reasons_for_negative_import_shares
  
