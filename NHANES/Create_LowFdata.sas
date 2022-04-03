dm log "clear";
dm output "clear";

ods preferences;
ods html close;
ods html ;

libname pride 'C:\Users\brs380\OneDrive - Harvard University\Migrated-P-Drive\NHANES\fped';
libname nhanes 'C:\Users\brs380\OneDrive - Harvard University\Migrated-P-Drive\NHANES\input_dietdata'; 

proc contents data=pride.hei_tert1118 varnum; run;

data lowF;
	set pride.hei_tert1118;
if RIAGENDR=2 AND INDFMPIR<=1.3;
if indfmpir =. then delete; *include only participants with known PIL information
run; 

/*Raw weighted Frequency of diet foods */
proc surveyfreq data=lowF;
	cluster sdmvpsu;
	strata sdmvstra;
	weight dietwt8yr;
	tables bb1-bb29;
run; 

data pride.nhanes1118_lowFdietdata;
	set lowF;
run; 
/*export lowF to nhanes_lowFadultdata1Dec2021.csv */
