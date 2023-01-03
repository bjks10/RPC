/********************************************/
/* NHANES - FPED 							*/
/* Task: Create tertile variables 			*/
/* Expand dataset to include 2017-18		*/
/********************************************/

dm log "clear";
dm output "clear";

ods preferences;
ods html close;
ods html ;

libname pride 'C:\Users\brs380\OneDrive - Harvard University\Migrated-P-Drive\NHANES\fped';
libname nhanes 'C:\Users\brs380\OneDrive - Harvard University\Migrated-P-Drive\NHANES\input_dietdata'; 

proc contents data=pride.HEI1112_avg varnum; run; 
proc contents data=pride.HEI1314_avg varnum; run; 
proc contents data=pride.HEI1516_avg varnum; run; 
proc contents data=pride.HEI1718_avg varnum; run; 


data hei1118;
	set pride.hei1112_avg pride.hei1314_avg pride.hei1516_avg pride.hei1718_avg;
run; 

proc sort data=hei1118; by seqn; run; 

proc contents data=pride.tert_nhaneslowF_subset1118 varnum; run; 

data tert1118;
	set pride.tert_nhaneslowF_subset1118; 

drop  _FREQ_;

RUN; 


proc sort data=tert1118; by seqn; run; 

data tert_hei1118_lowF;
	merge hei1118(in=inhei) tert1118(in=intert);
	by SEQN;
	if intert and inhei;
run; 


/*data check*/
proc means data=tert_hei1118clean n nmiss min mean max;
	var RIDAGEYR indfmpir;
run; 


data pride.hei_tert1118_lowF;
	set tert_hei1118_lowF;
run; 

/*export file to csv file*/
