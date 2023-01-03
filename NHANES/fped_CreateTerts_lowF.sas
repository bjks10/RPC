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


libname fped 'C:\Users\brs380\OneDrive - Harvard University\Migrated-P-Drive\NHANES\fped';
libname nhanes 'C:\Users\brs380\OneDrive - Harvard University\Migrated-P-Drive\NHANES\input_dietdata'; 

/*manually import demographic data from NHANES 2011-18: demo_g, demo_h, demo_i, demo_j*/

data fped.demos1118;
	set demo_g(in=a) demo_h(in=b) demo_i(in=c) demo_j(in=d);
	wtint8yr=wtint2yr/4;
	wtmec8yr=wtmec2yr/4;
if a then cycle='2011-2012';
if b then cycle='2013-2014';
if c then cycle='2015-2016';
if d then cycle='2017-2018';
keep SEQN cycle sdmvpsu sdmvstra wtint8yr wtmec8yr wtint2yr wtmec2yr riagendr ridageyr dmdmartl dmdeduc2 indfmpir ridreth3; 

run; 

proc contents data=fped.demos1118; run;

data demo1118;
	SET fped.demos1118;
	if ridageyr>=20;
run; 
data fped_dravg1116;
	set fped.fped_dravg1112 fped.fped_dravg1314  fped.fped_dravg1516;
array fg {37} fg1-fg37;
array foodgroup{37} foodgroup1--foodgroup37;

do i=1 to 37;
	fg{i}=foodgroup{i};
end;
drop i foodgroup1-foodgroup37;
run; 

proc contents data=fped.fped_dravg1112; run; 
proc contents data=fped.fped_dravg1314; run;
proc contents data=fped.fped_dravg1516; run; 
proc contents data=fped.fped_dravg1718; run; 

data fped_dravg1718;
*reorder variables to match 11-16;
retain seqn _Freq_ foodgroup2-foodgroup4 foodgroup1 foodgroup6 foodgroup8 foodgroup9 foodgroup7
foodgroup11 foodgroup12 foodgroup10 foodgroup13 foodgroup5 foodgroup14 foodgroup16 foodgroup17-
foodgroup15 foodgroup20-foodgroup25 foodgroup19 foodgroup26-foodgroup29 foodgroup18 foodgroup30-foodgroup37 nrecall;
set pride.fped_dravg1718;

array fg{37} fg1-fg37;
array foodgroup{37} foodgroup2--foodgroup37;

do i=1 to 37;
	fg{i} = foodgroup{i};
end;
drop i foodgroup2--foodgroup37; 
run;
proc contents data=fped_dravg1718 varnum; run; 
 
data fped1118;
	set fped_dravg1116 fped_dravg1718;
	drop fg4 fg8 fg11 fg13 fg17 fg24 fg29 fg33; 
run;

proc contents data=fped1118; run;
proc sort data=fped1118; by seqn; run;
proc sort data=demo1118; by seqn; run;  

/* Pull dietary weights from NHANES website */
data dietwt_ghij;
	set nhanes.dr1tot_g(keep=seqn wtdrd1 wtdr2d) 
		nhanes.dr1tot_h(keep=seqn wtdrd1 wtdr2d)
		nhanes.dr1tot_i(keep=seqn wtdrd1 wtdr2d)
		nhanes.dr1tot_j(keep=seqn wtdrd1 wtdr2d);

run; 
proc sort data=dietwt_ghij; by seqn; run; 

data fped_nhanesadult;
	merge demo1118 fped1118(in=fped) dietwt_ghij;
by seqn;

if ridageyr >=20; *keep only adults over 20;
*if ridreth3 in (1,2,3,4,6); *drop mixed race;

if nrecall=1 then dietwt=wtdrd1;
	else if nrecall=2 then dietwt=wtdr2d;
dietwt8yr = dietwt/4; 
if fped;
run; 
proc contents data=fped_nhanesadult ; run;


data fped_adult;
	set fped_nhanesadult;
	array fg {29} fg1--fg37;
	array bb {29} bb1-bb29;
 do i=1 to 29;
 	if fg{i}=0 then bb{i}=.;
		else bb{i}=fg{i};
 end; 
linkvar=1;
run; 

/*IDENTIFY PERSONS WHO ARE PREGNANT OR BREASTFEEDING */

proc contents data=fped.rhq_g varnum; run; 

data fped.rhq1118;
	set fped.rhq_g fped.rhq_h fped.rhq_i fped.rhq_j;

	if rhd143=1 then preg=1;
		else preg=0;
	if rhq200=1 then bfeed=1;
		else bfeed=0;
	if preg=1 or bfeed=1 then bfpreg=1;
		else bfpreg=0;

keep seqn bfeed preg bfpreg;
run; 

proc freq data=fped.rhq1118;
	tables bfeed preg bfpreg;
run;

proc sort data=fped.rhq1118; by seqn; run; 

proc sort data= fped_adult; by seqn; run;
/*Exclusion criteria 1: Female Adults */
data fped_rhq_ex1;
	merge fped.rhq1118 fped_adult(in=a); 
	by seqn; 
	if a;
	if riagendr=2;
run; 
/*Exclusion criteria 2: Living at or below 130% poverty income level */

data fped_rhq_ex2;
	set fped_rhq_ex1; 
	if indfmpir = . then delete;
	if indfmpir <=1.3 ; 
run;
/*examine breastfeed/lactating mom consumption of legumes */ 
proc sort data=fped_rhq_ex2; by bfpreg; run; 
proc surveymeans data=fped_rhq_ex2 percentile=(33,66) nobs nmiss plots=none;
	strata sdmvstra;
	cluster sdmvpsu;
	weight dietwt8yr;
	var bb10 bb22;
	by bfpreg;
	*ods output Quantiles=nhanesquant_lowF;
run; 
/*Exclusion criteria 3: Remove women who are pregnant or breastfeeding */
data fped_rhq_ex3;
	set fped_rhq_ex2;
if bfpreg=1 or bfpreg=. then delete;
run;

/*Calculate tertiles based on those not lactating or breastfeeding dataset */ 

proc surveymeans data=fped_rhq_ex3 percentile=(33,66) nobs nmiss plots=none;
	strata sdmvstra;
	cluster sdmvpsu;
	weight dietwt8yr;
	var bb1-bb29;
	ods output Quantiles=nhanesquant_lowF;
run; 


/*

proc surveymeans data=fped_adult percentile=(33,66) nobs nmiss plots=none;
	strata sdmvstra;
	cluster sdmvpsu;
	weight dietwt8yr;
	var bb1-bb29;
	ods output Quantiles=nhanesquant;
run; 
*/

data wtprctl33;
	set nhanesquant_lowF;
if quantile ^= 0.33 then delete;
keep varname estimate;
run; 

data wtprctl66;
	set nhanesquant_lowF;
if quantile ^= 0.66 then delete;
keep varname estimate;
run; 


proc transpose data=wtprctl33 out=widepct33 prefix=pctl33_;
    id varname;
    var estimate;
run;


proc transpose data=wtprctl66 out=widepct66 prefix=pctl66_;
    id varname;
    var estimate;
run;


data widepct33;
	set widepct33;
	linkvar=1;
run; 

data widepct66;
	set widepct66;
	linkvar=1;
run; 

data tert_nhaneslowf;
	merge fped_rhq_ex3 widepct33 widepct66;
	by linkvar;
	array fg {29} fg1--fg37;
	array pctl33_bb {29} pctl33_bb1-pctl33_bb29;
	array pctl66_bb {29} pctl66_bb1-pctl66_bb29;
	array bb {29} bb1-bb29;
 do i=1 to 29;
 	if fg{i}=0 then bb{i}=1;
		else if fg{i} >0 and fg{i} <= pctl33_bb{i} then bb{i}=2;
			else if fg{i}> pctl33_bb{i} and fg{i} <= pctl66_bb{i} then  bb{i}=3;
				else if fg{i} > pctl66_bb{i} then bb{i}=4;
end;
drop pctl33_bb1-pctl33_bb29 pctl66_bb1-pctl66_bb29 fg1--fg37 i linkvar  _NAME_;

run; 
proc contents data=tert_nhaneslowf varnum; run; 
/*
data tert_nhanes;
	merge fped_adult widepct33 widepct66;
	by linkvar;
	array fg {29} fg1--fg37;
	array pctl33_bb {29} pctl33_bb1-pctl33_bb29;
	array pctl66_bb {29} pctl66_bb1-pctl66_bb29;
	array bb {29} bb1-bb29;
 do i=1 to 29;
 	if fg{i}=0 then bb{i}=1;
		else if fg{i} >0 and fg{i} <= pctl33_bb{i} then bb{i}=2;
			else if fg{i}> pctl33_bb{i} and fg{i} <= pctl66_bb{i} then  bb{i}=3;
				else if fg{i} > pctl66_bb{i} then bb{i}=4;

end;
drop pctl33_bb1-pctl33_bb29 pctl66_bb1-pctl66_bb29 fg1--fg37 i linkvar  _NAME_;

run; 
*/
proc sort data=tert_nhanesLOWF; by seqn; run; 
proc contents data=tert_nhanesLOWF varnum; run; 

proc freq data=tert_nhanesLOWF;
	tables cycle* ridreth3/nopercent;
run; 

proc surveyfreq data=tert_nhanesLOWF;
	cluster sdmvpsu;
	strata sdmvstra;
	weight dietwt8yr;
	tables bb1-bb29;
run; 

data fped.tert_nhaneslowF_subset1118;
	set tert_nhanesLOWF;
run; 
