/***********************************************************************************;   
SAS Program:    MakeHEI.sas


Purpose:     SAS program to calculate HEI-2015 scores
                
Data In:    ...\Output\FPED_DR1TOT_XXXX.sas7bdat

Output:     
            

************************************************************************************/
dm log "clear";
dm output "clear";

ods preferences;
ods html close;
ods html ;

/*put demo_g with same files as */

* Create macro var, iPath, with directory of this program.;
*%let iPath=%sysfunc(tranwrd(%sysget(SAS_EXECFILEPATH),\%sysget(SAS_EXECFILENAME),));

%let iPath=C:\Users\brs380\OneDrive - Harvard University\Migrated-P-Drive\NHANES;
* Libname of folder where FPED_DR1TOT_1112 resides.;
libname Input "&iPath\fped"; /* Put all the FPED data here */ 
libname nhanes "&iPath\input_dietdata";


%let yrs=1718;
%let txtsvy=2017-2018;
%let svy=J;

data demo_&svy;
	set input.demo_&svy;
run;

/* 1. Create a folder on your computer “home folder”, and save the FPED data, NHANES data, Demographic data, and the required HEI-2015 macro in it. Specify the path to the folder. */

*%let home = C:\Users\brs380\Documents\FPED_NHANES; /*In this Example, the “home” folder is in C Drive, within Documents, and is called FPED_NHANES. */


/* 2. Libnames here specify the input files. */


*libname NH "&home\NH";  
*libname FPED "&home\FPED"; /*In this Example, the FPED data are in a folder called “FPED”, and the NHANES and Demographic data are in a folder called “NH”, all saved within the “home” folder. These are SAS datasets. */


/* skip -- 3.	Create a folder in the "Home" folder, where the output file, containing HEI-2015 component and total scores for each respondent, for the intake day, are to be exported. Specify the name of the folder. */

/* 4. Read in required HEI-2015 scoring macro. This macro must be saved within the home folder. */

%include "&iPath\hei2015_score_macro.sas";

/*NOTE: Once you have all the steps above, all you need to do is run the SAS program below. Unless you used different names for your datasets and folders, no other action is required from you. */
title 'HEI-2015 scores for NHANES &txtsvy. day 1, AGE >= 20, RELIABLE DIETS, Include Pregnant and Lactating Women';

/*Step 1: locate the required datasets and variables */

*part a: get FPED data per day;
data FPED1;
 set Input.fped_dr1tot_&yrs.;
run;

data FPED2;
	set Input.fped_dr2tot_&yrs.;
run;

*part b: get individual total nutrient intake if reliable recall status;
data NUTRIENT1 (keep=SEQN WTDRD1 DR1TKCAL DR1TSFAT DR1TALCO DR1TSODI DR1DRSTZ DR1TMFAT DR1TPFAT);
  set nhanes.DR1TOT_&svy;
  if DR1DRSTZ=1; /*reliable dietary recall status*/
run;

data NUTRIENT2 (keep=SEQN WTDR2D DR2TKCAL DR2TSFAT DR2TALCO DR2TSODI DR2DRSTZ DR2TMFAT DR2TPFAT);
	SET nhanes.DR2TOT_&svy;
	IF DR2DRSTZ=1;
RUN;

*part c: get demographic data for persons aged 20 and older;
data DEMO (keep=SEQN RIDAGEYR RIAGENDR SDDSRVYR SDMVPSU SDMVSTRA);
  set DEMO_&svy;
  if RIDAGEYR >= 20;
/*
  if      RIDRETH3 = 3            then rac = 1;   * Non-Hisp White;
    else if RIDRETH3 = 4            then rac = 2;   * Non-Hisp Black;
    else if RIDRETH3 = 6            then rac = 3;   * Non-Hisp Asian;
    else if RIDRETH3 in(1 2)        then rac = 4;   * Hispanic      ;
    else                                 rac = 5;   * Other Race    ;

if      0 <= INDFMPIR <= 1.3    then povbin = 1;   * Under 131% pov;
	else if INDFMPIR >1.3 then povbin=2; *above 131% poverty;

	if rac=1 and pov=1 then racpov=11;
		if rac=1 and pov in(2:4) then racpov =12;
		if rac=2 and pov=1 then racpov=21;
		if rac=2 and pov in(2:4) then racpov=22;
		if rac=3 and pov=1 then racpov=31;
		if rac=3 and pov in (2:4) then racpov=32;
		if rac=4 and pov=1 then racpov=41;
		if rac=4 and pov in (2:4) then racpov=42; */
run;


proc sort data=FPED1;
  by SEQN;
run;
proc sort data=FPED2;
  by SEQN;
run;
proc sort data=NUTRIENT1;
  by SEQN;
run;
proc sort data=NUTRIENT2;
  by SEQN;
run;
proc sort data=DEMO;
  by SEQN;
run;

data COHORT1;
  merge NUTRIENT1 (in=N) DEMO (in=D) FPED1;
  by SEQN;
  if N and D;
run;


data COHORT2;
  merge NUTRIENT2 (in=N) DEMO (in=D) FPED2;
  by SEQN;
  if N and D;
run;
/*Step 3: Creates additional required variables: FWHOLEFRT, MONOPOLY, VTOTALLEG, VDRKGRLEG, PFALLPROTLEG and PFSEAPLANTLEG */


data COHORT1;
  set COHORT1;
  by SEQN;

  FWHOLEFRT=DR1T_F_CITMLB+DR1T_F_OTHER;

  MONOPOLY=DR1TMFAT+DR1TPFAT;

  VTOTALLEG=DR1T_V_TOTAL+DR1T_V_LEGUMES;
  VDRKGRLEG=DR1T_V_DRKGR+DR1T_V_LEGUMES;

  PFALLPROTLEG=DR1T_PF_MPS_TOTAL+DR1T_PF_EGGS+DR1T_PF_NUTSDS+DR1T_PF_SOY+DR1T_PF_LEGUMES; 
  PFSEAPLANTLEG=DR1T_PF_SEAFD_HI+DR1T_PF_SEAFD_LOW+DR1T_PF_NUTSDS+DR1T_PF_SOY+DR1T_PF_LEGUMES;
run;


data COHORT2;
  set COHORT2;
  by SEQN;

  FWHOLEFRT=DR2T_F_CITMLB+DR2T_F_OTHER;

  MONOPOLY=DR2TMFAT+DR2TPFAT;

  VTOTALLEG=DR2T_V_TOTAL+DR2T_V_LEGUMES;
  VDRKGRLEG=DR2T_V_DRKGR+DR2T_V_LEGUMES;

  PFALLPROTLEG=DR2T_PF_MPS_TOTAL+DR2T_PF_EGGS+DR2T_PF_NUTSDS+DR2T_PF_SOY+DR2T_PF_LEGUMES; 
  PFSEAPLANTLEG=DR2T_PF_SEAFD_HI+DR2T_PF_SEAFD_LOW+DR2T_PF_NUTSDS+DR2T_PF_SOY+DR2T_PF_LEGUMES;
run;


/*Step 4: Apply the HEI-2015 scoring macro. */

%HEI2015 (indat= COHORT1, 
          kcal= DR1TKCAL, 
	  vtotalleg= VTOTALLEG, 
	  vdrkgrleg= VDRKGRLEG, 
	  f_total= DR1T_F_TOTAL, 
	  fwholefrt= FWHOLEFRT, 
	  g_whole= DR1T_G_WHOLE, 
	  d_total= DR1T_D_TOTAL, 
          pfallprotleg= PFALLPROTLEG, 
	  pfseaplantleg= PFSEAPLANTLEG, 
	  monopoly= MONOPOLY, 
	  satfat= DR1TSFAT, 
	  sodium= DR1TSODI, 
	  g_refined= DR1T_G_REFINED, 
	  add_sugars= DR1T_ADD_SUGARS, 
	  outdat= HEI2015_20xxd1); 


%HEI2015 (indat= COHORT2, 
          kcal= DR2TKCAL, 
	  vtotalleg= VTOTALLEG, 
	  vdrkgrleg= VDRKGRLEG, 
	  f_total= DR2T_F_TOTAL, 
	  fwholefrt= FWHOLEFRT, 
	  g_whole= DR2T_G_WHOLE, 
	  d_total= DR2T_D_TOTAL, 
          pfallprotleg= PFALLPROTLEG, 
	  pfseaplantleg= PFSEAPLANTLEG, 
	  monopoly= MONOPOLY, 
	  satfat= DR2TSFAT, 
	  sodium= DR2TSODI, 
	  g_refined= DR2T_G_REFINED, 
	  add_sugars= DR2T_ADD_SUGARS, 
	  outdat= HEI2015_20xxd2); 

/*Step 5: Displays and saves the results. */ 

*part a: this program saves one HEI-2015 score for each individual, based on one 24HR;


data HEI2015Rd1 (keep=SEQN DRTKCAL  HEI2015C1_TOTALVEG HEI2015C2_GREEN_AND_BEAN HEI2015C3_TOTALFRUIT HEI2015C4_WHOLEFRUIT 
      HEI2015C5_WHOLEGRAIN HEI2015C6_TOTALDAIRY HEI2015C7_TOTPROT HEI2015C8_SEAPLANT_PROT HEI2015C9_FATTYACID HEI2015C10_SODIUM
      HEI2015C11_REFINEDGRAIN HEI2015C12_SFAT HEI2015C13_ADDSUG HEI2015_TOTAL_SCORE); 
  set HEI2015_20xxd1; 
  DRTKCAL=DR1TKCAL;
  run; 
 data HEI2015Rd2 (keep=SEQN DRTKCAL HEI2015C1_TOTALVEG HEI2015C2_GREEN_AND_BEAN HEI2015C3_TOTALFRUIT HEI2015C4_WHOLEFRUIT 
      HEI2015C5_WHOLEGRAIN HEI2015C6_TOTALDAIRY HEI2015C7_TOTPROT HEI2015C8_SEAPLANT_PROT HEI2015C9_FATTYACID HEI2015C10_SODIUM
      HEI2015C11_REFINEDGRAIN HEI2015C12_SFAT HEI2015C13_ADDSUG HEI2015_TOTAL_SCORE); 
  set HEI2015_20xxd2; 
  DRTKCAL=DR2TKCAL;
  run; 
*part b: calculates an unweighted mean across all individuals in group; 
 
proc means n nmiss min max mean data=HEI2015Rd1; 
run; 
 
 
proc means n nmiss min max mean data=HEI2015Rd2; 
run; 
 

/* Combine day 1 and day 2 into a single data file*/

data input.HEI2015_year&yrs.;
	set HEI2015Rd1(in=indayone)
		HEI2015Rd2(in=indaytwo);
if indayone then day=1;
if indaytwo then day=2;

run; 

proc contents data =INPUT.hei2015_year&yrs. varnum;  run; 
proc sort data=INPUT.hei2015_year&yrs.; by seqn; run; 

Proc summary data=INPUT.HEI2015_year&yrs. nway;
    class seqn ;
    var HEI2015C1_TOTALVEG--DRTKCAL;
    output out = input.HEI&yrs._avg mean=;
run;

