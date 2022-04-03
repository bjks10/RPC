/***********************************************************************************;   
SAS Program:    MakeFoodGroups_nhanes.sas


Purpose:    Create NHANES survey cycle USDA 20**-20++ Food 
            Patterns Equivalents and create tables
                
Data In:    ...\FPED_DR1TOT_xxyy.sas7bdat
Data In:    ...\FPED_DR2TOT_xxyy.sas7bdat


Output:     ...\fped_dravgxxyy.sas7bdat       
Food Patterns Equivalents Database per 100 grams of FNDDS 20xx-20yy foods

************************************************************************************/


libname fped 'C:\Users\brs380\OneDrive - Harvard University\Migrated-P-Drive\NHANES\fped';
%let yrs=1718;

data adultdata1;
    set fped.FPED_DR1TOT_&yrs. (where = ((RIDAGEYR >= 20) and
                                (DRABF ne 1) and DR1DRSTZ = 1));

    all = 1;
    sex = RIAGENDR;

  
    if 20 <= RIDAGEYR <= 29    then ag1 = 4;
    else if 30 <= RIDAGEYR <= 39    then ag1 = 5;
    else if 40 <= RIDAGEYR <= 49    then ag1 = 6;
    else if 50 <= RIDAGEYR <= 59    then ag1 = 7;
    else if 60 <= RIDAGEYR <= 69    then ag1 = 8;
    else if RIDAGEYR >= 70          then ag1 = 9;

    if      RIDRETH3 = 3            then rac = 1;   * Non-Hisp White;
    else if RIDRETH3 = 4            then rac = 2;   * Non-Hisp Black;
    else if RIDRETH3 = 6            then rac = 3;   * Non-Hisp Asian;
    else if RIDRETH3 in(1 2)        then rac = 4;   * Hispanic      ;
    else                                 rac = 5;   * Other Race    ;

    if      INDFMIN2 in(1:5 13)     then inc = 1;   * Under $20k    ;
    else if INDFMIN2 in(6:10)       then inc = 2;   * $20 - $75k    ;
    else if INDFMIN2 in(14 15)      then inc = 3;   * $75k and Over ;
    else                                 inc = 4;   * Other         ;

    if      0 <= INDFMPIR <= 1.3    then pov = 1;   * Under 131% pov;
    else if 1.3 < INDFMPIR <= 3.50  then pov = 2;   * 131-350% pov  ;
    else if INDFMPIR > 3.50         then pov = 3;   * Over 350% pov ;
    else                                 pov = 4;   * Other/missing         ;
    
	if      0 <= INDFMPIR <= 1.3    then povbin = 1;   * Under 131% pov;
    else if INDFMPIR > 1.3         then povbin = 2;   * Over 131% pov ;
    else                                 povbin = 3;   * Other/Missing         ;
   * Apply shortened labels, these will appear in table headings.;
    label DR1T_G_TOTAL          = "Total Grain";
    label DR1T_G_WHOLE          = "Whole Grains";
    label DR1T_G_REFINED        = "Refined Grains";
    label DR1T_V_TOTAL          = "Total Vegetables";
    label DR1T_V_DRKGR          = "Dark Green";
    label DR1T_V_REDOR_OTHER    = "Other Red Orange";
    label DR1T_V_STARCHY_TOTAL  = "Total Starchy";
    label DR1T_V_STARCHY_POTATO = "Potatoes";
    label DR1T_V_STARCHY_OTHER  = "Other Starchy";
    label DR1T_V_REDOR_TOTAL    = "Total Red and Orange";
    label DR1T_V_REDOR_TOMATO   = "Tomatoes";
    label DR1T_V_OTHER          = "Other";
    label DR1T_F_TOTAL          = "Total Fruit";
    label DR1T_F_CITMLB         = "Citrus, Melons and Berries";
    label DR1T_F_OTHER          = "Other Fruit";
    label DR1T_F_JUICE          = "Fruit Juice";
    label DR1T_D_TOTAL          = "Total Dairy";
    label DR1T_D_MILK           = "Fluid Milk";
    label DR1T_D_YOGURT         = "Yogurt";
    label DR1T_D_CHEESE         = "Cheese";
    label DR1T_PF_TOTAL         = "Total Protein";
    label DR1T_PF_MPS_TOTAL     = "Total Meat, Poultry, and Seafood";
    label DR1T_PF_MEAT          = "Meat";
    label DR1T_PF_ORGAN         = "Organ";
    label DR1T_PF_CUREDMEAT     = "Cured Meats";
    label DR1T_PF_POULT         = "Poultry";
    label DR1T_PF_SEAFD_HI      = "Seafood High n-3";
    label DR1T_PF_SEAFD_LOW     = "Seafood Low n-3";
    label DR1T_PF_EGGS          = "Eggs";
    label DR1T_PF_SOY           = "Soybean Products";
    label DR1T_PF_NUTSDS        = "Nuts and Seeds";
    label DR1T_V_LEGUMES        = "Legumes as Vegetable";
    label DR1T_PF_LEGUMES       = "Legumes as Protein";
    label DR1T_OILS             = "Oils";
    label DR1T_SOLID_FATS       = "Solid Fats";
    label DR1T_ADD_SUGARS       = "Added Sugar";
    label DR1T_A_DRINKS         = "Alcoholic Drinks";
run;


data adultdata2;
    set fped.FPED_DR2TOT_&yrs. (where = ((RIDAGEYR >= 20) and
                                (DRABF ne 1) and DR2DRSTZ = 1));

    all = 1;
    sex = RIAGENDR;

  
    if 20 <= RIDAGEYR <= 29    then ag1 = 4;
    else if 30 <= RIDAGEYR <= 39    then ag1 = 5;
    else if 40 <= RIDAGEYR <= 49    then ag1 = 6;
    else if 50 <= RIDAGEYR <= 59    then ag1 = 7;
    else if 60 <= RIDAGEYR <= 69    then ag1 = 8;
    else if RIDAGEYR >= 70          then ag1 = 9;

    if      RIDRETH3 = 3            then rac = 1;   * Non-Hisp White;
    else if RIDRETH3 = 4            then rac = 2;   * Non-Hisp Black;
    else if RIDRETH3 = 6            then rac = 3;   * Non-Hisp Asian;
    else if RIDRETH3 in(1 2)        then rac = 4;   * Hispanic      ;
    else                                 rac = 5;   * Other Race    ;

    if      INDFMIN2 in(1:5 13)     then inc = 1;   * Under $20k    ;
    else if INDFMIN2 in(6:10)       then inc = 2;   * $20 - $75k    ;
    else if INDFMIN2 in(14 15)      then inc = 3;   * $75k and Over ;
    else                                 inc = 4;   * Other         ;

    if      0 <= INDFMPIR <= 1.3    then pov = 1;   * Under 131% pov;
    else if 1.3 < INDFMPIR <= 3.50  then pov = 2;   * 131-350% pov  ;
    else if INDFMPIR > 3.50         then pov = 3;   * Over 350% pov ;
    else                                 pov = 4;   * Other/missing         ;
    
	if      0 <= INDFMPIR <= 1.3    then povbin = 1;   * Under 131% pov;
    else if INDFMPIR > 1.3         then povbin = 2;   * Over 131% pov ;
    else                                 povbin = 3;   * Other/Missing         ;

   run;



********************************************************************;
*                                                                   ;
*  Create formats for the group variables.                          ;
*                                                                   ;
********************************************************************;
proc format;
    value agef
        4   = "20 - 29.............."
        5   = "30 - 39.............."
        6   = "40 - 49.............."
        7   = "50 - 59.............."
        8   = "60 - 69.............."
        9   = "    70 and over......";

    value sexf 
        1   = "Males:"
        2   = "Females:"
        3   = "Males and females:";

    value racf
        1   = "Non-Hispanic White:"
        2   = "Non-Hispanic Black:"
        3   = "Non-Hispanic Asian:"
        4   = "Hispanic/Latino:";

    value incf
        1   = "$0 - $24,999:"
        2   = "$25,000 - $74,999:"
        3   = "$75,000 and higher:"
        4   = "All Individuals:";

    value povf
        1   = "Under 131% poverty:"
        2   = "131-350% poverty:"
        3   = "Over 350% poverty:"
        4   = "All Individuals:";
	value racpovf
		11 = "NH White-Below POV"
		12 = "NH White-Above POV"
		21 = "NH Black-Below POV"
		22 = "NH Black-Above POV"
		31 = "NH Asian-Below POV"
		32 = "NH Asian-Above POV"
		41 = "Hispanic/Latino-Below POV"
		42 = "Non-Hispanic/Latino-Above POV";
	value povbin
		1	= "Below 131% poverty"
		2	= "Above 131% poverty";
run;
proc contents data=adultdata1 varnum;run;
data day1;
	set adultdata1;

array drf {37} DR1T_F_TOTAL--DR1T_A_DRINKS; *CHECK THE ORDER OF FOODS FROM CONTENTS;
array foodgroup {37} foodgroup1-foodgroup37;
	do i=1 to 37;
		foodgroup{i}=drf{i};
	end;

	drop i DR1T_F_TOTAL--DR1T_A_DRINKS;
run; 

data day2;
	set adultdata2;

array drf {37} DR2T_F_TOTAL--DR2T_A_DRINKS; *CHECK THE ORDER OF FOODS FROM CONTENTS;
array foodgroup {37} foodgroup1-foodgroup37;
	do i=1 to 37;
		foodgroup{i}=drf{i};
	end;
	drop i DR2T_F_TOTAL--DR2T_A_DRINKS;
run; 

data adultdata_all;
	set day1 (in=indayone)
		day2 (in=indaytwo);
if indayone=1 then day=1;
if indaytwo=1 then day=2;
	 * Apply shortened labels, these will appear in table headings.;
    label foodgroup1          = "Total Fruit";
	label foodgroup2          = "Citrus, Melons and Berries";
    label foodgroup3          = "Other Fruit";
    label foodgroup4        = "Fruit Juice";
	label foodgroup5          = "Total Vegetables";
    label foodgroup6          = "Dark Green Vegetables";
	label foodgroup7 = "Total Red/Orange Vegetables";
    label foodgroup8    = "Tomatoes";
    label foodgroup9  = "Other Red/Orange Vegetables";
    label foodgroup10   = "Total Starchy Vegetables";
    label foodgroup11  = "Potatoes";
    label foodgroup12    = "Other Starchy Vegetables";
    label foodgroup13          = "Other Vegetables";
    label foodgroup14         = "Legumes (vegetables)";
	label foodgroup15          = "Total Grain";
    label foodgroup16          = "Whole Grain";
    label foodgroup17          = "Refined Grain";
    label foodgroup18          = "Total protein";
	label foodgroup19         = "Total Meat, Poultry, Seafood";
    label foodgroup20           = "Meat, ns";
    label foodgroup21         = "Cured meats";
    label foodgroup22         = "Organ meat";
    label foodgroup23         = "Poultry";
    label foodgroup24     = "Seafood, high n-3";
    label foodgroup25          = "Seafood, low n-3";
    label foodgroup26     = "Eggs";
    label foodgroup27         = "Soybean products";
    label foodgroup28      = "Nuts and seeds";
    label foodgroup29     = "Legumes (protein)";
    label foodgroup30       = "Total Dairy";
    label foodgroup31           = "Fluid milk";
    label foodgroup32        = "Yogurt";
    label foodgroup33        = "Cheese";
    label foodgroup34             = "Oils";
    label foodgroup35       = "Solid Fats";
    label foodgroup36       = "Added Sugar";
    label foodgroup37         = "Alcoholic Drinks";

run;
proc contents data=adultdata_all varnum; run; 

proc sort data=adultdata_all; by seqn; run; 
/*take average of two recalls*/
Proc summary data=adultdata_all nway;
    class seqn ;
    var foodgroup1--foodgroup37;
    output out = fg_avg mean=;
run;

proc sort data = adultdata_all; by seqn; run; 


data fped.fped_dravg&yrs.;
	set fg_avg;
	nrecall=_FREQ_;
	drop _TYPE_;
run; 


proc contents data=fped.fped_dravg1516 varnum; run; 
proc contents data=fped.fped_dravg1718 varnum; run; 
