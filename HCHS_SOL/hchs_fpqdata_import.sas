/************************************************************/
/* Hispanic Community Health Survey/Study of Latinos		*/
/* Food Propensity Questionnaire: EXPLORATORY ANALYSIS		*/
/* PROGRAMMER: Briana Stephenson							*/
/* PI: Jianwen Cai											*/
/* Senior Personnel: Daniela Sotres-Alvarez					*/
/* 															*/
/* FPQ data - create consumption levvel vars				*/
/************************************************************/

dm log 'clear';
ods html close; 
ods preferences;
ods html newfile=none;

%let path = P:\HCHS_SOL\Data;
*%let home = P:\Stephenson_Analysis;
%let PW = LOS_SHCH;




%include "&path\hchsFPE_MACROVARS.sas";

*libname bjks "&home";
libname HCHS "&path";

data subp_inv4;
	set hchs.part_derv_inv4(pw=&pw);

	*create subpopulations;
	if center='B' and bkgrd1_c7=0 then subpop=1; *bronx dominican;
	if center='B' and bkgrd1_c7=4 then subpop=2; *bronx puerto rican;
	if center='C' and bkgrd1_c7=1 then subpop=3; *chicago central american;
	if center='C' and bkgrd1_c7=3 then subpop=4; *chicago mexican;
	if center='C' and bkgrd1_c7=4 then subpop=5; *chicago puerto rican;
	if center='M' and bkgrd1_c7=1 then subpop=6; *miami central american;
	if center='M' and bkgrd1_c7=2 then subpop=7; *miami cuban;
	if center='M' and bkgrd1_c7=5 then subpop=8; *miami south american;
	if center='S' and bkgrd1_c7=3 then subpop=9; *san diego mexicans;
	label subpop='Site-Ethnicity';
	if cmiss(of subpop) then delete;
keep id subpop;
run; 

proc sort data=subp_inv4; by id; run;
proc freq data=hchs.fpe_inv4(pw=&pw);
	tables fpe81f fpe81ap fpe80f fpe80ap fpe80bp fpe80cp fpe8f fpe8ap fpe11f fpe11ap fpe12f fpe12ap fpe38f fpe38ap fpe46f fpe46ap fpe53f fpe53ap fpe56f fpe56ap fpe57f fpe57ap fpe99f fpe99ap fpe114f fpe115a fpe115b fpe115c fpe115d;
run;

data fpe_inv4;
	set hchs.fpe_inv4(pw=&pw);
	keep id fpe115a--fpe115d;
	*fruit drink;
fpe5F=fpe5f*(1-fpe5ap); *regular fruit drink;
fpe5ap=fpe5f*fpe5ap;	*diet fruit drink;
	*soda;
fpe8f=fpe8f*(1-fpe8ap); *regular soda;
fpe8ap=fpe8f*fpe8ap; *diet soda;
	*cooked cereal;
fpe11f=fpe11f*fpe11ap; *oatmeal;
fpe11ap=fpe11f*(1-fpe11ap); *other cooked cereal;
	*cold cereal;
fpe12f=fpe12f*(1-fpe12ap);
fpe12ap=fpe12f*fpe12ap;
	*lettuce salad;
fpe38f=fpe38f*(1-fpe38ap); *other lettuce salad;
fpe38ap=fpe38f*fpe38ap; *dark green leaf salad;
	*tortillas tacos;
fpe46f=fpe46f*(1-fpe46ap); *other tortillas/tacos;
fpe46ap=fpe46f*fpe46ap; *corn tortillas/tacos;
	*rice;
fpe53f=fpe53f*(1-fpe53ap); *other rice;
fpe53ap=fpe53f*fpe53ap; *whole grain rice;
	*BREAD;
	*sandwich bread;
fpe56f=fpe56f*(1-fpe56ap); *non-white bread;
fpe56ap=fpe56f*fpe56ap; *white sandwich bread;
	*breads/rolls;
fpe57f=fpe57f*(1-fpe57ap); *nonwhite bread;
fpe57ap=fpe57f*fpe57ap; *white bread;
	*soup;

fpe80ap=fpe80f*fpe80ap; *bean soup;
fpe80bp=fpe80f*fpe80bp; *cream soup;
fpe80cp=fpe80f*fpe80cp; *tomato/veg soup;
fpe80f=fpe80f-(fpe80ap+fpe80bp+fpe80cp); *all other soup;
	if fpe80f<0 then fpe80f=0;
	*pizza;
fpe81f=fpe81f*(1-fpe81ap);
fpe81ap=fpe81f*fpe81ap;
	*pie;
fpe99f=fpe99f*(1-fpe99ap); *other pies;
fpe99ap=fpe99f*fpe99ap; *fruit pie;
run;
proc contents data=fpe_inv4 varnum; run;
proc freq data=fpe_inv4;
	tables fpe81f fpe81ap fpe80f fpe80ap fpe80bp fpe80cp fpe8f fpe8ap fpe11f fpe11ap fpe12f fpe12ap fpe38f fpe38ap fpe46f fpe46ap fpe53f fpe53ap fpe56f fpe56ap fpe57f fpe57ap fpe99f fpe99ap fpe114f fpe115a fpe115b fpe115c fpe115d;
run;

data fpe_combo;
	set fpe_inv4;
fpe57ap=fpe57ap+fpe56ap; *combine white breads;

fpe56f=fpe56f+fpe57f; *combine whole-grain breads;

fpe59F=fpe59f+fpe68f; *roast beef;

run; 

%let totfoods=132;
/*Reorder foodvariables in numeric order */

data fpe_reord;
	set fpe_combo;
	array fpe(*) _NUMERIC_;
	array fpq(*) fpq1-fpq&totfoods;
	do i = 1 to &totfoods;	
		fpq{i}=fpe{i};
	end;
	drop fpe115A--fpe115d;
	label &fpq;
run;

proc contents data=fpe_reord varnum; run; 


/********************************/
/*Create 5 consumption levels	*/
/*(1)no consumption 			*/
/*(2)<= 2-3/month				*/
/*(3)3-4/week					*/
/*(4)almost daily	 			*/
/*(5)daily+		 				*/
/********************************/

data fp_pent;
	set fpe_reord;
	array fpq(*) fpq1-fpq&totfoods;
	array fp_bin(*) fp_bin1-fp_bin&totfoods;
	do i = 1 to &totfoods;	
		if fpq[i]=0 then fp_bin[i] = 1; *no consumption;
		else if (fpq[i]>0 and fpq[i] <= (2.5/(365.25/12))) then fp_bin[i]=2; *2-3 months or less;
		else if (fpq[i]>(2.5/(365.25/12)) and fpq[i] <=(3.5/7)) then fp_bin[i]=3; *3-4 week - 2.5/month;
		else if (fpq[i]>(3.5/7) and fpq[i] <=1) then fp_bin[i]=4; *almost daily;
		else if fpq[i]>1  then fp_bin[i]=5; *more than daily;
   	end;
drop fpq1-fpq&totfoods fp_bin54 fp_bin65 fp_bin109; *drop duplicate of white bread, whole grain, roast beef;
label &fp_bin;
run;

proc contents data=fp_pent varnum; run; 
proc freq data=fp_pent;
	tables fp_bin1--fp_bin&totfoods;
run;
/******************************/
/*collapse FPQ foods as needed*/
/******************************/

data fp_bi;
	set fp_pent;
/*binary foods*/
array fp_bin[15] fp_bin1 fp_bin8 fp_bin14-fp_bin16 fp_bin19 fp_bin94 fp_bin102-fp_bin104  fp_bin114 fp_bin115 fp_bin130-fp_bin132;
array fp[15] fp1 fp8 fp14-fp16 fp19 fp94 fp102-fp104  fp114 fp115 fp130-fp132;
do i = 1 to 15;
	if fp_bin[i]<=2 then fp[i]=fp_bin[i];
	else if fp_bin[i]>2 then fp[i]=2;
end;

drop i;

run; 

data fp_tri;
	set fp_bi;
/*tertiary foods*/

array fp_bin[51] fp_bin9 fp_bin18 fp_bin33  fp_bin37-fp_bin40 fp_bin43 fp_bin44 fp_bin47 fp_bin51 fp_bin59-fp_bin64
			 		fp_bin66-fp_bin68 fp_bin71-fp_bin76 fp_bin78-fp_bin81 fp_bin84 fp_bin86
					fp_bin88-fp_bin90 fp_bin92 fp_bin93 fp_bin99
					fp_bin100 fp_bin105-fp_bin108 fp_bin110-fp_bin112 fp_bin117 fp_bin123
					fp_bin127-fp_bin129;
array fp[51] fp9 fp18 fp33 fp37-fp40 fp43 fp44 fp47 fp51 fp59-fp64 fp66-fp68 fp71-fp76 fp78-fp81 fp84 fp86 fp88-fp90 fp92 fp93 
					fp99 fp100 fp105-fp108 fp110-fp112 fp117 fp123 fp127-fp129;
do i=1 to 51;
	if fp_bin[i]<=3 then fp[i]=fp_bin[i];
		else if fp_bin[i]>3 then fp[i]=3;
end;
drop i;
run;

data fp_qua;
	set fp_tri;
/*quartenary foods*/
array fp_bin[53] fp_bin3-fp_bin6 fp_bin20-fp_bin32 fp_bin34-fp_bin36 fp_bin41 fp_bin42 fp_bin45 fp_bin46
					fp_bin48-fp_bin50 fp_bin52 fp_bin55-fp_bin58 fp_bin69 fp_bin70 fp_bin77 fp_bin82 fp_bin83 fp_bin85 fp_bin87 
					fp_bin91 fp_bin95-fp_bin98 fp_bin113 fp_bin118-fp_bin122 fp_bin124-fp_bin126 ;
array fp[53] fp3-fp6 fp20-fp32 fp34-fp36 fp41 fp42 fp45 fp46 fp48-fp50 fp52 fp55-fp58 fp69 fp70 fp77 fp82 fp83 fp85 fp87
				fp91 fp95-fp98 fp113 fp118-fp122 fp124-fp126;
do i = 1 to 53;
	if fp_bin[i]<=4 then fp[i]=fp_bin[i];
	else if fp_bin[i]>4 then fp[i]=4;
end; 
drop i;
run;

data fpe_bins;
	set fp_qua;
array fp_bin[10] fp_bin2 fp_bin7 fp_bin10-fp_bin13 fp_bin17 fp_bin53 fp_bin101 fp_bin116;
array fp[10] fp2 fp7 fp10-fp13 fp17 fp53 fp101 fp116;
do i=1 to 10;
	fp[i]=fp_bin[i];
end;
drop i; 
keep id fp1-fp&totfoods;
label &fp;
run;
proc contents data=fpe_bins varnum; run; 

proc freq data=fpe_bins;
	tables fp1-fp53 fp55-fp64 fp66-fp108 fp110-fp&totfoods;
run; 
/*reorder variables by similar food types*/

data fpq_bins;
	retain id fp2-fp6 fp102 fp103 fp116 fp9  fp117 fp7  fp8 fp10-fp17 fp19-fp28 fp119-fp127 fp124-fp127 fp29-fp35 fp106 fp36-fp43
			fp45-fp50 fp108 fp51-fp53 fp110 fp118 fp104 fp105 fp18 fp107 fp44 fp55-fp64 fp66-fp75 fp128 fp76 fp114 fp77-fp94
			fp115 fp95-fp101 fp111-fp113 fp129-fp132 fp1;
	set fpe_bins;
	

run;

proc contents data=fpq_bins varnum; run; 

proc sort data=fpq_bins; by id ; run;
/*combine subpop id to fpq data*/

data fpq_subpop;
	merge fpq_bins(in=a) subp_inv4;
	by id;
	if a;

	if cmiss(of _all_) then delete;
	
run;
/*save permanent dataset*/

data hchs.fpq_bins;
	set fpq_subpop;
run; 

proc export 
  data=hchs.fpq_bins 
  dbms=xlsx 
  outfile="c:\temp\hchs_fpqdata_export.xlsx" 
  replace;
run;
