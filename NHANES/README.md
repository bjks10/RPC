# RPC APPLICATION TO NHANES DIETARY RECALL DATA
The following procedures and supporting files are needed to successfully run RPC model using 24-hour dietary intake data from the National Health and Nutrition Examination Surveys, pooling four survey cycles: 2011-2012, 2013-2014, 2015-2016, 2017-2018.
Software necessary for implementation:

	* SAS
	* MATLAB
	* R


## STEP 1 - Download FPED data
FPED data can be downloaded directly from the [USDA/ARS](https://www.ars.usda.gov/northeast-area/beltsville-md-bhnrc/beltsville-human-nutrition-research-center/food-surveys-research-group/docs/fped-databases/ "USDA/ARS title") website.
You will need all the food pattern equivalents for foods in the WWEIA, NHANES 20XX-YY. Example attached pools 2011-2018 and can be replicated using the data files located in [nhanes-data](https://github.com/bjks10/NHANES_wtofm/tree/main/nhanes-data "nhanes-data title") subfolder. 

## STEP 2 - Make Food Groups from 24HR Diet Recall data
Run MakeFoodGroups_nhanes.sas  - run for each survey cycle separately. Make sure the order of variables match with the cycle run. Example runs 2011-12 (G), 2013-14 (H), 2015-16 (I), 2017-18 (J)can be directly found in [nhanes-data](https://github.com/bjks10/NHANES_wtofm/tree/main/nhanes-data "nhanes-data title") subfolder.
Note: The ordering of food group variables changed in 2017-2018 (J) cycle
	
	*Input files:	FPED_DR1TOT_XXYY.sas7bdat
			FPED_DR2TOT_XXYY.sas7bdat
			FPED_DR1IFF_XXYY.sas7bdat
			FPED_DR2IFF_XXYY.sas7bdat
	*Output files:	FPED_DRAVGXXYY.sas7bdat



## STEP 3 - Create Tertiles of Consumption
Run fped_CreateTerts.sas to combine survey cycles to a pooled dataset, drop "total" food labels, add diet weights. 
XXYY denote a single survey cycle (e.g. 1112 for 2011-2012 cycle). 
XXZZ is the full range of combined surveys (e.g. 1118 for 2011-2018 surveys pooled). 
[A] denotes the alpha code attached to each survey cycle: G=2011-2012, H=2013-2014, I=2015-2016, J=2017-2018.
	
	 *Input files: 	DEMO_[A].sas7bdat
			FPED_DRAVGXXYY.sas7bdat
			DR1TOT_[A].sas7bdat
	*Output files:	tert_nhanesXXZZ.sas7bdat


## STEP 4 - Calculate HEI2015 scores for all participants
Import SAS data files into SAS and run GenerateHEI.sas to calculate HEI2015 score for each survey cycle [X].
  	
	* Input files:	DEMO_X.sas7bdat
			FPED_DR1TOT_X.sas7bdat
			FPED_DR2TOT_X.sas7bdat
			hei2015_score_macro.sas
	* Output files:	HEIXXYY_avg.sas7bdat

## STEP 5 - Merge HEIscores to FPED tertile data
Import sas7bdat files into SAS and run Merge_FPEDHEI.sas to merge FPED-tertiles and HEI2015 score and demographic data to single dataset and export to CSV file.

  	* Input files:	HEIXXYY_avg.sas7bdat
			tert_nhanesXXZZ.sas7bdat
	* Output files: hei_tert1118.sas7bdat
			adultHEI_fped1118terts.csv

## STEP 6 - Generate target study population (Low-income adult women)
Import sas7bdat file into SAS and run Create_LowFdata.sas to subset population to female adult participants living at or below the 130% poverty income level.

	* Input files:	hei_tert1118.sas7bdat
	* Output files:	nhanes1118_lowFdietdata.sas7bdat
			nhanes_lowFadultdata1Dec2021.csv

## STEP 7 - Run Robust Profile Clustering (RPC) model
Import CSV file into MATLAB and run RPC_nhanesLOWF.m for global and local diet patterns of low-income adult women, where local patterns are specific to a participant's racial/ethnic background. 

	* Input files:		nhanes_lowFadultdata1Dec2021.csv
	* Output files:		RPCnhanesLOWF_MCMCmed.mat
				NHANESLowFrpc_assign1Dec2021.csv
	* Output figures: 	nhanesRPC_dendrogram.png
				noconsum_lowadultpat.fig
				highconsum_lowadultpat.fig
				theta0_lowFpattern.fig
				nu_lowwomen.fig

## STEP 8 - nhanesRPC_tables.R
Import saved output files into R. 

	* Input files:		NHANESLowFrpc_assign1Dec2021.csv
				nhanesadult_cvdriskHEIFHR1118.RData - created in nhanes_adultCVDHEI.R
				nhanes_lowFadultdata1Dec2021.csv 
	* Output files: 	NHANESLowF_indata_rpc1Dec2021.csv
