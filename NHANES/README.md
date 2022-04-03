# RPC APPLICATION TO NHANES DIETARY RECALL DATA
The following procedures and supporting files are needed to successfully run RPC model using 24-hour dietary intake data from the National Health and Nutrition Examination Surveys, pooling four survey cycles: 2011-2012, 2013-2014, 2015-2016, 2017-2018.
Software necessary for implementation


## STEP 1 - Download FPED data from: https://www.ars.usda.gov/northeast-area/beltsville-md-bhnrc/beltsville-human-nutrition-research-center/food-surveys-research-group/docs/fped-databases/ 
	You will need all the food pattern equivalents for foods in the WWEIA, NHANES 20XX-XX

## STEP 2 - Run MakeFoodGroups_nhanes.sas  - run for each survey cycle separately. Make sure the order of variables match with the cycle run. (e.g. order of food variables changed in 2017-2018 cycle) 
	*Input files:			FPED_DR1TOT_XXYY.sas7bdat
					FPED_DR2TOT_XXYY.sas7bdat
					FPED_DR1IFF_XXYY.sas7bdat
					FPED_DR2IFF_XXYY.sas7bdat

	* Output files:	FPED_DRAVGXXYY.sas7bdat

## STEP 3 - Run fped_CreateTerts.sas - Combine survey cycles to a pooled dataset, drop "total" food labels, add diet weights. XXYY denote a single survey cycle (e.g. 1112 for 2011-2012 cycle). XXZZ is the full range of combined surveys (e.g. 1118 for 2011-2018 surveys pooled).
	*Input files: 	DEMO_X.sas7bdat
				FPED_DRAVGXXYY.sas7bdat
				DR1TOT_X.sas7bdat

	*Output files:	tert_nhanesXXZZ.sas7bdat

## STEP 4 - Run GenerateHEI.sas - calculate HEI2015 score for each survey cycle.
	*Input files:		DEMO_X.sas7bdat
				FPED_DR1TOT_X.sas7bdat
				FPED_DR2TOT_X.sas7bdat
				hei2015_score_macro.sas

	*Output files:	HEIXXYY_avg.sas7bdat

## STEP 5 - Run MergeFPED_HEI.sas - Merge FPED-tertiles and HEI2015 score and demographic data to single dataset and export to CSV file
	*Input files:		HEIXXYY_avg.sas7bdat
				tert_nhanesXXZZ.sas7bdat

	*Output files:	hei_tert1118.sas7bdat
				adultHEI_fped1118terts.csv

## STEP 6 - Create_LowFdata.sas - Subset population to female adult participants living at or below the 130% poverty income level  
	Input files:		hei_tert1118.sas7bdat
	Output files:	nhanes1118_lowFdietdata.sas7bdat
				nhanes_lowFadultdata1Dec2021.csv

## STEP 7 - RPC_nhanesLOWF.m - Run robust profile clustering on female adults at or below 130% poverty level.
	*Input files:		nhanes_lowFadultdata1Dec2021.csv
	*Output files:	RPCnhanesLOWF_MCMCmed.mat
				NHANESLowFrpc_assign1Dec2021.csv
	*Output figures: nhanesRPC_dendrogram.png
				noconsum_lowadultpat.fig
				highconsum_lowadultpat.fig
				theta0_lowFpattern.fig
				nu_lowwomen.fig

## STEP 8 - nhanesRPC_tables.R
	*Input files:		NHANESLowFrpc_assign1Dec2021.csv
				nhanesadult_cvdriskHEIFHR1118.RData
				nhanes_lowFadultdata1Dec2021.csv 
	*Output files: 	NHANESLowF_indata_rpc1Dec2021.csv
