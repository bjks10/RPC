setwd("/Users/brianajoy/OneDrive\ -\ Harvard\ University/NHANES_SRP/Aim2")

# Supporting Libraries #
library(haven)
library(tidyverse)
library(survey)
library(psych)
library(readxl)
library(sas7bdat)

heiadult.nhanes <- read.sas7bdat('/Users/brianajoy/OneDrive\ -\ Harvard\ University/Migrated-P-Drive/NHANES/fped/hei_tert1118.sas7bdat')
nhanes_hei <- heiadult.nhanes[,c(1,4:17)]


#Calculate CVD Risk scores: 

#Data input: age, race, gender, marital status, education, hh size, bmi sbp hdl totchol bp_med smoker diabetes
demog2011 <-read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2011-2012/DEMO_G.XPT', col_select = c(SEQN, RIAGENDR,RIDAGEYR,RIDRETH3,DMDEDUC2,DMDHHSIZ,INDFMPIR,WTINT2YR,WTMEC2YR,SDMVPSU,SDMVSTRA))
demog2011 <- demog2011[ which(demog2011$RIDAGEYR>=20), ]

demog2013 <-read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2013-2014/DEMO_H.XPT', col_select = c(SEQN, RIAGENDR,RIDAGEYR,RIDRETH3,DMDEDUC2,DMDHHSIZ,INDFMPIR,WTINT2YR,WTMEC2YR,SDMVPSU,SDMVSTRA))
demog2013 <- demog2013[ which(demog2013$RIDAGEYR>=20 ), ]

demog2015 <-read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/DEMO_I.XPT', col_select = c(SEQN, RIAGENDR,RIDAGEYR,RIDRETH3,DMDEDUC2,DMDHHSIZ,INDFMPIR,WTINT2YR,WTMEC2YR,SDMVPSU,SDMVSTRA))
demog2015 <- demog2015[ which(demog2015$RIDAGEYR>=20 ), ]

demog2017 <-read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/DEMO_J.XPT', col_select = c(SEQN, RIAGENDR,RIDAGEYR,RIDRETH3,DMDEDUC2,DMDHHSIZ,INDFMPIR,WTINT2YR,WTMEC2YR,SDMVPSU,SDMVSTRA))
demog2017 <- demog2017[ which(demog2017$RIDAGEYR>=20 ), ]

demog.all <- rbind(demog2011,demog2013,demog2015,demog2017)
demog.all$WTINT8YR <- demog.all$WTINT2YR/4
demog.all$WTMEC8YR <- demog.all$WTMEC2YR/4
raceth <- c("Mexican","Other Hispanic", "NH White", "NH Black", "NH Asian", "Other/Mixed")
demog.all$RIDRETH3<- factor(demog.all$RIDRETH3, levels=c(1:4,6,7), labels=raceth)

demog.all$RIAGENDR <- factor(demog.all$RIAGENDR, levels=c(1,2), labels=c("male", "female"))
demog.all$EDU <- ifelse(demog.all$DMDEDUC2 < 3,"less than HS",ifelse(demog.all$DMDEDUC2==3,"HS/GED",ifelse(demog.all$DMDEDUC2>3,"At least some college","NA")))
demog.all$poverty <- ifelse(demog.all$INDFMPIR <=1.3, 1,0)
demog.all$poverty <- factor(demog.all$poverty,levels=c(0,1),labels=c("Above","At or Below"))
                            

#PULL NHANES BMI READING 
bmi_rdg11 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2011-2012/BMX_G.XPT',col_select = c(SEQN,BMXBMI))
bmi_rdg13 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2013-2014/BMX_H.XPT',col_select = c(SEQN,BMXBMI))
bmi_rdg15 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/BMX_I.XPT',col_select = c(SEQN,BMXBMI))
bmi_rdg17 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/BMX_J.XPT',col_select = c(SEQN,BMXBMI))

bmi.rdg <- rbind(bmi_rdg11,bmi_rdg13, bmi_rdg15, bmi_rdg17)
bmi.rdg$obese <- ifelse(bmi.rdg$BMXBMI>30,1,0)

#PULL NHANES BLOOD PRESSURE READINGS 
bp_rdg11 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2011-2012/BPX_G.XPT', col_select = c(SEQN, BPXSY1, BPXDI1,BPXSY2,BPXDI2,BPXSY3,BPXDI3,BPXSY4,BPXDI4))
bp_rdg13 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2013-2014/BPX_H.XPT', col_select = c(SEQN, BPXSY1, BPXDI1,BPXSY2,BPXDI2,BPXSY3,BPXDI3,BPXSY4,BPXDI4))
bp_rdg15 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/BPX_I.XPT', col_select = c(SEQN, BPXSY1, BPXDI1,BPXSY2,BPXDI2,BPXSY3,BPXDI3,BPXSY4,BPXDI4))
bp_rdg17 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/BPX_J.XPT', col_select = c(SEQN, BPXSY1, BPXDI1,BPXSY2,BPXDI2,BPXSY3,BPXDI3,BPXSY4,BPXDI4))

bp.all <- rbind(bp_rdg11,bp_rdg13, bp_rdg15,bp_rdg17)
bp.all$bpxsyavg <- rowMeans(cbind(bp.all$BPXSY1,bp.all$BPXSY2,bp.all$BPXSY3,bp.all$BPXSY4), na.rm=TRUE)
bp.all$bpxdiavg <-rowMeans(cbind(bp.all$BPXDI1,bp.all$BPXDI2,bp.all$BPXDI3,bp.all$BPXDI4), na.rm=TRUE)
bp.all$bp.bin <- ifelse(bp.all$bpxsyavg>140 | bp.all$bpxdiavg>90,1,0)
bp.rdg <- as.data.frame(cbind(bp.all$SEQN, bp.all$bpxsyavg, bp.all$bpxdiavg, bp.all$bp.bin))
colnames(bp.rdg) <- c("SEQN","SBP_avg","DBP_avg","HighBP_rdg")

#Taking bloodp ressure medication now 
bp_med11 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2011-2012/BPQ_G.XPT', col_select = c(SEQN, BPQ050A))
bp_med13 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2013-2014/BPQ_H.XPT', col_select = c(SEQN, BPQ050A))
bp_med15 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/BPQ_I.XPT', col_select = c(SEQN, BPQ050A))
bp_med17 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/BPQ_J.XPT', col_select = c(SEQN, BPQ050A))

bp_meds <- rbind(bp_med11,bp_med13, bp_med15,bp_med17)

#Pull HDL cholesterol and Total cholesterol
hdl_rdg11 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2011-2012/HDL_G.XPT',col_select = c(SEQN,LBDHDD))
hdl_rdg13 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2013-2014/HDL_H.XPT',col_select = c(SEQN,LBDHDD))
hdl_rdg15 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/HDL_I.XPT',col_select = c(SEQN,LBDHDD))
hdl_rdg17 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/HDL_J.XPT',col_select = c(SEQN,LBDHDD))

hdl.rdg <- rbind(hdl_rdg11,hdl_rdg13,hdl_rdg15,hdl_rdg17)

ldl_rdg11 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2011-2012/TRIGLY_G.XPT',col_select = c(SEQN,LBDLDL))
ldl_rdg13 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2013-2014/TRIGLY_H.XPT',col_select = c(SEQN,LBDLDL))
ldl_rdg15 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/TRIGLY_I.XPT',col_select = c(SEQN,LBDLDL))
ldl_rdg17 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/TRIGLY_J.XPT',col_select = c(SEQN,LBDLDL))

ldl.rdg <- rbind(ldl_rdg11,ldl_rdg13,ldl_rdg15,ldl_rdg17)


totchol_rdg11 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2011-2012/TCHOL_G.XPT',col_select = c(SEQN,LBXTC))
totchol_rdg13 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2013-2014/TCHOL_H.XPT',col_select = c(SEQN,LBXTC))
totchol_rdg15 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/TCHOL_I.XPT',col_select = c(SEQN,LBXTC))
totchol_rdg17 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/TCHOL_J.XPT',col_select = c(SEQN,LBXTC))

totchol.rdg <- rbind(totchol_rdg11,totchol_rdg13,totchol_rdg15,totchol_rdg17)

# NHANES - smoking status 
smoke_rdg11 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2011-2012/SMQ_G.XPT',col_select = c(SEQN,SMQ040))
smoke_rdg13 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2013-2014/SMQ_H.XPT',col_select = c(SEQN,SMQ040))
smoke_rdg15 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/SMQ_I.XPT',col_select = c(SEQN,SMQ040))
smoke_rdg17 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/SMQ_J.XPT',col_select = c(SEQN,SMQ040))

smoke.rdg <- rbind(smoke_rdg11,smoke_rdg13,smoke_rdg15, smoke_rdg17)
smoke.rdg$smoker <- ifelse(smoke.rdg$SMQ040==1,1,0)
smoke.rdg$smoker[is.na(smoke.rdg$smoker)] <- 0 



#NHANES - diabetes status 
glu_rdg11 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2011-2012/GLU_G.XPT',col_select = c(SEQN,LBXGLU))
glu_rdg13 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2013-2014/GLU_H.XPT',col_select = c(SEQN,LBXGLU))
glu_rdg15 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/GLU_I.XPT',col_select = c(SEQN,LBXGLU))
glu_rdg17 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/GLU_J.XPT',col_select = c(SEQN,LBXGLU))

glu.rdg <- rbind(glu_rdg11,glu_rdg13,glu_rdg15, glu_rdg17)
glu.rdg$glu.bin <- ifelse(glu.rdg$LBXGLU>126,1,0)

glu_qn11 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2011-2012/DIQ_G.XPT',col_select = c(SEQN,DIQ010,DIQ050,DIQ070))
glu_qn13 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2013-2014/DIQ_H.XPT',col_select = c(SEQN,DIQ010,DIQ050,DIQ070))
glu_qn15 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/DIQ_I.XPT',col_select = c(SEQN,DIQ010,DIQ050,DIQ070))
glu_qn17 <- read_xpt('https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/DIQ_J.XPT',col_select = c(SEQN,DIQ010,DIQ050,DIQ070))

glu.qn <- rbind(glu_qn11, glu_qn13, glu_qn15, glu_qn17)
glu.qn$glu.qflag <- ifelse(glu.qn$DIQ010==1 | glu.qn$DIQ050==1 | glu.qn$DIQ070==1, 1, 0)


glu.qnrdg <- merge(x=glu.rdg,y=glu.qn, by="SEQN", all.x=TRUE)
glu.qnrdg$diabetes <- ifelse(glu.qnrdg$glu.bin==1 | glu.qnrdg$glu.qflag==1,1,0)
glu.qnrdg$diabetes[is.na(glu.qnrdg$diabetes)] <- 0


## combine all datasets together ##
cvdrisk <-  bmi.rdg %>% full_join(bp.rdg, by = "SEQN") %>%
                            full_join(bp_meds, by = "SEQN") %>% full_join(hdl.rdg, by="SEQN") %>%
                            full_join(totchol.rdg, by ="SEQN") %>% full_join(smoke.rdg, by="SEQN") %>%
                            full_join(glu.qnrdg, by= "SEQN") %>%  full_join(ldl.rdg, by= "SEQN")

nhanes.cvdrisk <- merge(demog.all,cvdrisk,all.x=TRUE, all.y=FALSE)

#CHECK NA VALUES 
nhanes.cvdrisk$BPQ050A[is.na(nhanes.cvdrisk$BPQ050A)] <-0
nhanes.cvdrisk$diabetes[is.na(nhanes.cvdrisk$diabetes)] <-0


#Calculate AHA/ACC 2013 ASCVD risk scores
library(devtools)

install_github("vcastro/CVRisk")
library(CVrisk)


library(dplyr)
library(survey)
library(jtools)
#Framingham 2008 ASCVD risk score = ascvd_10y_frs
#Framingham 2008 ASCVD risk score simple (no lab) = ascvd_10y_frs_simple
nhanes_frscores <- compute_CVrisk(nhanes.cvdrisk, scores="ascvd_10y_frs",age = "RIDAGEYR", race = "RIDRETH3", gender ="RIAGENDR", bmi = "BMXBMI", 
                 sbp = "SBP_avg", hdl = "LBDHDD", totchol = "LBXTC", bp_med = "BPQ050A", 
                 smoker = "smoker", diabetes = "diabetes")

write.csv(nhanes_frscores,file='nhanes_frscores1118.csv')

nhanescvd.hei <- merge(nhanes_frscores,nhanes_hei,by="SEQN", all.x=TRUE,all.y=FALSE)

nhanescvd.hei$chol.flag <- ifelse(nhanescvd.hei$LBDLDL>150|nhanescvd.hei$LBXTC>200,1,0)
nhanescvd.hei$bp.flag <- ifelse(nhanescvd.hei$HighBP_rdg==1 | nhanescvd.hei$BPQ050A==1,1,0)

table(nhanescvd.hei$chol.flag)
table(nhanescvd.hei$bp.flag)
table(nhanescvd.hei$obese)
table(nhanescvd.hei$diabetes)
table(nhanescvd.hei$smoker)

table(nhanescvd.hei$chol.flag,nhanescvd.hei$poverty)
table(nhanescvd.hei$bp.flag,nhanescvd.hei$poverty)
table(nhanescvd.hei$obese,nhanescvd.hei$poverty)
table(nhanescvd.hei$diabetes,nhanescvd.hei$poverty)
table(nhanescvd.hei$smoker,nhanescvd.hei$poverty)

nhanescvd.hei$sumrisk <- rowSums(nhanescvd.hei[,c(16,24,31,48,49)], na.rm=TRUE)
nhanescvd.hei$sumrisk <- factor(nhanescvd.hei$sumrisk)
table(nhanescvd.hei$sumrisk)


nhanescvd.svy <- svydesign(id = ~SDMVPSU,
                           weights = ~WTMEC8YR,
                           strata=~SDMVSTRA,
                           nest = TRUE,
                           data = nhanescvd.hei)

#Overall population
sum(!is.na(nhanescvd.hei$HEI2015_TOTAL_SCORE))
svymean(~HEI2015_TOTAL_SCORE, nhanescvd.svy, na.rm=TRUE)
svyby(~HEI2015_TOTAL_SCORE,~poverty,nhanescvd.svy,svymean, na.rm=TRUE)
svyby(~HEI2015_TOTAL_SCORE,~RIDRETH3,nhanescvd.svy,svymean, na.rm=TRUE)
svyby(~HEI2015_TOTAL_SCORE,~poverty+RIDRETH3,nhanescvd.svy,svymean, na.rm=TRUE)

svymean(~ascvd_10y_frs,nhanescvd.svy, na.rm=TRUE)
svyby(~ascvd_10y_frs,~poverty,nhanescvd.svy,svymean, na.rm=TRUE)
svyby(~ascvd_10y_frs,~RIDRETH3,nhanescvd.svy,svymean, na.rm=TRUE)
svyby(~ascvd_10y_frs,~poverty+RIDRETH3,nhanescvd.svy,svymean, na.rm=TRUE)


# STRATIFY BY GENDER
svyby( ~HEI2015_TOTAL_SCORE,~RIAGENDR,nhanescvd.svy,svymean,na.rm=TRUE)
svyby( ~HEI2015_TOTAL_SCORE,~RIAGENDR+poverty,nhanescvd.svy,svymean,na.rm=TRUE)
svyby( ~HEI2015_TOTAL_SCORE,~RIAGENDR+RIDRETH3,nhanescvd.svy,svymean,na.rm=TRUE)


svyby( ~ascvd_10y_frs,~RIAGENDR+poverty,nhanescvd.svy,svymean,na.rm=TRUE)
svyby( ~ascvd_10y_frs,~RIAGENDR+RIDRETH3,nhanescvd.svy,svymean,na.rm=TRUE)

svycor(~ascvd_10y_frs+HEI2015_TOTAL_SCORE, nhanescvd.svy, na.rm=TRUE)
svycor(~bp.flag+bmi.flag+glu.flag+chol.flag, nhanescvd.svy, na.rm=TRUE)

svymean(~bp.flag+chol.flag+obese+diabetes+smoker, nhanescvd.svy, na.rm=TRUE)
svyby(~bp.flag+chol.flag+obese+diabetes+smoker, ~poverty,nhanescvd.svy,svymean, na.rm=TRUE)

svymean(~sumrisk,nhanescvd.svy,na.rm=TRUE)
svyby(~sumrisk,~poverty,nhanescvd.svy,svymean,na.rm=TRUE)
#svyby(~bp.flag+chol.flag+bmi.flag+glu.flag, ~RIDRETH3,nhanescvd.svy,svymean, na.rm=TRUE)

save(nhanescvd.hei,file="nhanesadult_cvdriskHEIFHR1118.RData")
write.csv(nhanescvd.hei,"nhanesadult_cvdriskHEIFHR1118.csv")


