setwd("/Users/brianajoy/OneDrive\ -\ Harvard\ University/NHANES_SRP/Aim3")

# Supporting Libraries #
library(haven)
library(tidyverse)
library(survey)
library(psych)

rpc.ci <- read.csv("NHANESLowFrpc_assign1Dec2021.csv")
rpc.ci <- rpc.ci %>% rename(
  SEQN = seqn_rpc
)

load('/Users/brianajoy/OneDrive\ -\ Harvard\ University/NHANES_SRP/Aim2/nhanesadult_cvdriskHEIFHR1118.RData')
nhanes.cvd <- nhanescvd.hei[,c(1:4,8:19,25,32,34:51)]
cvd.rpc <- merge(nhanes.cvd,rpc.ci,by="SEQN",all.x=FALSE,all.y=TRUE)

orig.nhanes <- read.csv("nhanes_lowFadultdata1Dec2021.csv")
fped.orig <- orig.nhanes[,c(1,18,23,37:65)]
fped.cvd.rpc <- merge(cvd.rpc,fped.orig, by="SEQN",all.x = TRUE, all.y = FALSE)

write.csv(fped.cvd.rpc,'NHANESLowF_indata_rpc1Dec2021.csv')
fped.cvd.rpc$MARSTAT <- ifelse(fped.cvd.rpc$DMDMARTL==1 | fped.cvd.rpc$DMDMARTL==6,1,ifelse(fped.cvd.rpc$DMDMARTL %in% c(2,3,4),2,ifelse(fped.cvd.rpc$DMDMARTL==5,3,NA)))
fped.cvd.rpc$MARSTAT <-factor(fped.cvd.rpc$MARSTAT,levels=c(1,2,3), labels=c("Married/LWP","Widowed-Divorced-Separated","Never Married"))

fped.cvd.rpc$agegroup <- ifelse(fped.cvd.rpc$RIDAGEYR>20 & fped.cvd.rpc$RIDAGEYR<35,1,
                                ifelse(fped.cvd.rpc$RIDAGEYR>=35 & fped.cvd.rpc$RIDAGEYR<50,2,
                                       ifelse(fped.cvd.rpc$RIDAGEYR>=50 & fped.cvd.rpc$RIDAGEYR<65,3,
                                              ifelse(fped.cvd.rpc$RIDAGEYR>=65,4, NA))))
fped.cvd.rpc$agegroup <- factor(fped.cvd.rpc$agegroup, levels=c(1:4),labels=c("20-34","35-49","50-64","65+"))

table(fped.cvd.rpc$RIDAGEYR)
table(fped.cvd.rpc$RIDRETH3)

table(fped.cvd.rpc$MARSTAT)
table(fped.cvd.rpc$EDU)

table(fped.cvd.rpc$bp.flag)
table(fped.cvd.rpc$diabetes)
table(fped.cvd.rpc$obese)
table(fped.cvd.rpc$chol.flag)
table(fped.cvd.rpc$smoker)


nhanesrpc.svy <- svydesign(id = ~SDMVPSU,
                           weights = ~dietwt8yr,
                           strata=~SDMVSTRA,
                           nest = TRUE,
                           data = fped.cvd.rpc)


# Generate Table 1
svymean(~RIDAGEYR+agegroup+RIDRETH3+MARSTAT+EDU,nhanesrpc.svy,na.rm=TRUE)
svymean(~HEI2015_TOTAL_SCORE+DRTKCAL+DMDHHSIZ+BMXBMI,nhanesrpc.svy,na.rm=TRUE)
svymean(~bp.flag+diabetes+obese+chol.flag+smoker,nhanesrpc.svy,na.rm=TRUE)

#Generate Table 2
svyby(~RIDAGEYR+agegroup+HEI2015_TOTAL_SCORE+DRTKCAL+DMDHHSIZ+INDFMPIR+MARSTAT+EDU+bp.flag+
        diabetes+obese+BMXBMI+chol.flag+smoker, ~RIDRETH3,nhanesrpc.svy,svymean,na.rm=TRUE)

#generate Table 3
svyby(~RIDAGEYR+agegroup+HEI2015_TOTAL_SCORE+DRTKCAL+DMDHHSIZ+INDFMPIR+RIDRETH3+MARSTAT+EDU+bp.flag+
        diabetes+obese+BMXBMI+chol.flag+smoker,~rpc_ci,nhanesrpc.svy,svymean,na.rm=TRUE)
svymean(~factor(rpc_ci),nhanesrpc.svy,na.rm=TRUE)



