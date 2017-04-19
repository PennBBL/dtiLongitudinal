### Merge DTI Data

library(reshape2)

demographics <- read.csv("/data/joy/BBL/projects/dtilongitudinal/data/n9498_demographics_go1_go2_113015.csv")
demographics <- demographics[!is.na(demographics$ageAtGo1Scan), ]
demographics <- demographics[c("bblid","race2","sex","ageAtGo1Scan","ageAtGo2Scan","meduCnbGo1","meduCnbGo2","handedness")]

demographics$timepoint <- 1

QA <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/neuroimaging/dti/n1601_dti_qa_20170301.csv")
QA$timepoint <- 1

finalDat <- merge(demographics, QA, by=c("bblid","timepoint"))

clinical <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/clinical/n1601_go1_diagnosis_dxpmr_20161014.csv")
clinical <- clinical[-which(clinical$goassessDxpmr6 == ""), ]

finalDat <- merge(finalDat, clinical, by=c("bblid","scanid"))

##Health Excludes

health1601 <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/health/n1601_health_20161214.csv")
health1601 <- health1601[,-(2:3)]
finalDat <- merge(finalDat, health1601, by="bblid")

finalDat <- finalDat[finalDat$dti64Exclude == 0,]
finalDat <- finalDat[finalDat$dti64Tsnr > 6.47, ]
finalDat <- finalDat[finalDat$ageAtGo1Scan > 11*12, ]
finalDat <- finalDat[which(finalDat$healthExcludev2 == 0), ]

write.csv(finalDat, "/data/joy/BBL/projects/dtilongitudinal/data/n1601_dtiAnalysis_finalListbblidscan.csv")
