### Merge DTI Data

library(reshape2)

demographics <- read.csv("/data/joy/BBL/projects/dtilongitudinal/data/n9498_demographics_go1_go2_113015.csv")
demographics <- demographics[!is.na(demographics$ageAtGo2Scan), ]

demographics <- demographics[c("bblid","race2","sex","ageAtGo1Scan","ageAtGo2Scan","meduCnbGo1","meduCnbGo2","handedness")]


demographics1 <- data_long <- melt(demographics[c("bblid","race2","sex","ageAtGo1Scan","ageAtGo2Scan","handedness")],
                  # ID variables - all the variables to keep but not split apart on
                  id.vars=c("bblid", "sex","race2","handedness"),
                  # The source columns
                  measure.vars=c("ageAtGo1Scan","ageAtGo2Scan"),
                  # Name of the destination column that will identify the original
                  # column that the measurement came from
                  variable.name="timepoint",
                  value.name="AgeAtScan")

demographics2 <- data_long <- melt(demographics[c("bblid","race2","sex","meduCnbGo1","meduCnbGo2","handedness")],
                                   # ID variables - all the variables to keep but not split apart on
                                   id.vars=c("bblid", "sex","race2","handedness"),
                                   # The source columns
                                   measure.vars=c("meduCnbGo1","meduCnbGo2"),
                                   # Name of the destination column that will identify the original
                                   # column that the measurement came from
                                   variable.name="timepoint",
                                   value.name="meduCnb")

demographics1$timepoint <- as.character(demographics1$timepoint)
demographics2$timepoint <- as.character(demographics2$timepoint)

demographics1$timepoint[which(demographics1$timepoint == "ageAtGo1Scan")] <- 1
demographics1$timepoint[which(demographics1$timepoint == "ageAtGo2Scan")] <- 2
demographics2$timepoint[which(demographics2$timepoint == "meduCnbGo1")] <- 1
demographics2$timepoint[which(demographics2$timepoint == "meduCnbGo2")] <- 2

demographics <- merge(demographics1, demographics2, by=c("bblid","race2","sex","handedness","timepoint"))

QA <- read.csv("/data/joy/BBL/projects/dtilongitudinal/data/n2416_dti_qa_20170301.csv")

finalDat <- merge(demographics, QA, by=c("bblid","timepoint"))

clinical <- read.csv("/data/joy/BBL/projects/dtilongitudinal/data/pnc_diagnosis_categorical_20170412.csv", na.strings="")

clinical$dx_t1_t2 <- as.factor(paste0(clinical$dx_t1_psychosis,"_",clinical$dx_t2_psychosis))
clinical$dx_t1_t2[which(clinical$dx_t1_t2 == "NA_NA")] <- NA

clinical$dx_t1_t2_clean <- as.factor(paste0(clinical$dx_t1_psychosis_clean,"_",clinical$dx_t2_psychosis_clean))
clinical$dx_t1_t2_clean[which(clinical$dx_t1_t2_clean == "NA_NA")] <- NA

clinical$dx_pscat_t2_psyispro <- as.character(clinical$dx_pscat_t2)
clinical$dx_pscat_t2_psyispro[which(clinical$dx_pscat_t2_psyispro == "psy")] <- "pro"



### I will use dx_expanded_psyispro_t1_t2 to build squeaky cleans vs. PSPS
clinical$dx_expanded_t1_t2 <- paste0(clinical$goassessDxpmr6,"_",clinical$dx_pscat_t2)
clinical$dx_expanded_psyispro_t1_t2 <- paste0(clinical$goassessDxpmr6,"_",clinical$dx_pscat_t2_psyispro)

clinical$dx_expanded_t1_t2[which(clinical$dx_expanded_t1_t2 == "NA_NA")] <- NA
clinical$dx_expanded_psyispro_t1_t2[which(clinical$dx_expanded_psyispro_t1_t2 == "NA_NA")] <- NA


finalDat <- merge(finalDat, clinical, by="bblid", all.x=T)

##Health Excludes

health1601 <- read.csv("/data/joy/BBL/projects/dtilongitudinal/data/n1601_health_20161214.csv")
health1601 <- health1601[,-(2:3)]
finalDat <- merge(finalDat, health1601, by="bblid")

health404 <- read.csv("/data/joy/BBL/projects/dtilongitudinal/data/n404_medication_exclusions.csv")

finalDat <- merge(finalDat, health404, all.x=T, by=c("bblid","timepoint","scanid"))

finalDat <- finalDat[finalDat$dti64Exclude == 0,]
finalDat <- finalDat[!is.na(finalDat$dx_t1_t2), ]
finalDat <- finalDat[finalDat$dti64Tsnr > 6.47, ]
finalDat <- finalDat[finalDat$AgeAtScan > 11*12, ]
#finalDat <- finalDat[which(finalDat$healthExcludev2 == 0), ]
#finalDat <- finalDat[-which(finalDat$psychoactiveMedMedical_20160930 == 1), ]
#finalDat <- finalDat[!(finalDat$timepoint == 2 & finalDat$dx_t1_t2 == "Healthy" & finalDat$psychoactiveMedPsych_20160930 == 1), ]

finalDat <- finalDat[which(finalDat$dx_t1_t2_clean == "PS_PS" | 
                             finalDat$dx_t1_t2_clean == "PS_NONPS" |
                             finalDat$dx_t1_t2_clean == "NONPS_PS" |
                             finalDat$dx_t1_t2_clean == "NONPS_NONPS"), ]

write.csv(finalDat, "/data/joy/BBL/projects/dtilongitudinal/data/n404_dtiAnalysis_finalListbblidscan.csv")
