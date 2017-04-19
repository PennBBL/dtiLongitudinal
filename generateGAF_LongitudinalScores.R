##Generate GAF Longitudinal Values


## Clinical Go1 GAF Score
clinicalGo1Gaf <- read.csv("/data/joy/BBL/projects/dtilongitudinal/data/n1601_goassess_gaf_scores_20170308.csv")
diagnosis <- read.csv("/data/joy/BBL/studies/pnc/n1601_dataFreeze/clinical/n1601_go1_diagnosis_dxpmr_20161014.csv")

clinicalGo1Gaf <- merge(clinicalGo1Gaf, diagnosis, by=c("bblid","scanid"))

clinicalGo1Gaf <- clinicalGo1Gaf[clinicalGo1Gaf$goassessDxpmr6 == "TD", ]

mean(clinicalGo1Gaf$gaf001)
sd(clinicalGo1Gaf$gaf001)


##Clinical Go2 GAF Score
clinicalGo2Gaf <- read.csv("/data/joy/BBL/projects/dtilongitudinal/data/n404_sops_gaf_20170308.csv")
diagnosis2 <- read.csv("/data/joy/BBL/projects/dtilongitudinal/data/pnc_diagnosis_categorical_20170412.csv")

diagnosis2$dx_t1_t2 <- as.factor(paste0(diagnosis2$dx_t1_psychosis,"_",diagnosis2$dx_t2_psychosis))
diagnosis2$dx_t1_t2[which(diagnosis2$dx_t1_t2 == "NA_NA")] <- NA

clinicalGo2Gaf <- merge(clinicalGo2Gaf, diagnosis2, by=c("bblid"))

clinicalGo2Gaf <- clinicalGo2Gaf[clinicalGo2Gaf$dx_t1_t2 == "NONPS_NONPS", ]
clinicalGo2Gaf <- clinicalGo2Gaf[!is.na(clinicalGo2Gaf$gaf_c), ]

mean(clinicalGo2Gaf$gaf_c)
sd(clinicalGo2Gaf$gaf_c)


## Clinical Go1 GAF Score
Go1Gaf <- read.csv("/data/joy/BBL/projects/dtilongitudinal/data/n1601_goassess_gaf_scores_20170308.csv")
Go1Gaf$gaf001 <- (Go1Gaf$gaf001 - mean(clinicalGo1Gaf$gaf001)) / sd(clinicalGo1Gaf$gaf001)
Go1Gaf$timepoint <- 1

## Clinical Go1 GAF Score
Go2Gaf <- read.csv("/data/joy/BBL/projects/dtilongitudinal/data/n404_sops_gaf_20170308.csv")
Go2Gaf$gaf_c <- (Go2Gaf$gaf_c - mean(clinicalGo2Gaf$gaf_c)) / sd(clinicalGo2Gaf$gaf_c)
Go2Gaf$timepoint <- 2

Go1Gaf <- Go1Gaf[c("bblid","timepoint","gaf001") ]
Go2Gaf <- Go2Gaf[c("bblid","timepoint","gaf_c") ]

names(Go1Gaf)[3] <- "gaf_z"
names(Go2Gaf)[3] <- "gaf_z"

GafLong <- rbind(Go1Gaf, Go2Gaf)

write.csv(GafLong, "/data/joy/BBL/projects/dtilongitudinal/data/n2416_GAF-zScores_longitudinal.csv", row.names=F)

## Raw GAF Scores
Go1Gaf <- read.csv("/data/joy/BBL/projects/dtilongitudinal/data/n1601_goassess_gaf_scores_20170308.csv")
#Go1Gaf$gaf001 <- (Go1Gaf$gaf001 - mean(clinicalGo1Gaf$gaf001)) / sd(clinicalGo1Gaf$gaf001)
Go1Gaf$timepoint <- 1

## Clinical Go1 GAF Score
Go2Gaf <- read.csv("/data/joy/BBL/projects/dtilongitudinal/data/n404_sops_gaf_20170308.csv")
#Go2Gaf$gaf_c <- (Go2Gaf$gaf_c - mean(clinicalGo2Gaf$gaf_c)) / sd(clinicalGo2Gaf$gaf_c)
Go2Gaf$timepoint <- 2

Go1Gaf <- Go1Gaf[c("bblid","timepoint","gaf001") ]
Go2Gaf <- Go2Gaf[c("bblid","timepoint","gaf_c") ]

names(Go1Gaf)[3] <- "gaf_raw"
names(Go2Gaf)[3] <- "gaf_raw"

GafLong <- rbind(Go1Gaf, Go2Gaf)

write.csv(GafLong, "/data/joy/BBL/projects/dtilongitudinal/data/n2416_GAF-rawScores_longitudinal.csv", row.names=F)


