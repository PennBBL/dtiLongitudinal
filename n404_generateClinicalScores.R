### Merge DTI Data

library(reshape2)

demographics <- read.csv("/data/joy/BBL/projects/dtilongitudinal/data/n9498_demographics_go1_go2_113015.csv")
demographics <- demographics[!is.na(demographics$ageAtGo2Scan), ]

demographics <- demographics[c("bblid","race2","sex","ageAtGo1Scan","ageAtGo2Scan","feduCnbGo1","feduCnbGo2","handedness")]

demographics2 <- data_long <- melt(demographics[c("bblid","race2","sex","feduCnbGo1","feduCnbGo2","handedness")],
                                   # ID variables - all the variables to keep but not split apart on
                                   id.vars=c("bblid", "sex","race2","handedness"),
                                   # The source columns
                                   measure.vars=c("feduCnbGo1","feduCnbGo2"),
                                   # Name of the destination column that will identify the original
                                   # column that the measurement came from
                                   variable.name="timepoint",
                                   value.name="feduCnb")

demographics2$timepoint <- as.character(demographics2$timepoint)
demographics2$timepoint[which(demographics2$timepoint == "feduCnbGo1")] <- 1
demographics2$timepoint[which(demographics2$timepoint == "feduCnbGo2")] <- 2

demographics2 <- demographics2[,c("bblid","timepoint","feduCnb") ]

gamraw <- read.csv("/data/joy/BBL/projects/dtilongitudinal/data/n2416_GAF-rawScores_longitudinal.csv")
gamZ <- read.csv("/data/joy/BBL/projects/dtilongitudinal/data/n2416_GAF-zScores_longitudinal.csv")

clinical <- merge(demographics2, gamZ, by=c("bblid","timepoint"))
clinical <- merge(clinical,gamraw , by=c("bblid","timepoint"))

sips2 <- read.csv("/data/joy/BBL/projects/dtilongitudinal/data/n404_go2CAPA_sumsips_sevsips_scores.csv")
sips1 <- read.csv("/data/joy/BBL/projects/dtilongitudinal/data/n1601_goassess_prime_sops_scores_20150723.csv")

sips1 <- sips1[c("bblid","goassessSevSipMean","goassessSumSips")]
sips2 <- sips2[c("bblid","SEV_SIP_MN","SUM_SIPS")]

names(sips1) <- c("bblid","sevSipsMean","SumSips")
names(sips2) <- c("bblid","sevSipsMean","SumSips")

sips1$timepoint <- 1
sips2$timepoint <- 2

sips <- rbind(sips1, sips2)

sops <- read.csv("/data/joy/BBL/projects/dtilongitudinal/data/n404_sops_gaf_20170308.csv")

sops$positive <- rowSums(sops[,2:6])
sops$negative <- rowSums(sops[,7:12])
sops$disorganized <- rowSums(sops[,13:16])
sops$general <- rowSums(sops[,17:20])

sops <- sops[,-22]
sops$timepoint <- 2

mmse <- read.csv("/data/joy/BBL/projects/dtilongitudinal/data/n404_mmse_total_score_GOxCAPARepositoryAna2017-03-08.csv")

mmse$timepoint <- 2

demo <- read.csv("/data/joy/BBL/projects/dtilongitudinal/data/n404_dtiAnalysis_finalListbblidscan.csv")

demo <- merge(demo, mmse, by=c("bblid","timepoint"), all.x=T)
demo <- merge(demo, sops, by=c("bblid","timepoint"), all.x=T)
demo <- merge(demo, sips, by=c("bblid","timepoint"), all.x=T)
demo <- merge(demo, clinical, by=c("bblid","timepoint"), all.x=T)

write.csv(demo, "/data/joy/BBL/projects/dtilongitudinal/data/n404_dtiAnalysis_finalListbblidscan_withClinical.csv", row.names=F)



