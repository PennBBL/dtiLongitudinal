library(tableone)

demo <- read.csv("/data/joy/BBL/projects/dtilongitudinal/data/n404_dtiAnalysis_finalListbblidscan_withClinical.csv")
vars <- c("race2","sex","handedness","AgeAtScan","meduCnb",
          "positive","negative",
          "disorganized","general",
          "sevSipsMean","SumSips",
          "feduCnb","gaf_z","gaf_raw","dti64MeanRelRMS","dti64Tsnr")

demo$sex <- as.factor(demo$sex)
demo$race2[which(demo$race2 == 3)] <- 2
demo$race2 <- as.factor(demo$race2)
demo$handedness <- as.factor(demo$handedness)
demo$AgeAtScan <- demo$AgeAtScan / 12




##Table 1
demo1 <- demo[demo$timepoint == 1, ]

## Create Table 1 stratified by trt
tableOne <- CreateTableOne(vars = vars, strata = c("dx_t1_t2_clean"), data = demo1)

t <- print(tableOne, argsApprox = c("sex","race2"), test=T, noSpaces=T)

write.csv(t, "/data/joy/BBL/projects/dtilongitudinal/output/dtiLongitudinal_table1_timepoint1.csv")

### Table 2
demo2 <- demo[demo$timepoint == 2, ]

## Create Table 2 stratified by trt
tableOne2 <- CreateTableOne(vars = vars, strata = c("dx_t1_t2_clean"), data = demo2)

t2 <- print(tableOne2, argsApprox = c("sex","race2"), test=T, noSpaces=T)

write.csv(t2, "/data/joy/BBL/projects/dtilongitudinal/output/dtiLongitudinal_table1_timepoint2.csv")


##Suplemment Tables

demo <- demo[which(demo$dx_expanded_psyispro_t1_t2 == "PS_pro" | demo$dx_expanded_psyispro_t1_t2 == "TD_none"), ]
demo$dx_expanded_psyispro_t1_t2 <- as.factor(as.character(demo$dx_expanded_psyispro_t1_t2))
demo1 <- demo[demo$timepoint == 1, ]

## Create Table 1 stratified by trt
tableOneSup <- CreateTableOne(vars = vars, strata = c("dx_expanded_psyispro_t1_t2"), data = demo1)
tSup <- print(tableOneSup, argsApprox = c("sex","race2"), test=T, noSpaces=T)

write.csv(tSup, "/data/joy/BBL/projects/dtilongitudinal/output/dtiLongitudinal_supplementtable1_timepoint1.csv")

##Timepoint 2
demo2 <- demo[demo$timepoint == 2, ]
## Create Table 1 stratified by trt
tableOne2Sup <- CreateTableOne(vars = vars, strata = c("dx_expanded_psyispro_t1_t2"), data = demo2)

t2Sup <- print(tableOne2Sup, argsApprox = c("sex","race2"), test=T, noSpaces=T)
write.csv(t2Sup, "/data/joy/BBL/projects/dtilongitudinal/output/dtiLongitudinal_supplementtable1_timepoint2.csv")

