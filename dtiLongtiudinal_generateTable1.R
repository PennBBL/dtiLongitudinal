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
demo$dx_t1_t2_clean <- as.factor(demo$dx_t1_t2_clean)



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


##Calculate Statistics 

contvars <- c("AgeAtScan","meduCnb",
          "positive","negative",
          "disorganized","general",
          "sevSipsMean","SumSips",
          "feduCnb","gaf_z","gaf_raw","dti64MeanRelRMS","dti64Tsnr")
catvars <- c("race2","sex","handedness")

group = "dx_t1_t2_clean"

data <- demo1

test.statistics <- function(catvars, contvars, group, data) {
  
  data[,which(names(data) == group)] <- as.factor(as.character(data[,which(names(data) == group)]))
  
  for (i in catvars) {
    if (all(is.na(data[which(names(data) == i)]))) {
      catvars <- catvars[-which(catvars == i )]
    }
  }
    
  for (i in contvars) {
    if (all(is.na(data[which(names(data) == i)]))) {
      contvars <- contvars[-which(contvars == i )]
    }
  }
  
  catStat <- t(mcmapply(function(x) {
    t <- chisq.test(data[, which(names(data) == x)],data[, which(names(data) == group)])
    c(x,t$p.value, t$statistic, t$parameter)
  }, catvars, SIMPLIFY=T))

  
  if (length(levels(data[,which(names(data) == group)])) > 2) {
    contStat <- t(mcmapply(function(x) {
      formula <- as.formula(paste(x, "~", group))
      t <- oneway.test(formula, data, var.equal=F)
      c(x,t$p.value, t$statistic, t$parameter)
    }, contvars, SIMPLIFY=T))
  } else if (length(levels(data[,which(names(data) == group)])) == 2) {
    contStat <- t(mcmapply(function(x) {
      formula <- as.formula(paste(x, "~", group))
      t <- t.test(formula, data, var.equal=F)
      c(x,t$p.value, t$statistic, t$parameter)
    }, contvars, SIMPLIFY=T))
  }
  
  row.names(contStat) <- NULL
  row.names(catStat) <- NULL
  
  return(list(catStat, contStat))
}

temp <- test.statistics(catvars, contvars, group, data = demo2)
