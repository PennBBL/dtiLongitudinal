## Create Function to get Statistics

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



##Declate objects for statistics function

contvars <- c("AgeAtScan","meduCnb",
              "positive","negative",
              "disorganized","general",
              "sevSipsMean","SumSips",
              "feduCnb","gaf_z","gaf_raw","dti64MeanRelRMS","dti64Tsnr")
catvars <- c("race2","sex","handedness")
group = "dx_t1_t2_clean"


#Load and Clean Demographic data

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

### Table 1 Statistics

demoStat1 <- test.statistics(catvars, contvars, group, data=demo1)
demoCat1 <- demoStat1[[1]]
demoCont1 <- demoStat1[[2]]

t <- as.data.frame(t)
t <- t[,-c(5:6)]
t[,5:9] <- NA

t[2:4,5:8] <- demoCat1
t[5:dim(t)[2],5:9] <- demoCont1
t <- t[,-5]

names(t)[5:8] <- c("pval","statistics","numdf","dendf")

t$pval2 <- as.character(round(as.numeric(t$pval), 3))
t$pval2[which(as.numeric(t$pval2) < 0.001)] <- "<0.001"


write.csv(t, "/data/joy/BBL/projects/dtilongitudinal/output/dtiLongitudinal_table1_timepoint1.csv")


### Table 2
demo2 <- demo[demo$timepoint == 2, ]

## Create Table 2 stratified by trt
tableOne2 <- CreateTableOne(vars = vars, strata = c("dx_t1_t2_clean"), data = demo2)
t2 <- print(tableOne2, argsApprox = c("sex","race2"), test=T, noSpaces=T)

### Table 2 Statistics

demoStat2 <- test.statistics(catvars, contvars, group, data=demo2)
demoCat2 <- demoStat2[[1]]
demoCont2 <- demoStat2[[2]]

t2 <- as.data.frame(t2)
t2 <- t2[,-c(5:6)]
t2[,5:9] <- NA

t2[2:4,5:8] <- demoCat2
t2[5:dim(t2)[1],5:9] <- demoCont2
t2 <- t2[,-5]

names(t2)[5:8] <- c("pval","statistics","numdf","dendf")

t2$pval2 <- as.character(round(as.numeric(t2$pval), 3))
t2$pval2[which(as.numeric(t2$pval2) < 0.001)] <- "<0.001"
write.csv(t2, "/data/joy/BBL/projects/dtilongitudinal/output/dtiLongitudinal_table1_timepoint2.csv")


##Suplemment Tables

demo <- demo[which(demo$dx_expanded_psyispro_t1_t2 == "PS_pro" | demo$dx_expanded_psyispro_t1_t2 == "TD_none"), ]
demo$dx_expanded_psyispro_t1_t2 <- as.factor(as.character(demo$dx_expanded_psyispro_t1_t2))
demo1 <- demo[demo$timepoint == 1, ]

## Create Table 1 stratified by trt
tableOneSup <- CreateTableOne(vars = vars, strata = c("dx_expanded_psyispro_t1_t2"), data = demo1)
tSup <- print(tableOneSup, argsApprox = c("sex","race2"), test=T, noSpaces=T)

### Table 1 Statistics Sup
group <- "dx_expanded_psyispro_t1_t2"
demoStatSup <- test.statistics(catvars, contvars, group, data=demo1)
demoCatSup <- demoStatSup[[1]]
demoContSup <- demoStatSup[[2]]

tSup <- as.data.frame(tSup)
tSup <- tSup[,-c(3:4)]
tSup[,3:7] <- NA

tSup[2:4,3:6] <- demoCatSup
tSup[5:dim(tSup)[1],3:6] <- demoContSup
tSup <- tSup[,-c(3,7)]

names(tSup)[3:5] <- c("pval","statistics","df")

tSup$pval2 <- as.character(round(as.numeric(tSup$pval), 3))
tSup$pval2[which(as.numeric(tSup$pval2) < 0.001)] <- "<0.001"

write.csv(tSup, "/data/joy/BBL/projects/dtilongitudinal/output/dtiLongitudinal_supplementtable1_timepoint1.csv")



##Timepoint 2
demo2 <- demo[demo$timepoint == 2, ]
## Create Table 1 stratified by trt
tableOne2Sup <- CreateTableOne(vars = vars, strata = c("dx_expanded_psyispro_t1_t2"), data = demo2)
t2Sup <- print(tableOne2Sup, argsApprox = c("sex","race2"), test=T, noSpaces=T)

### Table 1 Statistics Sup
group <- "dx_expanded_psyispro_t1_t2"
demoStat2Sup <- test.statistics(catvars, contvars, group, data=demo2)
demoCat2Sup <- demoStat2Sup[[1]]
demoCont2Sup <- demoStat2Sup[[2]]

t2Sup <- as.data.frame(t2Sup)
t2Sup <- t2Sup[,-c(3:4)]
t2Sup[,3:7] <- NA

t2Sup[2:4,3:6] <- demoCat2Sup
t2Sup[5:dim(t2Sup)[1],3:6] <- demoCont2Sup
t2Sup <- t2Sup[,-c(3,7)]

names(t2Sup)[3:5] <- c("pval","statistics","df")

t2Sup$pval2 <- as.character(round(as.numeric(t2Sup$pval), 3))
t2Sup$pval2[which(as.numeric(t2Sup$pval2) < 0.001)] <- "<0.001"


write.csv(t2Sup, "/data/joy/BBL/projects/dtilongitudinal/output/dtiLongitudinal_supplementtable1_timepoint2.csv")

