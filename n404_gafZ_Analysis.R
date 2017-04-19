library(parallel)
library(mgcv)
library(gamm4)

combineROI <- function(indata) {
  
  tmp <- indata
  
  ##############################################################################
  ################   Calculate average across left and right     ###############
  ##############################################################################
  
  tmp$mean_dti <- rowMeans(tmp[, grep("dti_", names(tmp))])
  
  ##############################################################################
  ################   Calculate average across left and right     ###############
  ##############################################################################
  
  index.roi <- grep("_l", names(tmp))
  
  for (i in index.roi) {
    #Compute mean across L-R
    tmp[,dim(tmp)[2] + 1] <- rowMeans(tmp[,c(i,i+1)], na.rm = T)
    
    #Create new name for column
    name.start <- as.vector(strsplit(names(tmp)[i], "_")[[1]])
    name.end <- paste0(name.start[4],"_",name.start[5])
    names(tmp)[dim(tmp)[2]] <- name.end
  }
  
  # Reorganize Dataset
  index <- c(grep("forceps", names(tmp)), grep("mean_", names(tmp)))
  temp <- tmp[, index]
  tmp <- cbind(tmp[, -index], temp)
  
  index <- c(grep("_l$", names(tmp)), grep("_r$", names(tmp)))
  
  tmp <- tmp[,-index]
  
  return(tmp)
}

model <- "~ s(AgeAtScan, k=4) + sex + race2+ gaf_z + dti64Tsnr"
var <- "dx_expanded_psyispro_t1_t2TD_none"
randomFor <- "~ (1|bblid)"
demo <- read.csv("/data/joy/BBL/projects/dtilongitudinal/data/n404_dtiAnalysis_finalListbblidscan_withClinical.csv")
#demo <- demo[which(demo$dx_expanded_psyispro_t1_t2 == "PS_pro" | demo$dx_expanded_psyispro_t1_t2 == "TD_none"), ]
demo$race2[which(demo$race2 == 3)] <- 2

### Analyze FA Data

path <- "/data/joy/BBL/studies/pnc/n2416_dataFreeze/neuroimaging/dti/n2416_DTI64/n2416_DTI_JHUTractFA.csv"
indata <- read.csv(path)
indata <- combineROI(indata)


data.analysis <- merge(demo, indata, by=c("bblid","scanid"))

covariates <- "~ s(AgeAtScan, k=4) + sex + race2 + gaf_z + dti64Tsnr "

gam_pval.fa = as.data.frame(matrix(NA, nrow=146, ncol=8))

for (i in 162:172) {
  
  var.name <- names(data.analysis)[i]
  m <- as.formula(paste(var.name, covariates,sep=""))
  gam1 <- gamm4(m, data=data.analysis, REML=T, random= ~ (1|bblid))
  
  n <- i 
  
  gam_pval.fa[n, 1] <- var.name
  gam_pval.fa[n, 2:(1 + dim(summary(gam1$gam)$p.table)[1])] <- summary(gam1$gam)$p.table[,4]
  gam_pval.fa[n, (2 + dim(summary(gam1$gam)$p.table)[1]):(1 
                                                          + dim(summary(gam1$gam)$p.table)[1] 
                                                          + dim(summary(gam1$gam)$s.table)[1]) ] <- summary(gam1$gam)$s.table[,4]
  
}

gam_pval.fa <- gam_pval.fa[!is.na(gam_pval.fa[,1]), ]
names(gam_pval.fa)[1] <- "ROI"
names(gam_pval.fa)[ 2:(1 + dim(summary(gam1$gam)$p.table)[1])] <- rownames(summary(gam1$gam)$p.table)
names(gam_pval.fa)[ (2 + dim(summary(gam1$gam)$p.table)[1]):(1 
                                                             + dim(summary(gam1$gam)$p.table)[1] 
                                                             + dim(summary(gam1$gam)$s.table)[1]) ] <- rownames(summary(gam1$gam)$s.table)

gam_corr.fa <- gam_pval.fa

for (i in 2:dim(gam_pval.fa)[2]) {
  gam_corr.fa[1:10,i] <- p.adjust(gam_pval.fa[1:10,i], method="bonferroni")
}

### Analyze MD Data

path <- "/data/joy/BBL/studies/pnc/n2416_dataFreeze/neuroimaging/dti/n2416_DTI64/n2416_DTI_JHUTractTR.csv"
indata <- read.csv(path)
indata <- combineROI(indata)


data.analysis <- merge(demo, indata, by=c("bblid","scanid"))

covariates <- "~ s(AgeAtScan, k=4) + sex + race2+ gaf_z + dti64Tsnr "

gam_pval.md = as.data.frame(matrix(NA, nrow=146, ncol=8))

for (i in 162:172) {
  
  var.name <- names(data.analysis)[i]
  m <- as.formula(paste(var.name, covariates,sep=""))
  gam1 <- gamm4(m, data=data.analysis, REML=T, random= ~ (1|bblid))
  
  n <- i 
  
  gam_pval.md[n, 1] <- var.name
  gam_pval.md[n, 2:(1 + dim(summary(gam1$gam)$p.table)[1])] <- summary(gam1$gam)$p.table[,4]
  gam_pval.md[n, (2 + dim(summary(gam1$gam)$p.table)[1]):(1 
                                                          + dim(summary(gam1$gam)$p.table)[1] 
                                                          + dim(summary(gam1$gam)$s.table)[1]) ] <- summary(gam1$gam)$s.table[,4]
  
}


gam_pval.md <- gam_pval.md[!is.na(gam_pval.md[,1]), ]
names(gam_pval.md)[1] <- "ROI"
names(gam_pval.md)[ 2:(1 + dim(summary(gam1$gam)$p.table)[1])] <- rownames(summary(gam1$gam)$p.table)
names(gam_pval.md)[ (2 + dim(summary(gam1$gam)$p.table)[1]):(1 
                                                             + dim(summary(gam1$gam)$p.table)[1] 
                                                             + dim(summary(gam1$gam)$s.table)[1]) ] <- rownames(summary(gam1$gam)$s.table)

gam_corr.md <- gam_pval.md

for (i in 2:dim(gam_pval.md)[2]) {
  gam_corr.md[1:10,i] <- p.adjust(gam_pval.md[1:10,i], method="bonferroni")
}



cor.test(data.analysis$gaf_z, data.analysis$mean_dti)


### Analyze AD Data

path <- "/data/joy/BBL/studies/pnc/n2416_dataFreeze/neuroimaging/dti/n2416_DTI64/n2416_DTI_JHUTractAD.csv"
indata <- read.csv(path)
indata <- combineROI(indata)


data.analysis <- merge(demo, indata, by=c("bblid","scanid"))

covariates <- "~ s(AgeAtScan, k=4) + sex + race2+ gaf_z + dti64Tsnr "

gam_pval.ad = as.data.frame(matrix(NA, nrow=146, ncol=8))

for (i in 162:172) {
  
  var.name <- names(data.analysis)[i]
  m <- as.formula(paste(var.name, covariates,sep=""))
  gam1 <- gamm4(m, data=data.analysis, REML=T, random= ~ (1|bblid))
  
  n <- i 
  
  gam_pval.ad[n, 1] <- var.name
  gam_pval.ad[n, 2:(1 + dim(summary(gam1$gam)$p.table)[1])] <- summary(gam1$gam)$p.table[,4]
  gam_pval.ad[n, (2 + dim(summary(gam1$gam)$p.table)[1]):(1 
                                                          + dim(summary(gam1$gam)$p.table)[1] 
                                                          + dim(summary(gam1$gam)$s.table)[1]) ] <- summary(gam1$gam)$s.table[,4]
  
}


gam_pval.ad <- gam_pval.ad[!is.na(gam_pval.ad[,1]), ]
names(gam_pval.ad)[1] <- "ROI"
names(gam_pval.ad)[ 2:(1 + dim(summary(gam1$gam)$p.table)[1])] <- rownames(summary(gam1$gam)$p.table)
names(gam_pval.ad)[ (2 + dim(summary(gam1$gam)$p.table)[1]):(1 
                                                             + dim(summary(gam1$gam)$p.table)[1] 
                                                             + dim(summary(gam1$gam)$s.table)[1]) ] <- rownames(summary(gam1$gam)$s.table)

gam_corr.ad <- gam_pval.ad

for (i in 2:dim(gam_pval.ad)[2]) {
  gam_corr.ad[1:10,i] <- p.adjust(gam_pval.ad[1:10,i], method="bonferroni")
}



### Analyze RD Data

path <- "/data/joy/BBL/studies/pnc/n2416_dataFreeze/neuroimaging/dti/n2416_DTI64/n2416_DTI_JHUTractRD.csv"
indata <- read.csv(path)
indata <- combineROI(indata)


data.analysis <- merge(demo, indata, by=c("bblid","scanid"))

covariates <- "~ s(AgeAtScan, k=4) + sex + race2+ gaf_z + dti64Tsnr "

gam_pval.rd = as.data.frame(matrix(NA, nrow=146, ncol=8))

for (i in 162:172) {
  
  var.name <- names(data.analysis)[i]
  m <- as.formula(paste(var.name, covariates,sep=""))
  gam1 <- gamm4(m, data=data.analysis, REML=T, random= ~ (1|bblid))
  
  n <- i 
  
  gam_pval.rd[n, 1] <- var.name
  gam_pval.rd[n, 2:(1 + dim(summary(gam1$gam)$p.table)[1])] <- summary(gam1$gam)$p.table[,4]
  gam_pval.rd[n, (2 + dim(summary(gam1$gam)$p.table)[1]):(1 
                                                          + dim(summary(gam1$gam)$p.table)[1] 
                                                          + dim(summary(gam1$gam)$s.table)[1]) ] <- summary(gam1$gam)$s.table[,4]
  
}


gam_pval.rd <- gam_pval.rd[!is.na(gam_pval.rd[,1]), ]
names(gam_pval.rd)[1] <- "ROI"
names(gam_pval.rd)[ 2:(1 + dim(summary(gam1$gam)$p.table)[1])] <- rownames(summary(gam1$gam)$p.table)
names(gam_pval.rd)[ (2 + dim(summary(gam1$gam)$p.table)[1]):(1 
                                                             + dim(summary(gam1$gam)$p.table)[1] 
                                                             + dim(summary(gam1$gam)$s.table)[1]) ] <- rownames(summary(gam1$gam)$s.table)

gam_corr.rd <- gam_pval.rd

for (i in 2:dim(gam_pval.rd)[2]) {
  gam_corr.rd[1:10,i] <- p.adjust(gam_pval.rd[1:10,i], method="bonferroni")
}



