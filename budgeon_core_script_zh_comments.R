#THIS SCRIPT RUNS THE TWO FUNCTIONS T DO THE FOUR STEP APPROACH ON ACTUAL DATA

library(lme4)
library(ggplot2)
library(truncnorm)
library(nlme)
library(MonoPoly)
library(reshape2)
library(minqa)
library(plyr)
#For parallel
library(doRNG)
library(foreach)
library(doParallel)

#THE FOLLOWING FILE READS IN THE ADNI DATA AND COMBINES WITH DX - all done already and now variables were created in SAS.
rm(list=ls(all=TRUE))

#THIS IS MY DATA SET - HAS TIME, ID and RESPONSE VARIABLE INCLUDED
data.complete <- read.csv('ADNI_combined_SAS_280416.csv',header=TRUE)

#MINIMUM NUMBER OF DATA POINTS TO USE (FOR LONGITUDINAL DATA - ENSURE ALL HAVE AT LEAST 3 - MAINLY FOR MIXED MODELS)
min.pts <- 3
  
data.use <- data.complete[data.complete$Frequency >= min.pts,] #SUBSETTING THE DATA TO HAVE ONLY THOSE WITH MINIMUM NUMBER OF POINTS

data.use <- data.use[order(data.use$RID, data.use$TIME),]

data.use$response <- data.use$SUMMARYSUVR_COMPOSITE_REFNORM #NAME OF RESPONSE VARIABLE

#CALCUALTING MEDIAN AND MEAN VALUES TO ADD ON PLOT
meanSUVR <- tapply(data.use$response[data.use$TIME==0], list(data.use$DX_bl[data.use$TIME==0],data.use$cut_comp[data.use$TIME==0]),mean)
medSUVR <- tapply(data.use$response[data.use$TIME==0], list(data.use$DX_bl[data.use$TIME==0],data.use$cut_comp[data.use$TIME==0]),median)

SUVR_HCneg <- medSUVR[2,1]
cut_SUVR <- data.use$cut_comp_value[1] #THIS IS THE VALUE OF 0.79 in PAPER - represented clinical cut off
SUVR_ADpos <- medSUVR[1,2]

#FITTED CURVE USING ALL DATA - CHOOSE OPTIONS THAT YOU WANT - SEE FUNCTION
reg.method <- "LMM"
p.degree <- 3
poly.type <- "nn" #DONT NEED TO CHANGE
optim <- "ALL" #DONT NEED TO CHANGE
rm.slopes <- "Yes"	
n.bs<-1000
error <- 0.1
seed.num <- 200
out.length <- 3001

#RUNS THE BOOTSTRAPPING AND THE OTHER FUNCTION
source("/Users/charleybudgeon/Dropbox/PhD/CODE FOR DISTRIBUTION/R_ADNI_bs.R")

#OBTAINING QUANTILES FOR BOOTSTRAP DATA
perc.5.BS <- apply(out.y.int.bs, 1, quantile, probs = 0.05,  na.rm = TRUE)
perc.2.5.BS <- apply(out.y.int.bs, 1, quantile, probs = 0.025,  na.rm = TRUE)
perc.50.BS <- apply(out.y.int.bs, 1, quantile, probs = 0.5,  na.rm = TRUE)
perc.95.BS <- apply(out.y.int.bs, 1, quantile, probs = 0.95,  na.rm = TRUE)
perc.97.5.BS <- apply(out.y.int.bs, 1, quantile, probs = 0.975,  na.rm = TRUE)

#COMBINING WITH X SEQUENCE 
data.conf <- cbind(xgr.seq, perc.5.BS,perc.2.5.BS,perc.50.BS,perc.95.BS,perc.97.5.BS)

#THIS IS WHAT I SAVED THE INDIVIDUAL OPTIONS AS - AFTER CHANGING THE OPTIONS ABOVE
save(list = ls(all.names = TRUE),file=paste0("/Users/charleybudgeon/Dropbox/PhD/CODE FOR DISTRIBUTION/ADNI_BS",reg.method,p.degree,rm.slopes,min.pts,".RData"))


###PLOTTING THE PROGRESSION CURVE

rm(list=ls(all=TRUE))

reg.method <- "LMM"
p.degree <- 3
rm.slopes <- "Yes"
n.bs <- 1000
min.pts <- 3

#THIS IS WHAT THE DATAFILES ARE SAVED AS
load(file=paste0("/Users/charleybudgeon/Dropbox/PhD/CODE FOR DISTRIBUTION/ADNI_BS",reg.method,p.degree,rm.slopes,min.pts,".RData"))

#Now to work out time taken between different points
#FROM MEAN SUVR HC negtive (0.7099) to SUVR cut off (0.79)
CU_SUVRHCneg.F
CU_SUVRcut.F
CU_SUVRADpos.F

#FITTED DATA
T3_T2_F <- CU_SUVRADpos.F - CU_SUVRcut.F
T2_T1_F <- CU_SUVRcut.F - CU_SUVRHCneg.F

#BOOTSTRAP
cut_offs <- as.data.frame(cut_offs)

T1_HCneg <- quantile(cut_offs$SUVRHCneg, c(0.025,0.5,0.975))
T2_Cut <- quantile(cut_offs$SUVRcut, c(0.025,0.5, 0.975))
T3_ADpos <- quantile(cut_offs$SUVRADpos, c(0.025,0.5,0.975))

T3_T2_BS <- cut_offs$SUVRADpos - cut_offs$SUVRcut
T2_T1_BS <- cut_offs$SUVRcut - cut_offs$SUVRHCneg

CI_T3_T2 <- quantile(T3_T2_BS, c(0.025,0.975))
CI_T2_T1 <- quantile(T2_T1_BS, c(0.025,0.975))


#FITTED
rate23_F <- (SUVR_ADpos -cut_SUVR)/(T3_T2_F)

#BS
rate23_BS <- (SUVR_ADpos  - cut_SUVR )/(T3_T2_BS)
CI_23rate <- quantile(rate23_BS, c(0.025,0.5,0.975))

#ADDING BOXPLOTS
#IMPORT ADNI DATA
data.comb3 <- read.csv('/Users/charleybudgeon/Dropbox/PhD/CODE FOR DISTRIBUTION/ADNI_combined_SAS_280416.csv',header=TRUE)

data.use <- data.comb3[data.comb3$Frequency >= min.pts,] #SUBSETTING THE DATA TO HAVE ONLY THOSE WITH MINIMUM NUMBER OF POINTS

data.use $SUVR_scale <- (data.use $SUMMARYSUVR_COMPOSITE_REFNORM-min(data.use $SUMMARYSUVR_COMPOSITE_REFNORM))/(max(data.use $SUMMARYSUVR_COMPOSITE_REFNORM)-min(data.use $SUMMARYSUVR_COMPOSITE_REFNORM)) 

data.use <- data.use[order(data.use$RID, data.use$TIME),]

meanSUVR <- tapply(data.use$SUMMARYSUVR_COMPOSITE_REFNORM[data.use$TIME==0], list(data.use$DX_bl[data.use$TIME==0]),mean)

data.use$DX_bl <- factor(data.use$DX_bl,
                         levels = c("HC", "MCI", "AD"))

par(mgp=c(2.5,1,0),mar=c(3,5,2,2))
boxplot(data.use$SUMMARYSUVR_COMPOSITE_REFNORM[data.use$TIME==0]~data.use$DX_bl[data.use$TIME==0],
        col=c("blue","dark green","red"),
        ylab=expression(paste("Neocortical Amyloid - ",beta," Burden (SUVR)")),
        cex.lab=2,cex.axis=2,pch=16,cex=0.8)

meanSUVR <- tapply(data.use$SUMMARYSUVR_COMPOSITE_REFNORM[data.use$TIME==0], list(data.use$DX_bl[data.use$TIME==0],data.use$cut_comp[data.use$TIME==0]),mean)
medSUVR <- tapply(data.use$SUMMARYSUVR_COMPOSITE_REFNORM[data.use$TIME==0], list(data.use$DX_bl[data.use$TIME==0],data.use$cut_comp[data.use$TIME==0]),median)
minSUVR <- tapply(data.use$SUMMARYSUVR_COMPOSITE_REFNORM[data.use$TIME==0], list(data.use$DX_bl[data.use$TIME==0],data.use$cut_comp[data.use$TIME==0]),min)
maxSUVR <- tapply(data.use$SUMMARYSUVR_COMPOSITE_REFNORM[data.use$TIME==0], list(data.use$DX_bl[data.use$TIME==0],data.use$cut_comp[data.use$TIME==0]),max)

#SUVR_HCneg <- meanSUVR[2,1]
SUVR_HCneg <- medSUVR[2,1]
cut_SUVR <- data.use$cut_comp_value[1]
#NEED TO USE MEDIAN FOR UPDATE OF PAPER
#SUVR_ADpos <- meanSUVR[1,2]
SUVR_ADpos <- medSUVR[1,2]

dat.ADNI <- data.use[data.use$TIME==0,]

dat.ADNI$DX_bl <- ordered(dat.ADNI$DX_bl, levels = c("HC", "MCI", "AD"))

par(fig=c(0.15,1,0,1),mar=c(4,2,2,2))
plot(xgr.seq~int.out$y.int,type="l",
     xlim=c(0,max(int.out$y.int)+10),
     ylim=c(min(xgr.seq), max(dat.ADNI$SUMMARYSUVR_COMPOSITE_REFNORM)-0.001),
     ylab="",xlab="Disease Progression (years)",
     lwd=2,col="black",main="",mgp=c(2,1,0),cex.lab=1.3,cex.axis=1.2)
lines(data.conf[,1]~data.conf[,3],lty=2)
lines(data.conf[,1]~data.conf[,6],lty=2)

abline(h=SUVR_HCneg,lty=3,col="blue")
abline(h=cut_SUVR,lty=3)
abline(h=SUVR_ADpos,lty=1,col="red")

segments(CU_SUVRHCneg.F,0, CU_SUVRHCneg.F, SUVR_HCneg,lty=3)
segments(CU_SUVRcut.F,0, CU_SUVRcut.F, cut_SUVR,lty=3)
segments(CU_SUVRADpos.F,0, CU_SUVRADpos.F, SUVR_ADpos,lty=3)

text(x=max(int.out$y.int)+6.3, y=SUVR_ADpos-0.03,paste0("Median SUVR AD+\n(",round(SUVR_ADpos,2),")"),cex=1.2)
text(x=max(int.out$y.int)+6.3, y=cut_SUVR-0.02,paste0("SUVR ",round(cut_SUVR,2)),cex=1.2)
text(x=max(int.out$y.int)+6.3, y=SUVR_HCneg-0.03,paste0("Median SUVR HC-\n(",round(SUVR_HCneg,2),")"),cex=1.2)

#EXCLUDES CI
text(x= (CU_SUVRcut.F - CU_SUVRHCneg.F)/2+ CU_SUVRHCneg.F, y= SUVR_HCneg-0.03, 
     paste0(round(T2_T1_F,2)," years"),cex=1.2)
arrows(CU_SUVRHCneg.F+0.1, SUVR_HCneg-0.045, CU_SUVRcut.F-0.1, SUVR_HCneg-0.045,length=0.1)
arrows(CU_SUVRcut.F-0.1, SUVR_HCneg-0.045, CU_SUVRHCneg.F+0.1, SUVR_HCneg-0.045,length=0.1)

#EXCUDES CI
text(x=(CU_SUVRADpos.F - CU_SUVRcut.F)/2+ CU_SUVRcut.F, y=cut_SUVR-0.03, paste0(round(T3_T2_F,2)," years"),cex=1.2)
arrows(CU_SUVRcut.F+0.1, cut_SUVR-0.045, CU_SUVRADpos.F-0.1, cut_SUVR-0.045,length=0.1)
arrows(CU_SUVRADpos.F-0.1, cut_SUVR-0.045, CU_SUVRcut.F+0.1, cut_SUVR-0.045,length=0.1)

#EXCLUDES CI
text(x= (CU_SUVRcut.F - CU_SUVRHCneg.F)/2+ CU_SUVRHCneg.F+2,y= SUVR_ADpos-0.03, 
     paste0(round(rate23_F,3)," SUVR/year"),cex=1.2)
arrows(CU_SUVRHCneg.F, SUVR_ADpos-0.01, CU_SUVRHCneg.F, cut_SUVR+0.01,length=0.1)
arrows(CU_SUVRHCneg.F, cut_SUVR+0.01, CU_SUVRHCneg.F, SUVR_ADpos-0.01,length=0.1)

par(fig=c(0,0.15,0,1), mar=c(4,2,2,1), new=TRUE)

boxplot(dat.ADNI$SUMMARYSUVR_COMPOSITE_REFNORM[dat.ADNI$cut_comp=='Negative']~dat.ADNI$DX_bl[dat.ADNI$cut_comp=='Negative'],
        ylim=c(min(xgr.seq),max(dat.ADNI$SUMMARYSUVR_COMPOSITE_REFNORM)-0.001), 
        col=c("blue","dark green","red"),outline=FALSE,axes=FALSE)
axis(1,at=c(1,2,3),labels=c("HC-","MCI-","AD-"),tick=FALSE,cex.axis=1)
mtext(expression(paste("Neocortical Amyloid - ",beta," Burden (SUVR)")),side=2,line=0.5,cex=1.3)

boxplot(dat.ADNI$SUMMARYSUVR_COMPOSITE_REFNORM[dat.ADNI$cut_comp=='Positive']~dat.ADNI$DX_bl[dat.ADNI$cut_comp=='Positive'],add=TRUE,col=c("blue","dark green","red"),outline=FALSE,axes=FALSE)

text(1,max(dat.ADNI$SUMMARYSUVR_COMPOSITE_REFNORM)+0.012,"HC+",cex=1)
text(2,max(dat.ADNI$SUMMARYSUVR_COMPOSITE_REFNORM)+0.012,"MCI+",cex=1)
text(3,max(dat.ADNI$SUMMARYSUVR_COMPOSITE_REFNORM)+0.012,"AD+",cex=1)














