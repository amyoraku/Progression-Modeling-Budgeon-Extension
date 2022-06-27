library(lme4)
library(ggplot2)
library(truncnorm)
library(nlme)
library(MonoPoly)
library(reshape2)
library(minqa)
library(plyr)
library(dplyr)

## Run functions_budgeon_zh_edit.R before this

adni_merge<-ADNIMERGE::adnimerge

# ZH data not provided, but I make the below assumptions about equiv variables:
# ZH TIME in CB dataset is equivalent to Time_Since_Baseline in AM dataset
# ZH ID is coded as RID in CB dataset and is equivalent to ID in AM dataset
# ZH response is manually coded shortly below in CB but is manually assigned as SUVR value

amy.data.all <- read.csv("~/Downloads/UCB_CL_no_cu_out_1.06.22.csv")

# Joined AM's data with ADNIMERGE for baseline DX
# Renamed AM's "ID" var to "RID" for compatibility with CB code/ADNI data
# dropped observations where Use.Base.Rate == "no" - this carries on to modeling
# did some type changes to deal with the issues with labelled columns in ADNIMERGE

amy.data.bridge <- dplyr::left_join(amy.data.all %>% dplyr::filter(Use.Base.Rate=="yes") %>% 
                                      dplyr::rename(RID=ID),adni_merge %>% dplyr::select(RID,DX.bl) %>% 
                                      dplyr::rename(DX_bl=DX.bl) %>% dplyr::transmute(RID=as.numeric(RID),DX_bl=as.factor(DX_bl)),by="RID")

## Two different definitions of DX at baseline - DX at BL from ADNIMERGE and DX at baseline observation
## Think DX at baseline observation (time=0, coded as DX_alt) here is correct but wanted to check

amy.data.bridge <- amy.data.bridge %>% dplyr::mutate(.,DX_bl=with(.,dplyr::case_when(
  (DX_bl %in% c("EMCI","LMCI")) ~ "MCI",
  (DX_bl %in% c("SMC","CN")) ~ "CN",
  (DX_bl == "AD") ~ "AD")
))

##DX_alt is correct
amy.data.bridge <- amy.data.bridge %>% dplyr::mutate(.,DX_alt=with(.,dplyr::case_when(
  (Dx %in% c("MCI","MCI_neg")) ~ "MCI",
  (Dx %in% c("CN_neg","CN")) ~ "CN",
  (Dx %in% c("AD","AD_neg") ~ "AD")
)))

# cut_comp categorizes subjects as Amyloid pos/neg
# not 100% sure about this - inferring from the tapply function below where medians are calculated for diff DX/positivity status combos
amy.data.bridge$cut_comp <- amy.data.bridge$SUMMARYSUVR_WHOLECEREBNORM_1.11CUTOFF

amy.data.use <- subset(amy.data.bridge, Use.Modeling == "yes")
amy.data.use <- amy.data.use %>% dplyr::mutate(RID=as.factor(RID)) 

# only 3+ time points
timept_3 <- subset(amy.data.use, Timepoint == 2)
keep <- timept_3$RID
amy.data.final <- subset(amy.data.use, RID %in% keep)

## changing name for compatibility with rest of code
data.use <- amy.data.final %>% dplyr::rename(TIME=Time_Since_Baseline)


ROI_modelling <- function(data, reg.method, p.degree, poly.type, optim, rm.slopes, n.bs, error, seed.num, out.length, spec_ROI){
  data.use$response <- data.use[ , spec_ROI]
  # sort the dataset by RID, then time
  data.use <- data.use[order(data.use$RID, data.use$TIME),] ## sorts dataset by RID, then time
  #CALCUALTING MEDIAN AND MEAN VALUES TO ADD ON PLOT
  meanSUVR <- tapply(data.use$response[data.use$TIME==0], list(data.use$DX_alt[data.use$TIME==0],data.use$cut_comp[data.use$TIME==0]),mean,na.rm=TRUE)
  medSUVR <- tapply(data.use$response[data.use$TIME==0], list(data.use$DX_alt[data.use$TIME==0],data.use$cut_comp[data.use$TIME==0]),median)
  SUVR_HCneg <- medSUVR[2,1]
  ## cut_SUVR <- data.use$cut_comp_value[1] #THIS IS THE VALUE OF 0.79 in PAPER - represented clinical cut off
  ## temporarily ignoring the cut_SUVR because I'm not sure we have a regional equivalent in Alison's case and keeping it with an arbitrary value causes downstream problems in the code
  ## one downstream problem: the cut_SUVR value is included in the vector used for numerical integration, which causes problems if it's outside of the bounds of the rest of the integral
  ## cut_SUVR <- 79
  SUVR_ADpos <- medSUVR[1,2]
  # collapsed R_ADNI_bs_alt and R_source_alt into one file - options here apply to original method and bootstrap
  set.seed(seed.num)
  reg.out <- reg.fn(data.use,
                    ID.uniq=unique(data.use$RID),
                    reg.method=reg.method,
                    rm.slopes=rm.slopes)
  xgr.seq <- unique(sort(c(SUVR_ADpos, SUVR_HCneg,seq(min(reg.out$Mean),max(reg.out$Mean),length=out.length))))
  #RUNNING POLYNOMIAL
  reg.mod2.out <- reg.mod2(data=reg.out, p.degree=p.degree)
  ## producing start/end points for interval in integration w/ a small adjustment for error
  if(min(reg.out$Mean)==0) {start.mod.F <- min(reg.out$Mean)-(error)
  } else if(min(reg.out$Mean)<0) {start.mod.F <- min(reg.out$Mean)*(1+error)
  } else {start.mod.F <- min(reg.out$Mean)*(1-error)}
  if(max(reg.out$Mean)==0) {finish.mod.F <- max(reg.out$Mean)+(error)
  } else if(max(reg.out$Mean)<0) {finish.mod.F <- max(reg.out$Mean)*(1-error)
  } else {finish.mod.F <- max(reg.out$Mean)*(1+error)}
  #}
  #RUNNING NON-NEGATIVE POLYNOMIAL
  ## ZH see Functions for more detail - function relies on NNegPol2(), makeRSS()
  poly.mod.out <- poly.mod(p.degree=p.degree,
                           coef.mod=reg.mod2.out$coef.mod,
                           start.mod=start.mod.F,
                           finish.mod=finish.mod.F,
                           x= reg.out$Mean,
                           y= reg.out$Slope,
                           poly.type=poly.type)
  int.out <- int.fn(xgr= xgr.seq,
                    root= poly.mod.out$root,
                    b.degree= poly.mod.out$b.degree,
                    numerical.int="T",
                    coef.mod= poly.mod.out$coef.mod)
  CU_SUVRADpos.F <- int.out$y.int[which(xgr.seq==SUVR_ADpos)]
  CU_SUVRHCneg.F <- int.out$y.int[which(xgr.seq==SUVR_HCneg)]
  ## CU_SUVRcut.F <- int.out$y.int[which(xgr.seq== cut_SUVR)]
  ## see above note re: cut_SUVR
  #set.seed(seed.num)
  out.y.int.bs <- matrix(nrow=length(xgr.seq), ncol=n.bs)
  #FOR CUT OFFS AND CIs
  cut_offs <- matrix(nrow=n.bs, ncol=3)
  colnames(cut_offs) <- c("SUVRADpos","SUVRHCneg","SUVRcut")
  # bootstrap section below - just runs same process as above with resampling
  for (i in 1:n.bs){
    # resampling for bootstrap
    ID.sample <- sample(unique(data.use$RID), length(unique(data.use$RID)), replace = TRUE)
    ID <- seq(1,length(ID.sample),1)
    ID.use <- cbind(ID.sample, ID)
    colnames(ID.use) <- c("RID","ID")
    ID.use.test<-as.data.frame(ID.use)
    ## Fixed bootstrap problem noted in original file
    data.bs.test<-dplyr::left_join(ID.use.test,data.use %>% dplyr::mutate(RID=as.numeric(RID)),by="RID")
    data.bs.test <- data.bs.test[order(data.bs.test$ID, data.bs.test$TIME),]
    data.bs.test <- dplyr::rename(data.bs.test, OLD_RID=RID, RID=ID)
    reg.out.BS <- reg.fn(dat=data.bs.test,
                         ID.uniq=unique(data.bs.test$RID),
                         reg.method=reg.method,
                         rm.slopes=rm.slopes)
    #RUNNING POLYNOMIAL
    reg.mod2.out.BS <- reg.mod2(data=reg.out.BS, p.degree=p.degree)
    if(min(reg.out.BS$Mean)==0) {start.mod.BS <- min(reg.out.BS$Mean)-(error)
    } else if(min(reg.out.BS$Mean)<0) {start.mod.BS <- min(reg.out.BS$Mean)*(1+error)
    } else {start.mod.BS <- min(reg.out.BS$Mean)*(1-error)}
    if(max(reg.out.BS$Mean)==0) {finish.mod.BS <- max(reg.out.BS$Mean)+(error)
    } else if(max(reg.out.BS$Mean)<0) {finish.mod.BS <- max(reg.out.BS$Mean)*(1-error)
    } else {finish.mod.BS <- max(reg.out.BS$Mean)*(1+error)}
    #RUNNING NON NEGATIVE POLYNOMIAL
    poly.mod.out.BS <- poly.mod(p.degree=p.degree,
                                coef.mod=reg.mod2.out.BS$coef.mod,
                                #start.mod=min(start.mod.F, start.mod.BS),
                                #finish.mod=max(finish.mod.F,finish.mod.BS),
                                start.mod=start.mod.F,
                                finish.mod=finish.mod.F,
                                x= reg.out.BS$Mean,
                                y= reg.out.BS$Slope,
                                poly.type=poly.type)
    int.out.BS <- int.fn(xgr= xgr.seq,
                         root= poly.mod.out.BS$root,
                         b.degree= poly.mod.out.BS$b.degree,
                         numerical.int="T",
                         coef.mod= poly.mod.out.BS$coef.mod)
    cut_offs[i,1] <- int.out.BS$y.int[which(xgr.seq==SUVR_ADpos)]
    cut_offs[i,2] <- int.out.BS$y.int[which(xgr.seq==SUVR_HCneg)]
    ##  cut_offs[i,3] <- int.out.BS$y.int[which(xgr.seq== cut_SUVR)]
    out.y.int.bs[,i] <- int.out.BS$y.int
    stop_code<-i
    print(i)
  }
  #OBTAINING QUANTILES FOR BOOTSTRAP DATA
  perc.5.BS <- apply(out.y.int.bs, 1, quantile, probs = 0.05,  na.rm = TRUE)
  perc.2.5.BS <- apply(out.y.int.bs, 1, quantile, probs = 0.025,  na.rm = TRUE)
  perc.50.BS <- apply(out.y.int.bs, 1, quantile, probs = 0.5,  na.rm = TRUE)
  perc.95.BS <- apply(out.y.int.bs, 1, quantile, probs = 0.95,  na.rm = TRUE)
  perc.97.5.BS <- apply(out.y.int.bs, 1, quantile, probs = 0.975,  na.rm = TRUE)
  #COMBINING WITH X SEQUENCE
  data.conf <- cbind(xgr.seq, perc.5.BS,perc.2.5.BS,perc.50.BS,perc.95.BS,perc.97.5.BS)
  ## Plotting wit ggplot
  plot_data <- data.frame(xgr.seq, int.out$y.int) #xgr.seq is on the y axis, int.out on x axis
  conf_1_cols <- c(1, 3)
  conf_2_cols <- c(1, 6)
  conf_1_data <- as.data.frame(data.conf[, conf_1_cols])
  conf_2_data <- as.data.frame(data.conf[, conf_2_cols])
  test_plot <- ggplot(plot_data, aes(int.out.y.int, xgr.seq)) +
    geom_line() +
    geom_line(aes(perc.2.5.BS, xgr.seq), data = conf_1_data, linetype = "dashed") +
    geom_line(aes(perc.97.5.BS, xgr.seq), data = conf_2_data, linetype = "dashed") +
    xlab("Disease Progression (years)") +
    ylab("") +
    theme_classic()
  test_plot
}

#ROI_modelling(data=data.use, reg.method="LMM", p.degree=3, poly.type="nn", optim="ALL", rm.slopes="Yes", n.bs=10, error=0.1, seed.num=200, out.length=3001, spec_ROI = "CTX_RH_CUNEUS_CL")
