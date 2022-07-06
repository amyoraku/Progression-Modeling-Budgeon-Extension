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

# AM remove the line below b/c don't need any data from the ADNIMERGE dataset
# adni_merge<-ADNIMERGE::adnimerge

# ZH data not provided, but I make the below assumptions about equiv variables:
# ZH TIME in CB dataset is equivalent to Time_Since_Baseline in AM dataset
# ZH ID is coded as RID in CB dataset and is equivalent to ID in AM dataset
# ZH response is manually coded shortly below in CB but is manually assigned as SUVR value

amy.data.all <- read.csv("~/Downloads/UCB_adj_SUVR_04.26.22.csv") # using database with all SUVR cases

# AM For this specific dataset, Use.Modeling is hard coded to select all cases that are amyloid positive and have 3+ time points
amy.data.bridge <- subset(amy.data.all, Use.Modeling == "yes")

# previous code here was redundant of above b/c filtering by Use.Modeling
amy.data.use <- amy.data.bridge
amy.data.use <- amy.data.use %>% dplyr::rename(TIME=Time_Since_Baseline)

# AM removed code for Dx at baseline that combined X_neg cases into X (want to keep distinction)
## Two different definitions of DX at baseline - DX at BL from ADNIMERGE and DX at baseline observation
## Think DX at baseline observation (time=0, coded as DX_alt) here is correct but wanted to check

##DX_alt is correct
#amy.data.bridge <- amy.data.bridge %>% dplyr::mutate(.,DX_alt=with(.,dplyr::case_when(
  #(Dx %in% c("MCI","MCI_neg")) ~ "MCI",
  #(Dx %in% c("CN_neg","CN")) ~ "CN",
  #(Dx %in% c("AD","AD_neg") ~ "AD")
#)))

# this is the dataframe that will be used for calculating thresholds
all.threshold.data <- subset(amy.data.all, Use.Base.Rate == "yes")
# changing name for compatibility with rest of code
all.threshold.data <- all.threshold.data %>% dplyr::rename(TIME=Time_Since_Baseline)

# AM removed the below lines for now -- changed the code in the function for threshold calculation so that it does not rely on cut.comp
# cut_comp categorizes subjects as Amyloid pos/neg
# not 100% sure about this - inferring from the tapply function below where medians are calculated for diff DX/positivity status combos
# amy.data.bridge$cut_comp <- amy.data.bridge$SUMMARYSUVR_WHOLECEREBNORM_1.11CUTOFF

# AM redundant now b/c of Use.Modeling above
# amy.data.use <- subset(amy.data.bridge, Use.Modeling == "yes")
# amy.data.use <- amy.data.use %>% dplyr::mutate(RID=as.factor(RID)) 

# AM also redundant b/c Use.Modeling above
# only 3+ time points
#timept_3 <- subset(amy.data.use, Timepoint == 2)
#keep <- timept_3$RID
#amy.data.final <- subset(amy.data.use, RID %in% keep)

ROI_modelling <- function(data, reg.method, p.degree, poly.type, optim, rm.slopes, n.bs, error, seed.num, out.length, spec_ROI){
  data.use$response <- data.use[ , spec_ROI]
  #CALCUALTING MEDIAN AND MEAN VALUES TO ADD ON PLOT
  all.threshold.data$threshold <- all.threshold.data[, spec_ROI]
  CU <- median(all.threshold.data$threshold[all.threshold.data$Dx=="CN"])
  MCI <- median(all.threshold.data$threshold[all.threshold.data$Dx=="MCI"])
  AD <- median(all.threshold.data$threshold[all.threshold.data$Dx=="AD"])
  #meanSUVR <- tapply(data.use$response[data.use$TIME==0], list(data.use$DX_alt[data.use$TIME==0],data.use$cut_comp[data.use$TIME==0]),mean,na.rm=TRUE)
  #medSUVR <- tapply(data.use$response[data.use$TIME==0], list(data.use$DX_alt[data.use$TIME==0],data.use$cut_comp[data.use$TIME==0]),median)
  #SUVR_HCneg <- medSUVR[2,1]
  ## cut_SUVR <- data.use$cut_comp_value[1] #THIS IS THE VALUE OF 0.79 in PAPER - represented clinical cut off
  ## temporarily ignoring the cut_SUVR because I'm not sure we have a regional equivalent in Alison's case and keeping it with an arbitrary value causes downstream problems in the code
  ## one downstream problem: the cut_SUVR value is included in the vector used for numerical integration, which causes problems if it's outside of the bounds of the rest of the integral
  ## cut_SUVR <- 79
  #SUVR_ADpos <- medSUVR[1,2]
  
  # collapsed R_ADNI_bs_alt and R_source_alt into one file - options here apply to original method and bootstrap
  set.seed(seed.num)
  reg.out <- reg.fn(data.use,
                    ID.uniq=unique(data.use$RID),
                    reg.method=reg.method,
                    rm.slopes=rm.slopes)
  xgr.seq <- unique(sort(c(CU,MCI, AD,seq(min(reg.out$Mean),max(reg.out$Mean),length=out.length)))) # change SUVR_X to CU, MCI, and AD
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
  CU_SUVRCU.F <- int.out$y.int[which(xgr.seq==CU)] # AM change to SUVRCU
  CU_SUVRMCI.F <- int.out$y.int[which(xgr.seq==MCI)] # AM change to SUVRMCI
  CU_SUVRAD.F <- int.out$y.int[which(xgr.seq==AD)] # AM add for AD
  ## CU_SUVRcut.F <- int.out$y.int[which(xgr.seq== cut_SUVR)]
  ## see above note re: cut_SUVR
  #set.seed(seed.num)
  out.y.int.bs <- matrix(nrow=length(xgr.seq), ncol=n.bs)
  #FOR CUT OFFS AND CIs
  cut_offs <- matrix(nrow=n.bs, ncol=3)
  colnames(cut_offs) <- c("SUVR_CU","SUVR_MCI","SUVR_AD") # AM change names of columns
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
    cut_offs[i,1] <- int.out.BS$y.int[which(xgr.seq==SUVR_AD)] # AM change name of SUVR_
    cut_offs[i,2] <- int.out.BS$y.int[which(xgr.seq==SUVR_CU)] # AM change name of SUVR_
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
    geom_hline(yintercept=c(CU, MCI, AD), color=c('blue','green', 'red')) + # AM adding lines for thresholds
    xlim(0,65) + # AM specific limits for this dataset
    ylim(0.5, 1.3) + # AM specific limits for this dataset
    xlab("Disease Progression (years)") +
    ylab("SUVR") +
    labs(title = spec_ROI) + # AM adding ROI as title of plot
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14)) # AM increasing size of fonts on plot
  test_plot
}

### Run modeling

# save model to extract info
model <- ROI_modelling(data=amy.data.use, 
              reg.method="LMM", 
              p.degree=3, 
              poly.type="nn", 
              optim="ALL", 
              rm.slopes="No", # change this depending on whether you want to include negative slopes or not
              n.bs=1000, 
              error=0.1, 
              seed.num=200, 
              out.length=3001, 
              spec_ROI = "CTX_ROSTRALANTERIORCINGULATE_SUVR") # AM change spec_ROI for each run

# plot graph
model$plot_env$test_plot

# thresholds
model$plot_env$CU
model$plot_env$MCI
model$plot_env$AD

# time corresponding to cut off values
model$plot_env$CU_SUVRCU.F
model$plot_env$CU_SUVRMCI.F
model$plot_env$CU_SUVRAD.F

# save means and slopes for future analysis
save <- model$plot_env$reg.out
write.csv(save, "CTX_ROSTRALANTERIORCINGULATE_reg_output.csv")

## to extract just the slopes (if full function fails)
amy.data.use$response <- amy.data.use$CTX_ROSTRALANTERIORCINGULATE_SUVR
reg.out <- reg.fn(amy.data.use,
                  ID.uniq=unique(amy.data.use$RID),
                  reg.method="LMM",
                  rm.slopes="No")

write.csv(reg.out, "CTX_ROSTRALANTERIORCINGULATE_reg_output.csv")

