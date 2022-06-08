#CB THIS USES THE FUNCTION BUT YOU NEED TO HAVE SOURCED THE DATA FIRST - SOME OF THESE TERMS ARE UNIQUE TO THE VARIABLES YOU ARE USING"
# ZH need to have run the function code -  Fn_sample_IVEC_ADNI.R/budgeon_functions_zh_comments.R - referred hereafter as "Functions"
# ZH need to have imported data in a compatible way first - use code in  R_source.R/budgeon_core_script_zh_comments.R - referred hereafter as "Source"

source("/Users/charleybudgeon/Dropbox/PhD/CODE FOR DISTRIBUTION/Fn_sample_IVEC_ADNI.R")

## ZH seed.num set in Source
set.seed(seed.num)

## ZH primary regression function for deriving midpoints/slopes - see Functions
## ZH options set in Source
reg.out <- reg.fn(data.use,
                    ID.uniq=unique(data.use$RID),
                    reg.method=reg.method,
                    rm.slopes=rm.slopes)

xgr.seq <- unique(sort(c(SUVR_ADpos, SUVR_HCneg, cut_SUVR,seq(min(reg.out$Mean),max(reg.out$Mean),length=out.length))))

#CB RUNNING POLYNOMIAL
reg.mod2.out <- reg.mod2(data=reg.out, p.degree=p.degree)

## ZH producing start/end points for interval in integration w/ a small adjustment for error
if(min(reg.out$Mean)==0) {start.mod.F <- min(reg.out$Mean)-(error)
} else if(min(reg.out$Mean)<0) {start.mod.F <- min(reg.out$Mean)*(1+error)
} else {start.mod.F <- min(reg.out$Mean)*(1-error)}

if(max(reg.out$Mean)==0) {finish.mod.F <- max(reg.out$Mean)+(error)
} else if(max(reg.out$Mean)<0) {finish.mod.F <- max(reg.out$Mean)*(1-error)
} else {finish.mod.F <- max(reg.out$Mean)*(1+error)}
#}

#CB RUNNING NON NEGATIVE POLYNOMIAL
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

#CB NEED TO WORK OUT TIME FOR WHEN y=cut offs
#ZH below produces cutoffs for plotting
CU_SUVRADpos.F <- int.out$y.int[which(xgr.seq==SUVR_ADpos)]
CU_SUVRHCneg.F <- int.out$y.int[which(xgr.seq==SUVR_HCneg)]
CU_SUVRcut.F <- int.out$y.int[which(xgr.seq== cut_SUVR)]


#CB BOOTSTRAP PART OF FUNCTION USE SEQUENECE FROM ABOVE
#ZH the bootstrap below iterates over the process in the code above/described in the paper and stores the results in out.y.int.bs

#CB set.seed(seed.num)

out.y.int.bs <- matrix(nrow=length(xgr.seq), ncol=n.bs)

#CB FOR CUT OFFS AND CIs
cut_offs <- matrix(nrow=n.bs, ncol=3)
colnames(cut_offs) <- c("SUVRADpos","SUVRHCneg","SUVRcut")


for (i in 1:n.bs){
  
## ZH the sample step is written incorrectly, which produces unrealistic confidence intervals
  ID.sample <- sample(unique(data.use$RID), length(unique(data.use$RID)), replace = TRUE)
  ID <- seq(1,length(ID.sample),1)
  ID.use <- cbind(ID.sample, ID)
  colnames(ID.use) <- c("RID","ID")
  ## ZH the merge function does not keep all rows from the left-hand dataframe, which leaves you with one observation/RID instead of several
  data.bs <- merge(data.use, ID.use, by="RID")
  data.bs <- data.bs[order(data.bs$ID, data.bs$TIME),]
  data.bs <- rename(data.bs, c("RID"="OLD_RID", "ID"="RID"))
  
  reg.out.BS <- reg.fn(dat=data.bs,
                       ID.uniq=unique(data.bs$RID),
                       reg.method=reg.method,
                       rm.slopes=rm.slopes)
  
  #CB RUNNING POLYNOMIAL
  reg.mod2.out.BS <- reg.mod2(data=reg.out.BS, p.degree=p.degree)
  
  	if(min(reg.out.BS$Mean)==0) {start.mod.BS <- min(reg.out.BS$Mean)-(error)
	} else if(min(reg.out.BS$Mean)<0) {start.mod.BS <- min(reg.out.BS$Mean)*(1+error)
	} else {start.mod.BS <- min(reg.out.BS$Mean)*(1-error)}

	if(max(reg.out.BS$Mean)==0) {finish.mod.BS <- max(reg.out.BS$Mean)+(error)
	} else if(max(reg.out.BS$Mean)<0) {finish.mod.BS <- max(reg.out.BS$Mean)*(1-error)
	} else {finish.mod.BS <- max(reg.out.BS$Mean)*(1+error)}

  
  #CB RUNNING NON NEGATIVE POLYNOMIAL
  poly.mod.out.BS <- poly.mod(p.degree=p.degree,
                              coef.mod=reg.mod2.out.BS$coef.mod,
                              #start.mod=min(start.mod.F, start.mod.BS),
                              #finish.mod=max(finish.mod.F,finish.mod.BS),
                              start.mod=start.mod.F,
                              finish.mod=finish.mod.F,
                              x= reg.out.BS$Mean,
                              y= reg.out.BS$Slope,
                              poly.type=poly.type)
  
  #CB NEED TO LOOK INTO USING OWN DATA WHEN BOOTSTRAPPING
  
  int.out.BS <- int.fn(xgr= xgr.seq,
                       root= poly.mod.out.BS$root,
                       b.degree= poly.mod.out.BS$b.degree,
                       numerical.int="T",
                       coef.mod= poly.mod.out.BS$coef.mod)
  
  cut_offs[i,1] <- int.out.BS$y.int[which(xgr.seq==SUVR_ADpos)]
  cut_offs[i,2] <- int.out.BS$y.int[which(xgr.seq==SUVR_HCneg)]
  cut_offs[i,3] <- int.out.BS$y.int[which(xgr.seq== cut_SUVR)]
  
  
  out.y.int.bs[,i] <- int.out.BS$y.int
  
  print(i)
}


