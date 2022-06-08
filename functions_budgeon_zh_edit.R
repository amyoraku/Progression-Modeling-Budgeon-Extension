## START OF REG.FN()
## Function creates regression means and slopes w/ linear model or linear mixed model
## Outputs a DF with subject-specific fitted means and slopes
reg.fn <- function(dat, ID.uniq, reg.method,rm.slopes){
  ##Creating a matrix to include the slopes and means for LM and LMM
  dat2 <- matrix(nrow=length(ID.uniq),ncol=2)
  colnames(dat2) <- c("Slope","Mean")
  
  if(reg.method=="LM"){
    for (i in 1:length(ID.uniq)){
      fm1 <- lm(dat$response[dat$RID==ID.uniq[i]]~ dat$TIME[dat$RID==ID.uniq[i]])
      ##slope for each individual
      dat2[i,1] <- coef(fm1)[2]
      ##fitted mean for each individual
      dat2[i,2] <- sum(range(fitted(fm1)))/2
    }
  }
  else {
    if(reg.method=="LMM"){
      ## The below appears to be unnecessary code
      lmm_mid <- matrix(nrow=length(ID.uniq),ncol=2)
      lmm_mid[,2] <- ID.uniq
      for (i in 1:length(ID.uniq)){
        lmm_mid[i,1] <- sum(range(dat$TIME[dat$RID==ID.uniq[i]]))/2				
      }
      
      ## switched optimizer because it fails to converge otherwise
      lmm_mod <- lmer(response~TIME+(TIME|RID),data=dat,control=lmerControl(optimizer="bobyqa"))
      ##slope for each individual
      dat2[,1] <- coef(lmm_mod)$RID[,2]
      ##fitted mean for each individual
      dat2[,2] <- coef(lmm_mod)$RID[,1] + coef(lmm_mod)$RID[,2]*lmm_mid[,1]
    }
    else{
      #print("Please specify method")
    }
  }
  dat2 <- na.omit(as.data.frame(dat2))
  if(rm.slopes=="Yes"){
    dat2.use <- subset(dat2, !Slope<0)
  }
  else {
    dat2.use <- dat2
  }
}
#END OF REG.FN() HERE


#A FUNCTION USED IN THE NEXT FUNCTION TO CONDUCT REGRESSION OF SLOPE AGAINST MEAN
## START OF SM.FIT()
sm.fit <- function(data, p){
  y<-data[,1]
  x<-data[,2]
  #linear model fitting returns the fitting linear model
  ms<-paste("y~", paste(paste("I(x^", 1:p,")"), collapse="+"))
  lin.mod<-lm(formula(ms))
  return(lin.mod)
}
#END OF SM.FIT()

# START OF REG.MOD2()
# FUNCTION STARTS HERE TO OBTAIN MODEL COEFFICIENTS (CUBIC/QUINTIC)
reg.mod2 <- function(data,p.degree){
  model.sm <-sm.fit(data,p.degree)
  
  
  coef.mod <- coef(model.sm)
  list(
    coef.mod=coef.mod
  )
}
#END OF REG.MOD2()

## NNeg2Pol is the implementation of the sum-of-squares polynomial parameterization described on page 5, fig. 8
## start.mod is the beginning of the interval (a)
## finish.mod is the end of the interval (b)

## START OF NNeg2Pol()
NNeg2Pol <- function(par, start.mod,finish.mod){
  n <- length(par)
  q <- n/2
  p1 <- par[1:q]
  p2 <- par[q+1:q]

## below is squaring the individual polynomial components for use in the below formula  
  p1sq <- convolve(p1, rev(p1), type="o")
  p2sq <- convolve(p2, rev(p2), type="o")
  
## below is equivalent to p(x) = (x − a)*(p1(x))^2 + (b − x)*(p2(x))^2
  convolve(p1sq, c(1, -start.mod), type="o") + convolve(p2sq, c(-1, finish.mod), type="o")
}
## END OF NNeg2Pol()

## START OF makeRSS()
# Wraps fitting steps in a functional closure, used for optimization step in poly.mod
makeRSS <- function(x, y, start.mod, finish.mod){
  
  force(x)
  force(y)
  force(start.mod)
  force(finish.mod)
  
  function(par){
    beta <- NNeg2Pol(par, start.mod=start.mod, finish.mod=finish.mod)
    fit <- evalPol(x, beta)
    sum((y-fit)^2)
  }
}
## END OF makeRSS()

## Function fits polynomial with optimized sum-of-squares parameterization
##START OF poly.mod()
poly.mod <- function(p.degree, coef.mod, start.mod, finish.mod, x, y, poly.type){
  if(poly.type=="unc"){
    root <- sort(Re(polyroot(coef.mod)))
    b.degree <- coef.mod[p.degree+1]
  }
  else {
    if(poly.type=="nn"){
      par <- c(0,0,0,0)
      nngcff <- NNeg2Pol(par, start.mod, finish.mod) ## I don't think this line serves a purpose anymore
      RSS <- makeRSS(x, y, start.mod, finish.mod)
      
      ## optimization routine for parameter values
      res <- optim(par, fn=RSS,control=list(maxit=1000))
      while(res$convergence==1) {
        res <- optim(res$par, fn=RSS,control=list(maxit=1000))}

      if(p.degree==3){
        ## produces non-negative sum-of-squares parameterization for optimized values
        bb <- NNeg2Pol(res$par, start.mod, finish.mod)
        coef.mod <- bb
      }
      else{
        if(p.degree==5){
          ## optimization routine for parameter values
          par <- c(res$par[1:2],0,res$par[3:4],0)
          res <- optim(par, fn=RSS, control=list(maxit=1000))
          while(res$convergence==1) {
            res <- optim(res$par, fn=RSS, control=list(maxit=1000))}

          ## produces non-negative sum-of-squares parameterization for optimized values
          bb <- NNeg2Pol(res$par, start.mod, finish.mod)
          coef.mod <- bb
        }
      }
      
      root <- sort(Re(polyroot(coef.mod)))
      b.degree <- coef.mod[p.degree+1]
    }
    else {
      print("Specify poly type to use")
    }
  }
  list(x=x, y=y, root=root,b.degree= b.degree,coef.mod=coef.mod, start.mod=start.mod, finish.mod=finish.mod,res=res)
}
#END OF poly.mod()
## LEFT WITH THE ROOTS OF THE MODEL FOR EITHER UNCONSTRAINED OR CONSTRAINED MODEL


#Produces integration of reciprocal of fitted polynomial from poly.mod
#START OF int.fun()
int.fn <- function(xgr, root, b.degree, numerical.int,coef.mod,pred.in){
  #numerical.int == "F" only works when b.degree==3
  if(numerical.int=="F"){
    if(b.degree < 0) stop ("b.degree negative rethink")
    A <- 1/((root[2]-root[1])*(root[3]-root[1])* b.degree)
    B <- 1/((root[1]-root[2])*(root[3]-root[2])* b.degree)
    C <- 1/((root[1]-root[3])*(root[2]-root[3])* b.degree)
    y.int <- A*log(xgr-root[1])+B*log(root[2]-xgr)+C*log(root[3]-xgr)+12
  }
  else {
    if(numerical.int=="T"){
      y.int <-vector(length=length(xgr))
      fn.num.int <- function(x) {1/evalPol(x, coef.mod)}
      for (i in seq_along(xgr)){
        y.int[i] <- integrate(fn.num.int,lower=xgr[1],upper=xgr[i])$value
      }
    }
    else {
    }
  }
  list(y.int=y.int)
}
#END OF int.fn()






