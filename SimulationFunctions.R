## 5 parameter sigmoid curve
pl.5 <- function (a,b,c,f,g, tau) {
  b.star <- b + (log((2 ^ (1/g)) -1) / a)
  y <- c + ((f-c) / (abs(1 + exp(-a*tau + a*b.star))) ^ g)
  return(y)
}


build.line <- function(length.subj, data, start.sim, eps, id.start, sample.unif.noise) {
  line.c <- function(slope, point, b, sample.unif.noise) {
    sample.unif.noise <- sample.unif.noise * slope
    return(((slope + sample.unif.noise) * point) + b)
  }
  
  subjects.list <- list()
  for(i in start.sim:nrow(data)) {
    point <- data["Domain"][i,]
    slope <- data["slope"][i,]
    b  <- data["Response"][i,]
    b  <- b - (slope * point)
    id <- rep(id.start, length(point))
    sim.Domain <- seq(point - (length.subj / 2), point + (length.subj / 2), by=eps)
    sample.unif.noise <- runif(1, -.1, .1)
    sub <- map(sim.Domain, line.c,
               slope=slope, b=b, sample.unif.noise)
    sub <- data.frame("Domain"   = sim.Domain,
                      "Response" = as.numeric(sub),
                      "id"       = id)
    subjects.list[[i]] <- sub
    id.start <- id.start + 1
  }
  subjects.list <- do.call(bind_rows, subjects.list)
  return(subjects.list)
}



build.simulated.sigmoid <- function(a,b,c,f,g, dom) {
  resp <- c()
  for(i in 1:length(dom)) {
    resp <- append(resp, pl.5(a, b, c, f, g, dom[i]))
  }
  
  slope <- c((resp[2:length(resp)] - resp[1:length(resp)-1])/dom[2], 0)
  df <- data.frame("Domain" = dom,
                   "Response" = resp,
                   "slope" = slope)
  
  return(df)
}


align.zero <- function(data) {
  splitdata <- split(data, data$id)
  comb.list <- list()
  for(i in 1:length(splitdata)) {
    sub <- splitdata[[i]]
    sub$Domain <- sub$Domain - min(sub$Domain)
    comb.list[[i]] <- sub
  }
  comb.list <- do.call(bind_rows, comb.list)
  return(comb.list)
}





area.under.curve <- function(data) {
  data <- data[order(data$Domain, decreasing = FALSE),]
  total.area <- c()
  for(i in 1:(nrow(data) - 1)) {
    y.1 <- data["Response"][i,]
    y.2 <- data["Response"][(i + 1),]
    h   <- data["Domain"][(i + 1),] - data["Domain"][i,]
    t.i <- ((y.1 + y.2) / 2) * h
    total.area <- append(total.area, t.i)
  }
  return(sum(total.area))
}


area.between.curves <- function(data.theor, data.expir, perc.error = NULL) {
  data.theor      <- data.theor[order(data.theor$Domain, decreasing = FALSE),]
  data.expir      <- data.expir[order(data.expir$Domain, decreasing = FALSE),]
  data.cutoff     <- min(c(max(data.theor$Domain), max(data.expir$Domain)))
  #data.theor      <- subset(data.theor, Domain <= data.cutoff)
  #data.expir      <- subset(data.expir, Domain <= data.cutoff)
  area.theor      <- area.under.curve(data.theor)
  area.expir      <- area.under.curve(data.expir)
  if(!is.null(perc.error)) {
    return(((abs(area.theor - area.expir)) / area.theor) * 100)
  } else {
    return(abs(area.theor - area.expir))
  }
}

estimated.pred <- function(sim.ids) {
  est.pred <- data.frame(matrix(ncol = 2))
  pred.seq <- seq(0, max(sim.ids$Domain), by=.1)
  for(i in 1:(length(pred.seq) - 1)) {
    a <- pred.seq[i]
    b <- pred.seq[i + 1]
    sim.sub <- subset(sim.ids, Domain >= a & Domain <= b)
    mn.resp <- mean(sim.sub$Response)
    row.i <- c(((a + b) / 2), mn.resp)
    est.pred <- rbind(est.pred, row.i)
  }
  colnames(est.pred) <- c("Domain", "Response")
  est.pred <- est.pred[2:nrow(est.pred),]
  return(est.pred)
}

AddNoise <- function(data) {
}


keep.simulated.subjs <- function(data, a, b) {
  j <- 1
  keep.frame <- list()
  data$id <- factor(data$id)
  split.data <- split(data, data$id)
  for(i in 1:length(split.data)) {
    subj <- split.data[[i]]
    if(min(subj$Response) >= a & max(subj$Response) <= b) {
      keep.frame[[j]] <- subj
      j <- j + 1
    }
  }
  keep.frame <- do.call(bind_rows, keep.frame)
  return(keep.frame)
}


construct.simulated.dataset <- function(a,b,c,f,g, dom, length.subj, start.sim, eps, id.start){
  sigmoid.df <- build.simulated.sigmoid(a,b,c,f,g, dom)
  sim.ids   <-  build.line(length.subj, sigmoid.df, start.sim, eps, id.start)
  sim.ids$id            <- factor(sim.ids$id)
  data.aligned          <- align.zero(sim.ids)
  data.aligned          <- subset(data.aligned, Response >= 0)
  data.aligned$id  <- factor(data.aligned$id)
  init.curve.plot  <- ggplot(data = sigmoid.df, aes(x=Domain, y=Response)) + geom_point(colour="red")  + xlim(0,35) + ylab("Simulated Response Value") + xlab("Disease Progression (Years)")
  curve.with.lines <- ggplot(subset(sim.ids, Response > 0), aes(x=Domain, y=Response)) + theme(legend.position = "none") + geom_line(aes(group=factor(id)), colour="grey") + geom_point(data = sigmoid.df, colour="red") + xlim(0,35) + ylab("Simulated Response Value") + xlab("Disease Progression (Time)") 
  time.since.bline <- ggplot(data = data.aligned, aes(x=Domain, y=Response)) + geom_line(data = data.aligned, aes(group=factor(id)), colour = "grey") +theme(legend.position = "none") + ylab("Simulated Response Value") + xlab("Years (From Baseline)") 
  colnames(data.aligned) <- c("Time_Since_Baseline",
                              "Simulated_Response", 
                              "ID")
  data.aligned$ID <- factor(data.aligned$ID)
  return(list("data" = data.aligned,
              "init_curve" = init.curve.plot,
              "curve_lines" = curve.with.lines,
              "time_bline" = time.since.bline,
              "sigmoid_df" = sigmoid.df,
              "sim.ids"    = sim.ids))
}
