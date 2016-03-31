library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


compare.fit.truth <- function(prm.sample, prm.name){
	
	truth <- get(prm.name)
	z <- prm.sample[[prm.name]]
	dd <- density(z,adjust = 0.4) 
	plot(dd,
		 main = prm.name,
		 lwd=3,
		 xlab="",
		 yaxt="n", ylab="",
		 xlim=range(truth,dd$x))
	m <- mean(z)
	md <- median(z)
	abline(v=c(m,md), col=c("blue","black"),lty=c(3,1))
	abline(v=truth,lty=2,col="red",lwd=6)
	q <- quantile(x = z, probs = c(0.025,0.1,0.9,0.975))
	mycol <-  rgb(0,0,0,0.1)
	polygon(x = c(q[1],q[4],q[4],q[1]), 
			y=c(0,0,max(dd$y),max(dd$y)), 
			border = NA,
			col =mycol)
	polygon(x = c(q[2],q[3],q[3],q[2]), 
			y=c(0,0,max(dd$y),max(dd$y)), 
			border = NA,
			col = mycol)
}

RESuDe.fit.stan <- function(model.filename, 
						  dat, 
						  n.iter, 
						  n.chains,
						  plot.trace = FALSE,
						  plot.pairs = FALSE,
						  plot.compTruth = FALSE
						  ) {
	
	fit <- stan(file     = model.filename, 
				data     = dat, 
				iter     = n.iter,
				cores    = parallel::detectCores() ,
				chains   = n.chains,
				control  = list(adapt_delta=0.8)
	)
	
	print(fit)
	
	prm <- extract(fit)
	np <- names(prm)
	np <- np[np!="Iout"]
	
	if(plot.trace) traceplot(fit, pars=np)
	if(plot.pairs) pairs(fit)
	
	np <- np[np!="lp__"]
	par(mfrow=c(2,2))
	if(plot.compTruth) for(x in np) compare.fit.truth(prm,x)
	
	return(list(fit = fit, prm.sample = prm))
}


