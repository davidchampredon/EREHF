library(rstan)
# rstan_options(auto_write = TRUE)
# options(mc.cores = parallel::detectCores())


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

plot.posterior <- function(prm.sample, prm.name){
	
	z <- prm.sample[[prm.name]]
	dd <- density(z,adjust = 0.4) 
	plot(dd,
		 main = prm.name,
		 lwd=3,
		 xlab="",
		 yaxt="n", ylab="",
		 xlim=range(dd$x))
	m <- mean(z)
	md <- median(z)
	abline(v=c(m,md), col=c("blue","black"),lty=c(3,1))
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
							n.cores = NULL,
							plot.compTruth = FALSE
) {
	
	if(is.null(n.cores)) nc <- parallel::detectCores()
	if(!is.null(n.cores)) nc <- n.cores
	
	fit <- stan(file     = model.filename, 
				data     = dat, 
				iter     = n.iter,
				cores    = nc,
				chains   = n.chains,
				control  = list(adapt_delta=0.8)
	)
	print(fit)
	
	prm <- extract(fit)
	np <- names(prm)
	np <- np[np!="Iout"]
	np <- np[np!="lp__"]
	nn <- sqrt(length(np))
	par(mfrow=c(round(nn,0),ceiling(nn)))
	if(plot.compTruth) {
		for(x in np) compare.fit.truth(prm,x)
	}
	if(!plot.compTruth) {
		for(x in np) plot.posterior(prm,x)
	}
	return(list(fit = fit, prm.sample = prm))
}


