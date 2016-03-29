

GI.dist <- function(t, GI.lag, GI.mean, GI.var){
	tvec <- 0:GI.lag
	GI.k <- GI.mean^2/GI.var
	GI.theta <- GI.var/GI.mean
	tmp <- tvec^GI.k * exp(-GI.theta*tvec)
	tmp2 <- t^GI.k * exp(-GI.theta*t)
	return(tmp2/sum(tmp))
}


generate.data <- function(pop.size, 
						  I.init,
						  R0, 
						  alpha, 
						  a, 
						  GI.lag, 
						  GI.mean, GI.var,
						  horizon, seed=123) {
	set.seed(seed)
	I <- vector()
	S <- vector()
	I[1] <- I.init
	S[1] <- pop.size - I.init
	
	for(t in 2:horizon){
		z <- 0
		for(j in 1:min(GI.lag,t-1)){
			z <- z + GI.dist(j, GI.lag, GI.mean, GI.var) * I[t-j]
		}
		I.tmp <- (S[t-1]/ pop.size)^(1+alpha) * R0 * exp(-a*t) * z 
		I[t] <- rpois(n=1, lambda =  min(I.tmp, S[t-1]) )
		S[t] <- S[t-1] - I[t]
		# print(paste(t,z,I[t]))
	}
	return(list(S=S, I=I))
}



sim.fcast <- function(pop.size, 
					  past.inc, # past.incidence
					  R0, 
					  alpha, 
					  a, 
					  GI.lag, 
					  GI.mean, GI.var,
					  fcast.horizon, 
					  seed=123){
	
	set.seed(seed)
	
	numobs <- length(past.inc)
	I <- past.inc
	S <- pop.size - cumsum(past.inc)

	for(t in (numobs+1):(numobs+fcast.horizon)){
		z <- 0
		for(j in 1:min(GI.lag,t-1)){
			z <- z + GI.dist(j, GI.lag, GI.mean, GI.var) * I[t-j]
		}
		I.tmp <- (S[t-1]/ pop.size)^(1+alpha) * R0 * exp(-a*t) * z 
		I[t] <- rpois(n=1, lambda =  min(I.tmp, S[t-1]) )
		S[t] <- S[t-1] - I[t]
	}
	return(list(S=S, I=I))
}