

GI_dist <- function(t, GI_span, GI_mean, GI_var){
	tvec <- 0:GI_span
	GI_k <- GI_mean^2/GI_var
	GI_theta <- GI_var/GI_mean
	tmp <- tvec^GI_k * exp(-GI_theta*tvec)
	tmp2 <- t^GI_k * exp(-GI_theta*t)
	return(tmp2/sum(tmp))
}


generate.data <- function(pop_size, 
						  I.init,
						  R0, 
						  alpha, 
						  a, 
						  GI_span, 
						  GI_mean, 
						  GI_var,
						  horizon,
						  seed=123) {
	set.seed(seed)
	
	if(length(I.init)==1){
		# used to generate synthetic data
		I <- vector()
		S <- vector()
		I[1] <- I.init
		S[1] <- pop_size - I.init
		numobs <- 1
	}
	if(length(I.init)>1){
		# Used when forecasting
		numobs <- length(I.init)
		I <- I.init
		S <- pop_size - cumsum(I.init)
	}
	
	for(t in (numobs+1):(numobs+horizon)){
		z <- 0
		for(j in 1:min(GI_span,t-1)){
			z <- z + GI_dist(j, GI_span, GI_mean, GI_var) * I[t-j]
		}
		I.tmp <- (S[t-1]/ pop_size)^(1+alpha) * R0 * exp(-a*t) * z 
		I[t] <- rpois(n=1, lambda =  min(I.tmp, S[t-1]) )
		S[t] <- max(0,S[t-1] - I[t])
	}
	return(list(S=S, I=I))
}


generate.data.wrap <- function(pop_size, 
							   I.init,
							   R0, 
							   alpha, 
							   a, 
							   GI_span, 
							   GI_mean, 
							   GI_var,
							   horizon,
							   last.obs,
							   do.plot = FALSE,
							   seed=123){
	
	
	# Generate full epidemic
	syn.data <- generate.data(pop_size = pop_size, 
							  I.init = I.init,
							  R0 = R0,
							  alpha = alpha,
							  a = a, 
							  GI_span = GI_span, 
							  GI_mean = GI_mean, 
							  GI_var = GI_var,
							  horizon = horizon,
							  seed = seed)
	syn.inc.full <- syn.data$I
	
	# Just take the start of epidemic
	# (will forecast the rest)
	syn.inc <- syn.inc.full[1:last.obs]
	numobs <- length(syn.inc)
	
	if(do.plot){
		# Generation interval distribution:
		par(mfrow=c(1,1))
		tt <- 0:GI_span
		plot(tt,GI_dist(t = tt,GI_span = GI_span,
						GI_mean = GI_mean, 
						GI_var = GI_var),
			 typ='l',lwd=6,col='blue',xlab="",ylab="",main="GI distrib")
		grid()
		
		# Incidence:
		par(mfrow=c(1,2))
		plot(syn.inc.full,typ='s', lwd=1, main="Incidence data",
			 xlab="time",ylab="")
		lines(syn.inc,typ='s', lwd=6)
		plot(syn.inc.full,typ='o',log="y", lwd=1, xlab="time",ylab="",cex=0.5)
		lines(syn.inc,typ='o', lwd=6) ; grid()
	}
	
	return(list(syn.inc.full = syn.inc.full,
				syn.inc = syn.inc))
}