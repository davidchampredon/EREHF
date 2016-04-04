

RESuDe.forecast <- function(prm, # <-- sampled parameter values after the fit
							fcast.horizon,
							pop_size, 
							kappa,
							alpha,
							GI_span, 
							GI_var,
							last.obs,
							syn.inc.full = NULL,
							do.plot = FALSE,
							seed=123){
	# Future time:
	tfut <- (last.obs):(last.obs+fcast.horizon)
	
	# Matrix storing all forecasts:
	nsamples <- length(prm[[1]])
	sf <- matrix(NA, nrow = nsamples, ncol = last.obs+fcast.horizon)
	
	I.init <- prm$Iout[1,1:last.obs]
	
	alpha_i <- alpha
	kappa_i <- kappa
	pop_size_i <- pop_size
	
	for(i in 1:nsamples){
		
		if(is.null(alpha)) alpha_i   <- prm$alpha[i]
		if(is.null(kappa)) kappa_i <- prm$kappa[i]
		if(is.null(pop_size)) pop_size_i <- prm$pop_size[i]
		
		R0_i      <- prm$R0[i]
		GI_mean_i <- prm$GI_mean[i]
		
		tmp <-  RESuDe.simulate(pop_size = pop_size_i,
								I.init   = I.init,
								R0       = R0_i,
								alpha    = alpha_i,
								kappa    = kappa_i,
								GI_span  = GI_span,
								GI_mean  = GI_mean_i,
								GI_var   = GI_var,
								horizon  = fcast.horizon,
								seed     = seed)
		sf[i,] <-tmp$I
		if(i%%1000==0) message(paste("forecasting:",i,"/",nsamples))
	}
	
	q <- c(0.025, 0.25, 0.5, 0.75, 0.975)
	
	fcast.cone <-t(apply(sf, 2, function(v){quantile(v, probs = q)}))
	
	if(do.plot){
		par(mfrow=c(1,1))
		ymx <- max(fcast.cone,syn.inc.full)
		plot(syn.inc.full,
			 main = "Forecast (median, 50%CI, 95%CI)",
			 ylim = c(1,ymx),
			 xlab="time",
			 log='y',
			 typ="o",
			 cex=1, pch=16,col='lightgray')
		points(I.init, pch=16)
		abline(v=last.obs, lty=2)
		grid()
		
		lines(fcast.cone[,3],lwd=1)
		lines(tfut,fcast.cone[tfut,3],lwd=2, col='blue')
		
		col.cone <- rgb(0,0,0.8,0.3)
		polygon(x = c(tfut,rev(tfut)), 
				y=c(fcast.cone[tfut,1]+1,rev(fcast.cone[tfut,5])+1),
				border = NA,
				col = col.cone
		)
		polygon(x = c(tfut,rev(tfut)), 
				y=c(fcast.cone[tfut,2]+1,rev(fcast.cone[tfut,4])+1),
				border = NA,
				col = col.cone
		)
	}
	
	return(list(sf = sf, fcast.cone = fcast.cone))
}