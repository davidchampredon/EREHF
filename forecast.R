
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
							CI1 = 50,
							CI2 = 95,
							seed=123){
	message("\n\nStarting forecast simulations ...\n\n")
	library(snowfall)
	sfInit(parallel = TRUE, 
		   cpu = parallel::detectCores())
	
	# Future time:
	tfut <- (last.obs):(last.obs+fcast.horizon)
	
	# Matrix storing all forecasts:
	nsamples <- length(prm[[1]])
	sf <- matrix(NA, nrow = nsamples, ncol = last.obs+fcast.horizon)
	
	I.init <- prm$Iout[1,1:last.obs]
	
	# Wrap function for snowfall parralel execution:
	simul.wrap <- function(i, 
						   I.init,
						   prm, 
						   fcast.horizon,
						   pop_size, 
						   kappa,
						   alpha,
						   GI_span, 
						   GI_var,
						   last.obs,
						   seed){
		
			alpha_i <- alpha
			kappa_i <- kappa
			pop_size_i <- pop_size
			GI_var_i <- GI_var
			
			if(is.null(alpha)) alpha_i   <- prm$alpha[i]
			if(is.null(kappa)) kappa_i <- prm$kappa[i]
			if(is.null(pop_size)) pop_size_i <- prm$pop_size[i]
			if(is.null(GI_var)) GI_var_i <- prm$GI_var[i]
			
			R0_i      <- prm$R0[i]
			GI_mean_i <- prm$GI_mean[i]
			
			tmp <-  RESuDe.simulate(pop_size = pop_size_i,
									I.init   = I.init,
									R0       = R0_i,
									alpha    = alpha_i,
									kappa    = kappa_i,
									GI_span  = GI_span,
									GI_mean  = GI_mean_i,
									GI_var   = GI_var_i,
									horizon  = fcast.horizon,
									seed     = seed)
			return(tmp)
		}
	

	### Parallel execution:
	sfExportAll()
	idx.apply <- 1:nsamples
	res <- sfSapply(idx.apply, simul.wrap,
					I.init=I.init,
					prm=prm, # <-- sampled parameter values after the fit
					fcast.horizon=fcast.horizon,
					pop_size=pop_size,
					kappa=kappa,
					alpha=alpha,
					GI_span=GI_span,
					GI_var=GI_var,
					last.obs=last.obs,
					seed=seed,
					simplify = FALSE)
	sfStop()
	
	for(i in idx.apply) sf[i,] <-res[[i]]$I
	
	# Quantiles incidence cone:
	q <-  0.5+c(-CI2/2/100,-CI1/2/100,0,CI1/2/100,CI2/2/100)
	fcast.cone <-t(apply(sf, 2, function(v){quantile(v, probs = q)}))
	
	# Peak incidence (level and timing):
	pk.inc     <- apply(X = sf,MARGIN = 1,FUN = max)
	pktime.inc <- apply(X = sf,MARGIN = 1,FUN = which.max)
	
	pk.inc.true <- NULL
	pktime.inc.true <- NULL
	final.size.true <- NULL
	if(!is.null(syn.inc.full)) {
		pk.inc.true <- max(syn.inc.full)
		pktime.inc.true <- which.max(syn.inc.full)
		final.size.true <- sum(syn.inc.full)
	}
	
	# Final size:
	final.size <- rowSums(sf)
	
	if(do.plot){
		par(mfrow=c(2,2))
		
		ymx <- max(fcast.cone,syn.inc.full)
		xx <- c(1:length(I.init),tfut)
		
		if(is.null(syn.inc.full)) ymx <- max(fcast.cone)
		
		plot(syn.inc.full,
			 main = paste0("Forecast (median, ",CI1,"%CI, ",CI2,"%CI)"),
			 xlim = c(0,length(xx)),
			 ylim = c(1,ymx),
			 xlab = "time",
			 ylab = "", 
			 las = 1,
			 log = 'y',
			 typ = "o",
			 cex = 0.6, pch=16, 
			 col='lightgray')
		points(I.init, pch=16)
		abline(v=last.obs, lty=2)
		grid()
		
		# Median incidence:
		lines(fcast.cone[,3],lwd=2)
		lines(tfut,fcast.cone[tfut,3],lwd=2, col='blue')
		
		# Credible intervals incidence:
		col.cone <- rgb(0,0,0.8,0.3)
		polygon(x = c(tfut,rev(tfut)), 
				y = c(fcast.cone[tfut,1]+1,rev(fcast.cone[tfut,5])+1),
				border = NA,
				col = col.cone
		)
		polygon(x = c(tfut,rev(tfut)), 
				y=c(fcast.cone[tfut,2]+1,rev(fcast.cone[tfut,4])+1),
				border = NA,
				col = col.cone
		)
		
		plot.density.ci <- function(dat, q, title="", dat.true=NULL,logscale=FALSE){
			dfs <- density(dat)
			qt <- quantile(dat,probs = q)
			plot(dfs,
				 main = title,
				 ylab = '', yaxt='n',
				 xlab = '',
				 log = ifelse(logscale,'x',''),
				 typ='l',
				 lwd=3)
			
			ys <- 0.8*max(dfs$y)
			segments(x0=qt[1],x1=qt[5],y0=ys,y1=ys,lwd=2, col='gray')
			segments(x0=qt[2],x1=qt[4],y0=ys,y1=ys,lwd=7, col='gray')

			# Median
			abline(v=qt[3])
			
			if(!is.null(dat.true)){
				abline(v=dat.true,col="red",lty=2,lwd=3)
				text(x=dat.true,y=max(dfs$y)*1.01,
					 labels = paste0("true value=",dat.true),
					 pos = 4,
					 col='red')
			}
		}
		
		plot.hist.ci <- function(dat, q, title="", dat.true=NULL,logscale=FALSE){
			
			xlab <- ''
			if(logscale){
				dat      <- log10(dat)
				if(!is.null(dat.true)) dat.true <- log10(dat.true)
				xlab <- 'log10 scale'
			} 
			
			qt <- quantile(dat,probs = q)
			h <- hist(dat, 
					  main = title,
					  xlab = xlab,
					  ylab='',yaxt='n',
					  border = rgb(0,0,0,0.15), col=rgb(0,0,0,0.1))
		
			ys <- 0.8*max(h$counts)
			segments(x0=qt[1],x1=qt[5],y0=ys,y1=ys,lwd=2, col='blue')
			segments(x0=qt[2],x1=qt[4],y0=ys,y1=ys,lwd=7, col='blue')
			
			# Median
			abline(v=qt[3], col='blue')
			
			if(!is.null(dat.true)){
				dat.true.label <- dat.true
				if(logscale) dat.true.label <- paste(round(dat.true,1),";",round(10^dat.true))
				abline(v=dat.true,col="red",lty=2,lwd=3)
				text(x=dat.true,y=max(h$counts)*1.01,
					 labels = paste0("true value=",dat.true.label),
					 pos = 4,
					 col='red')
			}
		}
		
		plot.hist.ci(dat = final.size,
						q = q,
						title = "Final Size",
						dat.true = final.size.true,
						logscale = TRUE)
		
		plot.hist.ci(dat = pktime.inc,
						q = q,
						title = "Peak Timing",
						dat.true = pktime.inc.true)
		
		plot.hist.ci(dat = pk.inc,
						q = q,
						title = "Peak incidence",
						dat.true = pk.inc.true,
						logscale = TRUE)
		
		
	}
	message("... forecast simulations done!\n\n")
	return(list(sf = sf, 
				fcast.cone = fcast.cone,
				final.size = final.size,
				pk.inc = pk.inc,
				pktime.inc = pktime.inc))
}