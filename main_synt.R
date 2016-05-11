source("RESuDe_FCT.R")
source("fit_stan.R")
source("forecast.R")

t1 <- as.numeric(Sys.time())
set.seed(1234)
pdf("plot_fcast_synt.pdf", width=12, height = 8)

# Generate synthetic data
# (parameters are read and asign from the file):
prmdata <- read.csv('prmdata.csv',header = FALSE)
pname <- trimws(as.character(prmdata[,1]))
for(i in 1:length(pname)) assign(pname[i], prmdata[i,2])

D <- RESuDe.simulate.wrap(pop_size, 
						  I.init,
						  R0, 
						  alpha, 
						  kappa, 
						  GI_span, 
						  GI_mean, 
						  GI_var,
						  horizon,
						  last.obs,
						  do.plot = TRUE,
						  seed=123)

dat.obs <- D$syn.inc
dat.obs <c(2,  9,  7 , 9 ,14, 13, 12, 15, 21, 16)
dat.full <- D$syn.inc.full

# forecasting horizon:
fcast.horizon <- horizon-last.obs

# Effective population bias:
bias.pop_size <- 1.0
# effective population bounds:
pop_mean <-  pop_size * bias.pop_size
pop_lsd  <-  2.0
pop_hi <- qlnorm(p=0.99, meanlog = log(pop_mean),sdlog = pop_lsd)
pop_lo <- qlnorm(p=0.01, meanlog = log(pop_mean),sdlog = pop_lsd)
pop_lo <- max(pop_lo,sum(dat.obs)+1)# <-- pop cannot be smaller than cumul incidence!
pop_hi <- round(pop_hi,0)
pop_lo <- round(pop_lo,0)
message(paste("\n\nEffective pop size:",pop_lo,"--",pop_mean,"--",pop_hi))

# Reporting rate:
rep.rate <- 0.4

# Define Stan's known data:
dat <- list(numobs   = last.obs,
			Iobs     = dat.obs,
			# rep.rate = rep.rate,
			R0_lo    = 0.7,
			R0_hi    = 10,
			GI_meanlo= 1,
			GI_meanhi= 15,
			GI_varhi = 4,
			GI_varlo = 2,
			alpha_hi = 4, #alpha*1.001,
			alpha_lo = 0, #alpha*0.999,
			kappa_hi = 0.1,
			kappa_lo = 0,
			pop_hi   = pop_hi,
			pop_lo   = pop_lo, 
			pop_mean = pop_mean,
			pop_lsd  = pop_lsd,
			alpha    = alpha,
			# GI_mean  = GI_mean,
			GI_span  = GI_span
)

# Fit RESuDe model using Stan:
mp <- read.csv("mcmc.csv", header = FALSE)
FIT <- RESuDe.fit.stan(model.filename = 'fit.stan', 
					   dat = dat, 
					   n.iter = mp[mp[,1]=='iter',2], 
					   n.chains = mp[mp[,1]=='nchains',2],
					   plot.compTruth = TRUE
) 
# Show diagnostic plots for Stan fit:
diagnostic <- mp[mp[,1]=='diagnostic',2]
# Diagnostic:
if(diagnostic){
	np <- names(FIT$prm.sample)
	np <- np[np!="Iout"]
	pairs(FIT$fit)	
	np <- np[np!="lp__"]
	traceplot(FIT$fit, pars=np, alpha=0.5,inc_warmup=TRUE)
}

# Forecast based on fit:
if(FALSE){
	FCAST <- RESuDe.forecast(prm = FIT$prm.sample, # <-- sampled parameter values after the fit
							 fcast.horizon,
							 pop_size = NULL, #pop_size.guess, 
							 kappa = NULL, 
							 alpha = NULL,
							 GI_var=NULL,
							 GI_span, 
							 last.obs,
							 syn.inc.full = dat.full,
							 do.plot = TRUE,
							 CI1 = 50, CI2=95,
							 seed=123)
}

if(TRUE){
fcast2 <- RESuDe.forecast.light(prm = FIT$prm.sample,
								fcast.horizon =fcast.horizon,
								kappa = 0, alpha = 0, 
								GI_span = GI_span,
								GI_mean = 3,
								GI_var = 3, 
								last.obs = last.obs,
								nmax = 200,
								syn.inc.full = dat.full, 
								do.plot = T,
								CI1 = 50,CI2=95,seed = 123)
}


dev.off()
t2 <- as.numeric(Sys.time())
message(paste("\n\n- - - - - Completed in",round((t2-t1)/60,1),"minutes - - - - - -\n\n"))
