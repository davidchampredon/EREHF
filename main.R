source("RESuDe_FCT.R")
source("fit_stan.R")
source("forecast.R")

t1 <- as.numeric(Sys.time())
set.seed(1234)

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

# forecasting horizon:
fcast.horizon <- 30

bias.pop_size <- 0.2
pop_size.guess <-  pop_size * bias.pop_size

# Show diagnostic plots for Stan fit:
diagnostic <- FALSE

# Define Stan's known data:
dat <- list(numobs   = last.obs,
			Iobs     = D$syn.inc,
			GI_var   = GI_var,
			#kappa    = kappa,
			#pop_size = pop_size.guess,
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

# Diagnostic:
if(diagnostic){
	np <- names(FIT$prm.sample)
	np <- np[np!="Iout"]
	pairs(FIT$fit)	
	np <- np[np!="lp__"]
	traceplot(FIT$fit, pars=np, alpha=0.5,inc_warmup=TRUE)
}

# Forecast based on fit:
FCAST <- RESuDe.forecast(prm = FIT$prm.sample, # <-- sampled parameter values after the fit
						 fcast.horizon,
						 pop_size = NULL, #pop_size.guess, 
						 kappa = NULL, 
						 alpha = alpha,
						 GI_span, 
						 GI_var,
						 last.obs,
						 syn.inc.full = D$syn.inc.full,
						 do.plot = TRUE,
						 seed=123)

t2 <- as.numeric(Sys.time())
message(paste("\n\n- - - - - Completed in",round((t2-t1)/60,1),"minutes - - - - - -"))
