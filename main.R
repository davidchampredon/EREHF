source("RESuDe_FCT.R")
source("fit_stan.R")
source("forecast.R")

# Generate synthetic data
pop_size = 1e4
I.init = 1
R0 = 2.5
alpha = 1.25
kappa = 0.0
GI_span = 20 
GI_mean = 3
GI_var = 3
horizon = 120

# 'Observed' data and forecasting horizon:
last.obs <- 21
fcast.horizon <- 50

# Simulate synthetic data:
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

# Define Stan's known data:
dat <- list(numobs   = last.obs,
			Iobs     = D$syn.inc,
			GI_var   = GI_var,
			#kappa    = kappa,
			pop_size = pop_size * 1.0 ,
			# alpha    = alpha,
			# GI_mean  = GI_mean,
			GI_span  = GI_span
)

# Fit RESuDe model using Stan:
mp <- read.csv("mcmc.csv", header = FALSE)
FIT <- RESuDe.fit.stan(model.filename = 'fit.stan', 
					   dat = dat, 
					   n.iter = mp[mp[,1]=='iter',2], 
					   n.chains = mp[mp[,1]=='nchains',2],
					   plot.trace = TRUE,
					   plot.pairs = FALSE,
					   plot.compTruth = TRUE
) 

# Forecast based on fit:
FCAST <- RESuDe.forecast(prm = FIT$prm.sample, # <-- sampled parameter values after the fit
					  fcast.horizon,
					  pop_size, 
					  kappa, 
					  GI_span, 
					  GI_var,
					  last.obs,
					  syn.inc.full = D$syn.inc.full,
					  do.plot = TRUE,
					  seed=123)

