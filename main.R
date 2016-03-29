source("generate_data.R")
source("fit_stan.R")
source("forecast.R")

# Generate synthetic data
pop_size = 1e6
I.init = 1
R0 = 2.5
alpha = 1.25
a = 0.00
GI_span = 30 
GI_mean = 3
GI_var = 3
horizon = 90
last.obs <- 15
fcast.horizon <- 20


D <- generate.data.wrap(pop_size, 
						I.init,
						R0, 
						alpha, 
						a, 
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
			a        = a,
			pop_size = pop_size*0.4,
			# alpha    = alpha,
			# GI_mean  = GI_mean,
			GI_span   = GI_span
)

# Fit ERE model using Stan:
FIT <- ERE.fit.stan(model.filename = 'fit.stan', 
					dat = dat, 
					n.iter = 1000, 
					n.chains = 3,
					plot.trace = FALSE,
					plot.pairs = FALSE,
					plot.compTruth = TRUE
) 

# Forecast based on fit:
FCAST <- ERE.forecast(prm = FIT$prm.sample, # <-- sampled parameter values after the fit
					  fcast.horizon,
					  pop_size, 
					  R0, 
					  alpha, 
					  a, 
					  GI_span, 
					  GI_mean, 
					  GI_var,
					  last.obs,
					  syn.inc.full = D$syn.inc.full,
					  do.plot = TRUE,
					  seed=123)

