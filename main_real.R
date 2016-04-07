source("RESuDe_FCT.R")
source("fit_stan.R")
source("forecast.R")

t1 <- as.numeric(Sys.time())
set.seed(1234)
pdf("plot_fcast_real.pdf", width=12, height = 8)

prmdata <- read.csv('prmdata_real.csv',header = FALSE)
pname <- trimws(as.character(prmdata[,1]))
for(i in 1:length(pname)) assign(pname[i], prmdata[i,2])

# Read data from database
# (parameters are read and asign from the file):
source("../Datsid/read_db.R")
country <-  "LIBERIA"
disease <- "ebola"
synthetic <- 2
db0 <- get.epi.ts(db.path = "../Datsid/a.db",country =country,disease = disease,synthetic = synthetic)
db <- subset(db0, eventtype=='incidence' & eventtype2=='confirmed' & socialstruct=='gen_pop')
db$t <- date.to.duration(db$reportdate)

dat.obs <- db$count[1:last.obs]
dat.full <- db$count

plot(db$t, db$count, typ='s',main=paste(country,disease,db$source[1],synthetic,sep=";"))
plot(dat.full,typ='s')
points(dat.obs,pch=16)


# forecasting horizon:
fcast.horizon <- horizon-last.obs

# effective population bounds:
pop_mean <-  pop_size 
pop_lsd  <-  2.0
pop_hi <- qlnorm(p=0.99, meanlog = log(pop_mean),sdlog = pop_lsd)
pop_lo <- qlnorm(p=0.01, meanlog = log(pop_mean),sdlog = pop_lsd)
pop_lo <- max(pop_lo,sum(dat.obs)+1)# <-- pop cannot be smaller than cumul incidence!
pop_hi <- round(pop_hi,0)
pop_lo <- round(pop_lo,0)
message(paste("\n\nEffective pop size:",pop_lo,"--",pop_mean,"--",pop_hi))


# Define Stan's known data:
dat <- list(numobs   = last.obs,
			Iobs     = dat.obs,
			R0_lo    = 0.7,
			R0_hi    = 10,
			GI_meanhi= 4,
			GI_meanlo= 0.6,
			GI_varhi = 0.5,
			GI_varlo = 0.4,
			alpha_hi = 3, #alpha*1.001,
			alpha_lo = 0, #alpha*0.999,
			kappa_hi = 0.05,
			kappa_lo = 0,
			pop_hi   = pop_hi,
			pop_lo   = pop_lo, 
			pop_mean = pop_mean,
			pop_lsd  = pop_lsd,
			GI_span  = GI_span
)

# Fit RESuDe model using Stan:
mp <- read.csv("mcmc.csv", header = FALSE)
FIT <- RESuDe.fit.stan(model.filename = 'fit.stan', 
					   dat = dat, 
					   n.iter = mp[mp[,1]=='iter',2], 
					   n.chains = mp[mp[,1]=='nchains',2],
					   plot.compTruth = FALSE
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
FCAST <- RESuDe.forecast(prm = FIT$prm.sample, # <-- sampled parameter values after the fit
						 fcast.horizon,
						 GI_span, 
						 last.obs,
						 syn.inc.full = dat.full,
						 pop_size = NULL, #pop_size.guess, 
						 kappa = NULL, 
						 alpha = NULL,
						 GI_var = NULL,
						 do.plot = TRUE,
						 CI1 = 50, CI2=95,
						 seed=123)

dev.off()
t2 <- as.numeric(Sys.time())
message(paste("\n\n- - - - - Completed in",round((t2-t1)/60,1),"minutes - - - - - -\n\n"))