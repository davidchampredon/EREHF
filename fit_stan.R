library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source("generate_data.R")



compare.fit.truth <- function(z,truth,title=""){
	
	dd <- density(z,adjust = 0.4) 
	plot(dd,
		 main = deparse(substitute(truth)),
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

# Generate synthetic data
pop.size = 1e6
I.init = 1
R0 = 2.5
alpha = 1.25
a = 0.00
GI.lag = 30 
GI.mean = 3
GI.var = 3
horizon = 90

par(mfrow=c(1,1))
tt <- 0:GI.lag
plot(tt,GI.dist(t = tt,GI.lag = GI.lag,
				GI.mean = GI.mean, 
				GI.var = GI.var),
	 typ='l',lwd=6,col='blue',xlab="",ylab="",main="GI distrib")
grid()

syn.data <- generate.data(pop.size = pop.size, 
						  I.init = I.init,
						  R0 = R0,
						  alpha = alpha,
						  a = a, 
						  GI.lag = GI.lag, 
						  GI.mean = GI.mean, 
						  GI.var = GI.var,
						  horizon = horizon,
						  seed = 12345)

syn.inc.full <- syn.data$I
last.obs <- 15
syn.inc <- syn.inc.full[1:last.obs]
numobs <- length(syn.inc)
par(mfrow=c(1,2))
plot(syn.inc.full,typ='s', lwd=1, main="Incidence data",
	 xlab="time",ylab="")
lines(syn.inc,typ='s', lwd=6)
plot(syn.inc.full,typ='o',log="y", lwd=1, xlab="time",ylab="",cex=0.5)
lines(syn.inc,typ='o', lwd=6) ; grid()



# Define Stan's known data:
inc.dat <- list(numobs   = numobs,
				Iobs     = syn.inc,
				# GI_mean  = GI.mean,
				GI_var   = GI.var,
				a        = a,
				# alpha    = alpha,
				pop_size = pop.size*0.4,
				GI_lag   = GI.lag
)


# Fit model with Stan:
fit <- stan(file     = 'fit.stan', 
			data     = inc.dat, 
			iter     = 7000,
			cores    = 4,
			chains   = 3,
			control  = list(adapt_delta=0.8)
)

#plot(fit)
traceplot(fit, pars=c("R0","alpha","GI_mean"))
print(fit)
#pairs(fit)

prm <- extract(fit)
names(prm)

par(mfrow=c(2,2))
compare.fit.truth(z = prm$alpha, truth = alpha)
compare.fit.truth(z = prm$R0, truth = R0)
compare.fit.truth(z = prm$GI_mean, truth = GI.mean)
# compare.fit.truth(z = prm$pop_size, truth = pop.size)
# compare.fit.truth(z = prm$GI_var, truth = GI.var)


# === FORECAST ====

fcast.horizon <- 20
tfut <- (last.obs):(last.obs+fcast.horizon)

nsamples <- length(prm[[1]])
sf <- matrix(NA, nrow = nsamples, ncol = last.obs+fcast.horizon)

for(i in 1:nsamples){
	alpha_i   <- prm$alpha[i]
	R0_i      <- prm$R0[i]
	GI_mean_i <- prm$GI_mean[i]
	
	tmp <-  sim.fcast(pop.size = pop.size,
						past.inc = prm$Iout[1,1:numobs],
						R0       = R0_i,
						alpha    = alpha_i,
						a        = a,
						GI.lag   = GI.lag,
						GI.mean  = GI_mean_i,
						GI.var   = GI.var,
						fcast.horizon =  fcast.horizon,
						seed = 12345)
	 sf[i,] <-tmp$I
	 if(i%%1000==0) message(paste("forecasting:",i,"/",nsamples))
}

q <- c(0.025, 0.25, 0.5, 0.75, 0.975)

fcast.cone <-t(apply(sf, 2, function(v){quantile(v, probs = q)}))
fcast.cone
ymx <- max(fcast.cone,syn.inc.full)

par(mfrow=c(1,1))
plot(syn.inc.full,
	 ylim = c(1,ymx),
	 xlab="time",
	 log='y',
	 typ="o",
	 cex=1, pch=16,col='lightgray')
points(syn.inc,pch=16)
abline(v=last.obs, lty=2)
grid()

lines(fcast.cone[,3],lwd=1)
lines(tfut,fcast.cone[tfut,3],lwd=2, col='blue')


col.cone <- rgb(0,0,1,0.3)
polygon(x = c(tfut,rev(tfut)), 
		y=c(fcast.cone[tfut,1]+1,rev(fcast.cone[tfut,5])),
		border = NA,
		col = col.cone
		)
polygon(x = c(tfut,rev(tfut)), 
		y=c(fcast.cone[tfut,2],rev(fcast.cone[tfut,4])),
		border = NA,
		col = col.cone
)

# lines(fcast.cone[,1])
# lines(fcast.cone[,2])
# lines(fcast.cone[,4])
# lines(fcast.cone[,5])
