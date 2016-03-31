source("RESuDe_FCT.R")


pop_size = 1e4
I.init = 1
R0 = 2.5
alpha = 0
kappa = 0.0
GI_span = 20 
GI_mean = 3
GI_var = 3
horizon = 120


sim0 <- RESuDe.simulate(pop_size, 
						I.init,
						R0, 
						alpha, 
						kappa, 
						GI_span, 
						GI_mean, 
						GI_var,
						horizon)

sim <- RESuDe.simulate(pop_size, 
					   I.init,
					   R0, 
					   alpha ,#= 2, 
					   kappa = 0.02, 
					   GI_span, 
					   GI_mean, 
					   GI_var,
					   horizon)


par(mfrow=c(1,1))
plot(sim0$I, typ="o", log="y", ylim=c(1,1000))
lines(sim$I, typ="l", col="red")
grid(lty=1)
