###
### Finding equivalent _continuous_ distributions
### to a Poisson (discrete) distribution
###

set.seed(1234)
n <- 1e5

# Poisson is defined by one parameter only: lambda
# E(Poisson) = lambda
# var(Poisson) = lambda
lambda <- 30
s.pois <- rpois(n, lambda)

# Gamma is defined by scale 'a'  and rate 'b'
# E(Gamma) = a/b
# Var(Gamma) = a/b^2
# Choosing a and b to match the 2 firsts moments of Poisson
# we must have:
# a = lambda and b=1
a <- lambda
b <- 1
s.gam <- rgamma(n,shape=a,rate=b)

# Lognormal is defined by 2 parameters mu and sigma:
# E(Logn) = exp(mu+sigma^2/2)
# Var(Logn) = (exp(sigma^2)-1)exp(2mu+sigma^2)
# Choosing a and b to match the 2 firsts moments of Poisson
# we must have:
# mu = log(lambda)-0.5*log(1+1/lambda)
# sigma = sqrt( log(1+1/lambda) )
mu <- log(lambda)-0.5*log(1+1/lambda)
sigma <- sqrt( log(1+1/lambda) )
s.ln <- rlnorm(n,meanlog = mu, sdlog = sigma)


summary(s.pois)
summary(s.gam)
summary(s.ln)

sd(s.pois)
sd(s.gam)
sd(s.ln)

d.pois <- density(s.pois)
d.gam  <- density(s.gam)
d.ln   <- density(s.ln)

plot(d.pois, 
	 xlim=range(d.pois$x,d.gam$x,d.ln$x), 
	 ylim=range(d.pois$y,d.gam$y), 
	 lwd=2,
	 typ='l', log='y')
lines(d.gam,col='red',lwd=2)
lines(d.ln,col='green',lwd=2)
