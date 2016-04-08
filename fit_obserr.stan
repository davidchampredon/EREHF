data {
	int numobs;  // number of observations
	int<lower=0> I_obs[numobs]; // Observed incidence 
	int<lower=0> GI_span;  // [0;GI_span] = Support of the generation interval distribution
	real rep_rate;
	real rep_var;
	// bounds for fitted parameters:
	real pop_lo;
	real pop_hi;
	real pop_mean;
	real pop_lsd;
	real R0_lo;
	real R0_hi;
	real GI_meanlo;
	real GI_meanhi;
	real GI_varlo;
	real GI_varhi;
	real alpha_lo;
	real alpha_hi;
	real kappa_lo;
	real kappa_hi;
}

parameters {
  real<lower=R0_lo,     upper=R0_hi>     R0;
  real<lower=GI_meanlo, upper=GI_meanhi> GI_mean;  
  real<lower=GI_varlo,  upper=GI_varhi>  GI_var;  
  real<lower=kappa_lo,  upper=kappa_hi>  kappa; // intervention rate
  real<lower=alpha_lo,  upper=alpha_hi>  alpha;  // heterogeneity
  real<lower=pop_lo,    upper=pop_hi>    pop_size;
  
  real<lower=0.001,     upper=0.999>     unit_beta[numobs];
}

transformed parameters {
  vector[numobs] I;
  vector[numobs] I_unobs;
  
  I[1] <- I_obs[1]; // right thing to do???
  for(t in 2:numobs){
  	I_unobs[t] <- unit_beta[t] * I[t];
  	I[t]       <- I_obs[t] + I_unobs[t];
  }
}

model {
	real z;
	vector[numobs] S;
	// vector[numobs] I; // unobserved ('latent') incidence
	vector[GI_span] GI_dist;
	real sGI;
	real I_tmp;
	vector[2] nextS;
	real GI_k;
	real GI_theta;
	
	real a;
	real b;
	
	// ==== Priors ====
	
	GI_mean  ~ uniform(GI_meanlo, GI_meanhi);
	GI_var   ~ uniform(GI_varlo, GI_varhi);
	R0       ~ uniform(R0_lo, R0_hi); 
	alpha    ~ uniform(alpha_lo,alpha_hi); 
	kappa    ~ uniform(kappa_lo, kappa_hi);
	pop_size ~ lognormal(log(pop_mean),pop_lsd);
	
	// Define generation interval distribution
	GI_k     <- GI_mean^2/GI_var;
	GI_theta <- GI_var/GI_mean;
	sGI      <- 0;
	for( t in 1:GI_span){
		sGI <- sGI + t^GI_k * exp(-GI_theta*t);
	}
	for( t in 1:GI_span){
		GI_dist[t] <- t^GI_k * exp(-GI_theta*t)/sGI;
	}
	
	// RENEWAL EQUATION:
	
	S[1] <- pop_size - I[1];
	for(t in 2:numobs){
		z <- 0;
		for(j in 1:min(GI_span,t-1)) {
			z <- z + GI_dist[j]*I[t-j];
		}
		I_tmp <- (S[t-1]/ pop_size)^(1+alpha) * R0 * exp(-kappa*t) * z ;
		
		print(t," ---> ",I_tmp);
		
		// True/Unobserved/Latent incidence.
		// Forced to be continuous because
		// Stan cannot handle discrete latent variables.
		// So, ideally we would like I[t]~poisson(I_tmp), but
		// approximated with continuous distribution
		// here gamma, with parameters chosen such that
		// mean and variance are the same as the poisson's
		
		I[t] ~ gamma(I_tmp,1.0);
		// I[t] ~ gamma( 5 ,1.0);
		
		// Observed incidence derived from the latent one.
		
		// converts mean and variance to Beta parameters (a and b):
		a <- rep_rate * ( (1/rep_rate-1)*rep_rate^2/rep_var - 1);
		b <- a*(1/rep_rate-1);
		
		unit_beta ~ beta(a,b);
		// I_obs[t] ~  normal( unit_beta * I[t], 0.001);   
		// increment_log_prob()
		
		
		// Cap the number of 
		// new infections to the 
		// remaining susceptibles:
		nextS[1] <- 0;
		nextS[2] <- S[t-1] - I[t];
		S[t] <- max(nextS);
	}
}

generated quantities {
  vector[numobs] Iout;
  for(i in 1:numobs) Iout[i] <- I_obs[i];
}