data {
	int numobs;  // number of observations
	int<lower=0> Iobs[numobs]; // Observed incidence 
	int<lower=0> GI_span;  // [0;GI_span] = Support of the generation interval distribution
	
	// bounds for fitted parameters:
	real pop_lo;
	real pop_hi;
	real pop_mean;
	real pop_lsd;
	real R0_lo;
	real R0_hi;
	
	real GI_mean;
	real GI_var;
	// real GI_meanlo;
	// real GI_meanhi;
	// real GI_varlo;
	// real GI_varhi;
	
	real alpha;
	real kappa;
	// real alpha_lo;
	// real alpha_hi;
	// real kappa_lo;
	// real kappa_hi;
}

parameters {
  real<lower=R0_lo,  upper=R0_hi>     R0;
  real<lower=pop_lo, upper=pop_hi>    pop_size;
  // real<lower=GI_meanlo, upper=GI_meanhi> GI_mean;  
  // real<lower=GI_varlo,  upper=GI_varhi>  GI_var;  
  // real<lower=kappa_lo,  upper=kappa_hi>  kappa; // intervention rate
  // real<lower=alpha_lo,  upper=alpha_hi>  alpha;  // heterogeneity
}

transformed parameters {
  
}

model {
	real z;
	vector[numobs] S;
	vector[GI_span] GI_dist;
	real sGI;
	real I_tmp;
	vector[2] nextS;
	real GI_k;
	real GI_theta;
	
	// ==== Priors ====
	
	R0       ~ uniform(R0_lo, R0_hi); 
	pop_size ~ lognormal(log(pop_mean),pop_lsd);
	// GI_mean  ~ uniform(GI_meanlo, GI_meanhi);
	// GI_var   ~ uniform(GI_varlo, GI_varhi);
	// alpha    ~ uniform(alpha_lo,alpha_hi); 
	// kappa    ~ uniform(kappa_lo, kappa_hi);
	
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
	
	S[1] <- pop_size - Iobs[1];
	for(t in 2:numobs){
		z <- 0;
		for(j in 1:min(GI_span,t-1)) {
			z <- z + GI_dist[j]*Iobs[t-j];
		}
		I_tmp <- (S[t-1]/ pop_size)^(1+alpha) * R0 * exp(-kappa*t) * z ;
		Iobs[t] ~  poisson( I_tmp );   //increment_log_prob(poisson_log(Iobs[t],I_tmp));
		// Cap the number of 
		// new infections to the 
		/// remaining susceptibles:
		nextS[1] <- 0;
		nextS[2] <- S[t-1] - Iobs[t];
		S[t] <- max(nextS);
	}
}

generated quantities {
  vector[numobs] Iout;
  for(i in 1:numobs) Iout[i] <- Iobs[i];
}