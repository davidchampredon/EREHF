data {
	int numobs;  // number of observations
	int<lower=0> Iobs[numobs]; // Observed incidence 
	int<lower=0> GI_span;  // [0;GI_span] = Support of the generation interval distribution
	
    real<lower=0> GI_var; // GI ~ Gamma ; Variance
	real<lower=0> pop_size;
	
}

parameters {
  real<lower=0> R0;
  real<lower=0> alpha;  // heterogeneity
  real<lower=0> GI_mean;  // GI ~ Gamma ; Mean
  real kappa; // intervention rate
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
	alpha ~ uniform(0,6); //lognormal(log(0.01),0.5);
	GI_mean ~ uniform(0,15);
	R0 ~ uniform(0.9,10); 
	kappa ~ uniform(0.0001,0.03);
	
	// Define generation interval distribution
	GI_k <- GI_mean^2/GI_var;
	GI_theta <- GI_var/GI_mean;
	sGI <- 0;
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
		for(j in 1:min(GI_span,t-1)){
			z <- z + GI_dist[j]*Iobs[t-j];
		}
		I_tmp <- (S[t-1]/ pop_size)^(1+alpha) * R0 * exp(-kappa*t) * z ;
		//Iobs[t] ~  poisson( I_tmp ); 
		
		increment_log_prob(poisson_log(Iobs[t],I_tmp));
		
		nextS[1] <- 0;
		nextS[2] <- S[t-1] - Iobs[t];
		S[t] <- max(nextS);
	}
}

generated quantities {
  vector[numobs] Iout;
  for(i in 1:numobs) Iout[i] <- Iobs[i];
}