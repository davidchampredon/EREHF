data {
	int numobs;  // number of observations
	int<lower=0> Iobs[numobs]; // Observed incidence 
	int<lower=0> GI_lag;
	real a; // intervention rate
  real<lower=0> GI_var;
	real<lower=0> pop_size;
}

parameters {
  real<lower=0> R0;
  real<lower=0> alpha;  // heterogeneity
  real<lower=0> GI_mean;
}

transformed parameters {
  
}

model {
	real z;
	vector[numobs] S;
	vector[GI_lag] GI_dist;
	real sGI;
	real I_tmp;
	vector[2] lambda;

	real GI_k;
	real GI_theta;
	
	// ==== Priors ====
	
	alpha ~ uniform(0,6); //lognormal(log(0.01),0.5);
	GI_mean ~ uniform(0,15); // lognormal(log(4),1); //uniform(0,5);
	R0 ~ uniform(0.9,10); //lognormal(log(4),0.5);
	// pop_size ~ uniform(990000,1010000); //lognormal(log(1e6),0.01);
	
	// Define generation interval distribution
	GI_k <- GI_mean^2/GI_var;
	GI_theta <- GI_var/GI_mean;
	
	sGI <- 0;
	for( t in 1:GI_lag){
		sGI <- sGI + t^GI_k * exp(-GI_theta*t);
	}
	for( t in 1:GI_lag){
		GI_dist[t] <- t^GI_k * exp(-GI_theta*t)/sGI;
	}
	
	S[1] <- pop_size - Iobs[1];
	
	for(t in 2:numobs){
		z <- 0;
		for(j in 1:min(GI_lag,t-1)){
			z <- z + GI_dist[j]*Iobs[t-j];
		}
		I_tmp <- (S[t-1]/ pop_size)^(1+alpha) * R0 * exp(-a*t) * z ;
		
		lambda[1] <- I_tmp;
		lambda[2] <- S[t-1];
		Iobs[t] ~  poisson( I_tmp ); //poisson( min(lambda) );
		
		S[t] <- S[t-1] - Iobs[t];
	}
}

generated quantities {
  vector[numobs] Iout;
  for(i in 1:numobs) Iout[i] <- Iobs[i];
}