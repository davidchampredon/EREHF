data {
	int numobs;  // number of observations
	int<lower=0> Iobs[numobs]; // Observed incidence 
	int<lower=0> GI_span;  // [0;GI_span] = Support of the generation interval distribution
	real alpha;  // heterogeneity
    real<lower=0> GI_var; // GI ~ Gamma ; Variance
	//real<lower=0> pop_size;
}

parameters {
  real<lower=0,upper=10.0> R0;
  //real<lower=0, upper=6> alpha;  // heterogeneity
  real<lower=0.1, upper=15> GI_mean;  // GI ~ Gamma ; Mean
  real<lower=0.0001,upper=0.03> kappa; // intervention rate
  real<lower=100, upper=1000000> pop_size;
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
	
	//alpha ~ uniform(0,6); //lognormal(log(0.01),0.5);
	GI_mean ~ uniform(0,15);
	R0 ~ uniform(0.9,10); 
	kappa ~ uniform(0.0001,0.03);
	pop_size ~ lognormal(log(1e5),0.6);
	
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
		
		// ---- DEBUG ----
		if(is_nan(I_tmp)){
			print("DEBUG :::","t=",t,
			" S=",S[t-1], " alpha=",alpha, " popsize=",pop_size, 
			" R0=",R0," kappa=",kappa, " z=",z);
			print("GI_mean=",GI_mean, " GI_theta=",GI_theta, " sGI=",sGI);
			for(j in 1:min(GI_span,t-1)){
				print("t=", t, " j=",j, " GI_dist[j]=",GI_dist[j] , " Iobs[t-j]=",Iobs[t-j]);
			}
		}
		// ---- END DEBUG ----	
		
		Iobs[t] ~  poisson( I_tmp );   //increment_log_prob(poisson_log(Iobs[t],I_tmp));
		
		nextS[1] <- 0;
		nextS[2] <- S[t-1] - Iobs[t];
		S[t] <- max(nextS);
	}
}

generated quantities {
  vector[numobs] Iout;
  for(i in 1:numobs) Iout[i] <- Iobs[i];
}