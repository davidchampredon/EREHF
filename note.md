# Notes

It seems there is no need to fit the population size. When try to fit it, MCMC fails. But when I enter on purpose a wrong value for pop_size, the forcast is still pretty good. 
 
Increasing parameter 'alpha'  (heterogeneity parameter) has the effect of 
 * lowering the peak incidence and bring slightly sooner the peak time, 
 * does _not_ affect the duration of the overall epidemic 
 * nor the exponential growth/decay. 

Fitting 'alpha' is pointless at the start of the epidemic, when the ratio S/N is still very close to 1. It is actually even harmfull for the fit, because alpha is basically unidentifiable until close to peak incidence.
 
Increasing parameter 'kappa' (intervention rate) has the effect of 
 * Slowing the initial growth from exponential to sub-exponential.
 * Reduce the peak incidence
 * push the peak incidence later
 * increase the duration of the epidemic