data {
  int<lower=0> N;    
  int<lower=0> y[N];                             // Outcome (int in this case) 
}
parameters {
  real<lower=0> omega;
  real<lower=0> mu1;
  // As rho increases, the cond. mean (omega) has a smaller influence on the uncond. mean
  // As kappa increases, the unconditional mean becomes more influenced by mu_y[t-1]  
  real rho;                                      // = log((alpha+beta) / (1-(alpha+beta)))
  real kappa;                                    // = log(alpha / beta)
  real<lower=0> theta;                           
}
transformed parameters {
  real beta;
  real alpha;
  vector[N] mu_y;
  
  beta = exp(rho) / ((1+exp(rho)) * (1+exp(kappa)));  
  alpha = beta * exp(kappa);

  // Residual trend component 
  mu_y[1] = mu1;
  for (t in 2:N) {
    mu_y[t] = (omega * (1-(alpha+beta)) 
               + alpha * y[t-1] 
               + beta * mu_y[t-1]);    
  }
}
model {
  omega ~ normal(2.25, 0.225);                   // Constant (unconditional expectation)
  mu1 ~ normal(1.5, 0.15);                       
  rho ~ normal(0, 1.6);                          // (Roughly flat on the [0, 1] scale)
  kappa ~ normal(0, 1.6);                        // (Roughly flat on the [0, 1] scale)
  theta ~ gamma(2, 0.1);                         // Overdispersion term
  y ~ neg_binomial_2(mu_y, theta);               // Likelihood
}

