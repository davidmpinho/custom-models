functions {
  vector additive_trend(int n, vector innovations, real scale) {
    /** 
     * Return a non-centered additive component with time-dependence (a 'trend').
     * If the innovations come from a standard normal -- normal(0, 1) -- this is 
     * equivalent to specifying: trend[t] ~ normal(trend[t-1], scale).   
     * 
     * @param n           Number of unique indexes
     * @param innovations Vector of size n with unit scale
     * @param scale       Scale of the distribution
     * @return            Scaled additive trend with mean equal to zero
    */
    vector[n] trend;
        
    trend = cumulative_sum(innovations);
    return (trend - mean(trend)) * scale;
  }  
  vector additive_factor(int n, vector innovations, real scale) {
    /** 
     * Return a non-centered additive component without time-dependence, 
     * assuming the first parameter is equal to zero. (The first parameter 
     * is the one that is most common.) This procedure is identical to 
     * specifying a non-centered hierarchical predictor with two differences:
     * I demean the innovations to avoid issues
     * of bimodality in those parameters; and I assume that the first parameter 
     * is equal to zero to also avoid issues of bimodality, this time arising 
     * from the interaction between alpha and that first parameter.
     * 
     * @param n           Number of unique indexes
     * @param innovations Vector of size n-1 with unit scale
     * @param scale       Scale of the distribution
     * @return            Scaled additive factor with mean equal to zero
    */
    vector[n] factor;
    
    factor[1] = 0;
    factor[2:n] = (innovations - mean(innovations)) * scale;
    return factor;
  }
} 
data {
  int<lower=0> N;    
  int<lower=0> N_holiday;    
  int<lower=0> y[N];                             // Outcome (int in this case) 
  int<lower=1,upper=366> day_year[N];            // Index for the day of year
  int<lower=1,upper=N_holiday+1> holiday[N];     // Index for holiday (1 means 'no holiday') 
}
parameters {
  real<lower=0> omega;
  real<lower=0> mu1;
  real rho;                                      // = log(alpha + beta)
  real kappa;                                    // = log(alpha / beta)
  real<lower=0> theta;                           
  real<lower=0> sigma_year;
  real<lower=0> sigma_holiday;
  vector[366] epsilon_year;
  vector[N_holiday] epsilon_holiday;
}
transformed parameters {
  real alpha;
  real beta;
  vector[N] mu;
  vector[N] phi;
  vector[N] mu_y;
  vector[366] phi_year;
  vector[N_holiday+1] phi_holiday;
  
  // Seasonal components -- equivalent to phi_year[t] ~ normal(phi_year[t-1], sigma_year)
  phi_year = additive_trend(366, epsilon_year, sigma_year);

  // Sporadic components (or anything that requires an index) 
  phi_holiday = additive_factor(N_holiday+1, epsilon_holiday, sigma_holiday);

  phi = exp(phi_year[day_year] + phi_holiday[holiday]);

  // Follows from rho being log(alpha+beta); and kappa being log(alpha/beta)
  alpha = exp(rho) / (1+exp(rho)) / (1+exp(kappa));
  beta = exp(rho) / (1+exp(rho)) * exp(kappa) / (1+exp(kappa));

  // Residual trend component 
  mu[1] = mu1;
  for (t in 2:N) {
    mu[t] = (omega * (1-(alpha+beta)) 
               + alpha * y[t-1] / phi[t-1]       // deseasoned y
               + beta * mu[t-1]);    
  }
  mu_y = mu .* phi;
}
model {
  omega ~ normal(2.25, 0.225);                   // Constant (unconditional expectation)
  mu1 ~ normal(1.5, 0.15);                       // First value of mu
  rho ~ normal(0, 1.6);                          // (Roughly flat on the [0, 1] scale)
  kappa ~ normal(0, 1.6);                        // (Roughly flat on the [0, 1] scale)
  theta ~ gamma(2, 0.1);                         // Overdispersion term
  sigma_year ~ normal(0, 0.02);                  // Scale parameters
  sigma_holiday ~ normal(0, 3);
  epsilon_year ~ std_normal();                   // Innovations
  epsilon_holiday ~ std_normal();
  y ~ neg_binomial_2(mu_y, theta);               // Likelihood
}

