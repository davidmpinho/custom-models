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
  int<lower=1,upper=31> day_month[N];            // Index for the day of month
  int<lower=1,upper=7> day_week[N];              // Index for the day of week                    
  int<lower=1,upper=N_holiday+1> holiday[N];     // Index for holiday (1 means 'no holiday') 
}
parameters {
  real<lower=0> theta;                           // Parameter 
  real alpha;                                    // Constant
  real<lower=0> sigma_resid;                     // Scale for non-seasonal component
  real<lower=0> sigma_year;                      // Scale for seasonal component 
  real<lower=0> sigma_month;            
  real<lower=0> sigma_week;
  real<lower=0> sigma_holiday;                   // Scale for holidays 
  vector[N] epsilon_resid;                       // Innovations
  vector[366] epsilon_year;
  vector[31] epsilon_month;
  vector[7] epsilon_week;
  vector[N_holiday] epsilon_holiday;            
}
transformed parameters {
  vector[N] mu;
  vector[N] phi;
  vector[N] mu_y;
  vector[366] phi_year;
  vector[31] phi_month;
  vector[7] phi_week;
  vector[N_holiday+1] phi_holiday;
  
  // Seasonal components -- equivalent to mu[t] ~ normal(mu[t-1], sigma_resid)
  mu = additive_trend(N, epsilon_resid, sigma_resid);    
  phi_year = additive_trend(366, epsilon_year, sigma_year);
  phi_month = additive_trend(31, epsilon_month, sigma_month);
  phi_week = additive_trend(7, epsilon_week, sigma_week);
  
  // Sporadic components (or anything that requires an index) 
  phi_holiday = additive_factor(N_holiday+1, epsilon_holiday, sigma_holiday);
  
  phi = (phi_year[day_year] + phi_month[day_month] 
           + phi_week[day_week] + phi_holiday[holiday]);
  mu_y = exp(alpha + mu + phi);
}
model {
  alpha ~ normal(log(2.25), 0.1);                // 2.25 is the unconditionnal mean of y
  sigma_resid ~ normal(0, 0.02);                 // Scale parameters (mildly informative)
  sigma_year ~ normal(0, 0.02);
  sigma_month ~ normal(0, 0.07);
  sigma_week ~ normal(0, 0.15);
  sigma_holiday ~ normal(0, 5);
  epsilon_resid ~ std_normal();                  // Innovations 
  epsilon_year ~ std_normal(); 
  epsilon_month ~ std_normal();
  epsilon_week ~ std_normal();
  epsilon_holiday ~ std_normal();
  theta ~ gamma(2, 0.1);                         // Overdispersion term 
  y ~ neg_binomial_2(mu_y, theta);               // Likelihood 
}

