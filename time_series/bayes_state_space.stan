functions {
  vector additive_trend(int size, vector innovations, real sd) {
    // Takes in realizations ("innovations") from a standard normal distribution,
    //   makes a "trend" by applying the cumulative sum to those innovations, 
    //   and then centers that trend and scaled it by the standard deviation ("sd")
    //   specified.
    vector[size] trend;
        
    trend = cumulative_sum(innovations);
    return (trend - mean(trend)) * sd;
  }  
} 
data {
  int<lower=0> N;    
  int<lower=0> y[N];                             // Outcome (int in this case) 
  int<lower=1,upper=366> day_year[N];            // Index for the day of year
  int<lower=1,upper=31> day_month[N];            // Index for the day of month
  int<lower=1,upper=7> day_week[N];              // Index for the day of week                    
}
parameters {
  real<lower=0> alpha;
  real<lower=0> sigma_resid;
  real<lower=0> sigma_year;
  real<lower=0> sigma_month;
  real<lower=0> sigma_week;
  vector[N] epsilon_resid;
  vector[366] epsilon_year;
  vector[31] epsilon_month;
  vector[7] epsilon_week;
}
transformed parameters {
  vector[N] mu;
  vector[N] phi;
  vector[366] phi_year;
  vector[31] phi_month;
  vector[7] phi_week;
  
  mu = additive_trend(N, epsilon_resid, sigma_resid);
  phi_year = additive_trend(366, epsilon_year, sigma_year);
  phi_month = additive_trend(30, epsilon_month, sigma_month);
  phi_week = additive_trend(7, epsilon_week, sigma_week);
  phi = phi_year[day_year] + phi_month[day_month] + phi_week[day_week];
}
model {
  alpha ~ normal(0, 0.1);                        // Constant
  sigma_resid ~ normal(0, 0.1);                  // Scale parameters
  sigma_year ~ normal(0, 0.1);
  sigma_month ~ normal(0, 0.1);
  sigma_week ~ normal(0, 0.1);
  epsilon_resid ~ std_normal();                  // Innovations -- equivalent to `mu ~ normal(0, sigma_resid)`
  epsilon_year ~ std_normal(); 
  epsilon_month ~ std_normal();
  epsilon_week ~ std_normal();
  y ~ poisson_log(alpha + mu + phi);             // Likelihood
}

