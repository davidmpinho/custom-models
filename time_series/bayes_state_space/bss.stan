functions {
  vector additive_trend(int size, vector innovations, real sd) {
    vector[size] trend;
        
    trend = cumulative_sum(innovations);
    return (trend - mean(trend)) * sd;
  }  
} 
data {
  int<lower=0> N;    
  int<lower=0> y[N]; 
  int<lower=1> day_week[N];
  int<lower=1> day_month[N];
  int<lower=1> day_year[N];
}
parameters {
  real<lower=0> alpha;
  real<lower=0> sd_resid;
  real<lower=0> sd_week;
  real<lower=0> sd_month;
  real<lower=0> sd_year;
  vector[N] resid_innovations;
  vector[366] year_innovations;
  vector[30] month_innovations;
  vector[7] week_innovations;
}
transformed parameters {
  vector[N] mu;
  vector[N] phi;
  vector[366] phi_year;
  vector[30] phi_month;
  vector[7] phi_week;
  
  phi_year = additive_trend(366, year_innovations, sd_year);
  phi_month = additive_trend(30, month_innovations, sd_month);
  phi_week = additive_trend(7, week_innovations, sd_week);
  phi = phi_week[day_week] + phi_year[day_year] + phi_month[day_month];
  mu = additive_trend(N, resid_innovations, sd_resid);
}
model {
  alpha ~ normal(0, 0.1);
  sd_resid ~ normal(0, 0.1);
  sd_week ~ normal(0, 0.1);
  sd_month ~ normal(0, 0.1);
  sd_year ~ normal(0, 0.1);
  resid_innovations ~ std_normal();
  week_innovations ~ std_normal();
  month_innovations ~ std_normal();
  year_innovations ~ std_normal(); 
  y ~ poisson_log(alpha + mu + phi);
}

