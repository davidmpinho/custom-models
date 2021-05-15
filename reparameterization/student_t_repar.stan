data {
  int<lower=0> N;
  vector[N] x;
  real y[N];
}
parameters {
  // theta = sqrt(nu / (nu-2) * sigma^2), sigma being the scale parameter of the student_t
  real<lower=0> theta;                    // Standard deviation of the student_t
  real log_nu;                            // Student-t's degrees of freedom: exp(log_nu) + 2
  real beta;
}
transformed parameters {
  vector[N] mu_y;
  real<lower=0> sigma; 
  
  sigma = theta / sqrt((exp(log_nu)+2) / exp(log_nu));               
  mu_y = beta * x;
}
model {
  theta ~ normal(0, 1);              
  log_nu ~ normal(2.57, 0.8);           // Roughly implies: exp(log_nu) ~ gamma(2, 0.1)
  beta ~ normal(0, 0.5);             
  y ~ student_t(exp(log_nu), mu_y, sigma);    
}
