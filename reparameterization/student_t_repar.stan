data {
  int<lower=0> N;
  vector[N] x;
  real y[N];
}
parameters {
  // theta = sqrt(nu / (nu-2) * sigma^2)
  // kappa = sqrt(nu / (nu-2)) - 1
  real<lower=0> theta;                    // Standard deviation of the student_t
  real<lower=0> kappa;                    // Standard deviation minus 1 of student_t when sigma=1,
  real beta;
}
transformed parameters {
  vector[N] mu_y;
  real<lower=0> sigma; 
  real<lower=2> nu; 
  
  sigma = theta / (kappa+1);
  nu = 2 * square(kappa+1) / (square(kappa+1) - 1);
  mu_y = beta * x;
}
model {
  theta ~ normal(0, 1);              
  kappa ~ lognormal(-2.68, 0.68);    // Very roughly similar to: nu ~ gamma(2, 0.1); do prior simulations
  beta ~ normal(0, 0.5);             // Slope for linear model
  y ~ student_t(nu, mu_y, sigma);    // Likelihood
}
