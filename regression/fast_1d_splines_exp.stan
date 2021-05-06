functions {
  vector additive_trend(int n, vector innovations, vector scale) {
    vector[n] trend;
        
    trend = cumulative_sum(innovations);
    return (trend - mean(trend)) .* scale;
  }  
  vector mult_trend(int n, vector innovations, real scale) {
    // Only difference is the 'exp'
    vector[n] trend;
        
    trend = cumulative_sum(innovations);
    return exp((trend - mean(trend)) * scale);
  }  
} 
data {
  // TODO: delete any code related to the student-t reparameterization
  int<lower=0> N;    
  int<lower=0> M;                                   // Number of predictors (on X)
  int<lower=0> N_p;                                 // Number of partitions 
  real y[N];                                        // Outcome 
  matrix[N, M] X;
  int<lower=1,upper=N_p> partition[N, M];           // Index for the day of year
  matrix[N_p+1, M] partition_knots;                 // These should be equally spaced
}
transformed data {
  int<lower=2,upper=N_p+1> partition_plus_1[N, M];
  real n_part_per_sd[M];  
  real partition_size[M];  
  matrix<lower=0,upper=1>[N, M] pct_p1; 
  matrix<lower=0,upper=1>[N, M] pct_p2; 
  
  for (j in 1:M) {
    partition_size[j] = (partition_knots[2, j] - partition_knots[1, j]);  
    pct_p1[ , j] = ((X[ , j] - partition_knots[partition[ , j], j]) 
                      / partition_size[j]); 
    pct_p2[ , j] = 1 - pct_p1[ , j];
    for (i in 1:N)
        partition_plus_1[i, j] = partition[i, j] + 1;
    n_part_per_sd[j] = sd(X[ , j]) / partition_size[j];      
  }
}
parameters {
  real alpha;
  real<lower=0> theta;                              // Parameter for student_t reparameterization 
  real<lower=0> kappa;                              // Parameter for student_t reparameterization 
  real<lower=0> sigma_partition[M];                 // Scale for seasonal component 
  real<lower=0> sigma_partition_e[M];               // Scale for errors of partition
  matrix[N_p+1, M] epsilon_partition;
  matrix[N_p+1, M] epsilon_partition_e;
}
transformed parameters {
  vector[N] mu_y;
  matrix[N_p+1, M] phi_partition;
  matrix[N_p+1, M] phi_partition_e;
  real<lower=0> sigma; 
  real<lower=1> nu; 
  
  mu_y = rep_vector(alpha, N);
  for (j in 1:M) {
    phi_partition_e[ , j] = mult_trend(N_p+1, epsilon_partition_e[ , j], 
                                       sigma_partition_e[j]/sqrt(n_part_per_sd[j]));
    {
      vector[N_p+1] tr_sigma_partition;
      real adj_sigma_partition = sigma_partition[j] / sqrt(n_part_per_sd[j]);
      
      tr_sigma_partition = (rep_vector(adj_sigma_partition, N_p+1) 
                              .* phi_partition_e[ , j]);
      phi_partition[ , j] = additive_trend(N_p+1, 
                                           epsilon_partition[ , j],
                                           tr_sigma_partition);
    }
    mu_y += (phi_partition[partition[ , j], j] .* pct_p2[ , j] 
               + phi_partition[partition_plus_1[ , j], j] .* pct_p1[ , j]);
  }
  sigma = theta / (kappa+1);
  nu = 2 * square(kappa+1) / (square(kappa+1) - 1);
}
model {
  alpha ~ normal(0, 1);                             // Constant
  theta ~ normal(0, 1);                             // The standard deviation of the residuals
  kappa ~ normal(-2.68, 0.68);                      // Identical to setting nu ~ gamma(2, 0.1)  
  sigma_partition ~ normal(0, 0.5);                 // This is (mostly) invariant to the number of partitions
  sigma_partition_e ~ normal(0, 0.2); 
  for (j in 1:M) {
    epsilon_partition[ , j] ~ std_normal();
    epsilon_partition_e[ , j] ~ std_normal();
    phi_partition[ , j] ~ normal(0, 5);
    phi_partition_e[ , j] ~ normal(0, 5);
  }
  y ~ student_t(nu, mu_y, sigma);               // Likelihood 
}
