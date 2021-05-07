functions {
  vector additive_trend(int n, vector innovations, real scale) {
    vector[n] trend;
        
    trend = cumulative_sum(innovations);
    return (trend - mean(trend)) * scale;
  }  
} 
data {
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
    pct_p1[ , j] = ((X[ , j] - partition_knots[partition[ , j], j]
                    ) / partition_size[j]); 
    pct_p2[ , j] = 1 - pct_p1[ , j];
    for (i in 1:N)
        partition_plus_1[i, j] = partition[i, j] + 1;
    n_part_per_sd[j] = sd(X[ , j]) / partition_size[j];      
  }
}
parameters {
  real alpha;
  real<lower=0> sigma; 
  real<lower=0> sigma_partition[M];              
  matrix[N_p+1, M] epsilon_partition;
}
transformed parameters {
  vector[N] mu_y;
  matrix[N_p+1, M] phi_partition;
  
  mu_y = rep_vector(alpha, N);
  for (j in 1:M) {
    phi_partition[ , j] = additive_trend(N_p+1, 
                                         epsilon_partition[ , j],
                                         sigma_partition[j]/sqrt(n_part_per_sd[j]));
    mu_y += (phi_partition[partition[ , j], j] .* pct_p2[ , j] 
               + phi_partition[partition_plus_1[ , j], j] .* pct_p1[ , j]);
  }
}
model {
  alpha ~ normal(0, 0.5);                        // Constant
  sigma ~ normal(0, 0.5);
  sigma_partition ~ normal(0, 0.5);              // This is (approximately) invariant to the number of partitions
  for (j in 1:M) {
    epsilon_partition[ , j] ~ std_normal();
  }
  y ~ normal(mu_y, sigma);                       // Likelihood 
}

