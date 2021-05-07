functions {
  vector additive_trend(int n, vector innovations, real scale) {
    vector[n] trend;
        
    trend = cumulative_sum(innovations);
    return (trend - mean(trend)) * scale;
  }  
} 
data {
  int<lower=1> N;    
  int<lower=1> M;                                   // Number of predictors (on X)
  int<lower=1> N_p;                                 // Number of partitions (which is n_knots-1)
  real y[N];                                        // Outcome
  matrix[N, M] X;                                   // Predictors modeled with the piecewise linear function
}
transformed data {
  matrix[N_p+1, M] knots;                 
  int<lower=1,upper=N_p> partition[N, M];           // Index for the parameters to be weighted
  int<lower=2,upper=N_p+1> partition_plus_1[N, M];  
  real n_part_per_sd[M];  
  real partition_size[M];  
  matrix<lower=0,upper=1>[N, M] pct_p1; 
  matrix<lower=0,upper=1>[N, M] pct_p2; 
  
  for (j in 1:M) {
    partition_size[j] = (max(X[ , j]) - min(X[ , j]) + 0.002) / N_p;  // 0.002 is for the floating-point error
    knots[ , j] =  rep_vector(min(X[ , j]) - 0.001, N_p+1);
    for (p in 2:(N_p+1)) {  
      knots[p, j] += partition_size[j] * (p-1);
    }
    {
      int last_k;
      for (i in 1:N) {
        last_k = 2;
        while (X[i, j] > knots[last_k, j]) {
          last_k += 1;
        }
        partition[i, j] = last_k - 1;
        partition_plus_1[i, j] = last_k;
      }
    }
    pct_p1[ , j] = ((X[ , j] - knots[partition[ , j], j]
                    ) / partition_size[j]); 
    pct_p2[ , j] = 1 - pct_p1[ , j];
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
  sigma_partition ~ normal(0, 0.5);              // This is (mostly) invariant to the number of partitions
  for (j in 1:M) {
    epsilon_partition[ , j] ~ std_normal();      // For the non-centered parameterization 
  }
  y ~ normal(mu_y, sigma);                       // Likelihood 
}
