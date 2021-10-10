// Improvements for the future: have phi also be a vector; this allows me to have variables
// with different levels of knots (which I could maybe input in the data block)
functions {
  vector additive_trend(int n, vector innovations, real scale) {
    vector[n] trend;

    trend = cumulative_sum(innovations);
    return (trend - mean(trend)) * scale;
  }
}
data {
  int<lower=0> N;                                // Number of observations
  int<lower=0> M_h;                              // Number of variables to be modelled hierarchically (assuming exchangeability)
  int<lower=0> M_s;                              // Number of variables to be modelled with first-degree splines
  int<lower=1> indexes_h[N, M_h];                // Indexes
  int<lower=0> data_s[N, M_s];                   // Data (integers)
  int<lower=0,upper=1> missing_h[N, M_h];        // When 'indexes_h' has a missing observation (1) or not (0)
  int<lower=0,upper=1> missing_s[N, M_s];        // When 'data_s' has a missing observation (1) or not (0)
  int<lower=0,upper=1> y[N];                     // Outcome 
}
transformed data {
  int<lower=1> N_p = 50;                         // Number of partitions (which is n_knots-1) when data_as_index = 0
  matrix[N_p+1, M_s] knots;
  int<lower=1,upper=N_p> partition[N, M_s];      // Index for the parameters to be weighted
  int<lower=2,upper=N_p+1> partition_plus_1[N, M_s];
  real partition_size[M_s];
  real n_part_per_sd[M_s];
  matrix<lower=-0.001,upper=1.001>[N, M_s] pct_p1;    // Upper/lower limit is due to floating-point error
  matrix<lower=-0.001,upper=1.001>[N, M_s] pct_p2;
  matrix[N, M_s] X;
  int<lower=0,upper=1> has_miss_h[M_h];
  int<lower=0,upper=1> has_miss_s[M_s];
  int<lower=1> indexes_h_offset[N, M_h] = indexes_h;                // Indexes
  int<lower=1> N_beta_h[M_h];

  for (j in 1:M_h) {
    N_beta_h[j] = max(indexes_h[ , j]);
    if (j > 1) {
      for (i in 1:N) {
        indexes_h_offset[i, j] += sum(N_beta_h[1:j-1]);
      }
    }
    has_miss_h[j] = sum(missing_h[ , j]) >= 1;
  }
  for (j in 1:M_s) {
    X[ , j] = to_vector(data_s[ , j]);
    has_miss_s[j] = sum(missing_s[ , j]) >= 1;
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
  vector[sum(N_beta_h)] beta_raw_h;         // Raw parameter for hierarchical model
  vector[sum(has_miss_h)] nu_h;             // Parameters for missing variables
  vector[sum(has_miss_s)] nu_s;             // Parameters for missing variables
  real<lower=0> sigma_h[M_h];
  real<lower=0> sigma_s[M_s];
  matrix[N_p+1, M_s] epsilon_s;
}
transformed parameters {
  vector[N] theta_hat;
  matrix[N_p+1, M_s] phi_s;
  vector[sum(N_beta_h)] beta_h;             // beta_raw_h * sigma_h

  theta_hat = rep_vector(alpha, N);
  {
    int offset_nu_s = 0;
    int offset_nu_h = 0;
    int offset_beta_h = 0;
    vector[N] theta_hat_s;
    vector[N] theta_hat_h;

    for (j in 1:M_s) {
      phi_s[ , j] = additive_trend(N_p+1,
                                   epsilon_s[ , j],
                                   sigma_s[j] / sqrt(n_part_per_sd[j]));
      theta_hat_s = (phi_s[partition[ , j], j] .* pct_p2[ , j]
                     + phi_s[partition_plus_1[ , j], j] .* pct_p1[ , j]);
      if (has_miss_s[j]) {
        for (i in 1:N) {
          if (missing_s[i, j]) {
              theta_hat_s[i] = nu_s[j - offset_nu_s];
          }
        }
      }
      offset_nu_s += -(has_miss_s[j] - 1);
      theta_hat += theta_hat_s;
    }
    for (j in 1:M_h) {
      if (j == 1) {
        beta_h[1:N_beta_h[j]] = beta_raw_h[1:N_beta_h[j]] * sigma_h[j];
      }
      else {
        {
          int prev_max = sum(N_beta_h[1:(j-1)]);
          beta_h[(prev_max+1):(prev_max+N_beta_h[j])] = beta_raw_h[(prev_max+1):(prev_max+N_beta_h[j])] * sigma_h[j];
        }
      }
      theta_hat_h = beta_h[indexes_h_offset[ , j]];
      if (has_miss_h[j]) {
        for (i in 1:N) {
          if (missing_h[i, j]) {
            theta_hat_h[i] = nu_h[j - offset_nu_h];
          }
        }
      }
      offset_nu_h += -(has_miss_h[j] - 1);
      theta_hat += theta_hat_h;
    }
  }
}
model {
  alpha ~ normal(-3, 3);                        // Expected value of ~5%
  beta_raw_h ~ std_normal();
  sigma_h ~ exponential(1/0.10);
  sigma_s ~ exponential(1/0.10);
  nu_h ~ normal(0, 0.05);
  nu_s ~ normal(0, 0.05);
  for (j in 1:M_s) {
    epsilon_s[ , j] ~ std_normal();
  }
  y ~ bernoulli_logit(theta_hat);
}
