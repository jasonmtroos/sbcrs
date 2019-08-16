data {
  int<lower = 0> n_obs;
  int<lower = 0> n_groups;
  int<lower = 0> n_types;
  int<lower = 1> group[n_obs];
  int<lower = 1> type[n_obs];
  vector[n_obs] y;
}
parameters {
  matrix[n_groups, n_types] mu;
  vector<lower = 0>[n_types] sigma;
  real nu;
}
model {
  sigma ~ exponential(1); 
  for (g in 1:n_groups) {
    mu[g,] ~ normal(0, 1);
  }
  for (o in 1:n_obs) {
    int g;
    int t;
    g = group[o];
    t = type[o];
    y[o] ~ normal(mu[g, t], sigma[t]);
  }
  nu ~ normal(0, 1); // just a nuisance scalar parameter
}
