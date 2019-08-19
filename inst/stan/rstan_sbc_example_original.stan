data {
  int<lower = 1> N;
  real<lower = 0> a;
  real<lower = 0> b;
  int<lower = 0, upper = N> y;
}
parameters {
  real<lower = 0, upper = 1> pi;
}
model {
  target += beta_lpdf(pi | a, b);
  target += binomial_lpmf(y | N, pi);
}
generated quantities {
  real pi_ = beta_rng(a, b);

  vector[N] log_lik;
  for (n in 1:y) log_lik[n] = bernoulli_lpmf(1 | pi);
  for (n in (y + 1):N) log_lik[n] = bernoulli_lpmf(0 | pi);
}
