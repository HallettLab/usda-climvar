// Beverton-Holt growth model

data{
  int<lower = 1> N;
  int Fecundity[N];
  vector[N] trhi;
  vector[N] inter;
}

parameters{
  //real<lower = 0> lambda;
  real<lower = 0> lambda;
  real alpha_trhi;
  real alpha_inter;
}

model{
  // create a vector of predictions
  vector[N] F_hat;

  // set priors
  alpha_trhi ~ normal(0, 1);
  alpha_inter ~ normal(0, 1);
  //lambda ~ gamma(0.001, 0.001);
  lambda ~ normal(100, 10);


  // implement the biological model
  for(i in 1:N){
    F_hat[i] = lambda*trhi[i] / (1 + alpha_trhi*trhi[i] + alpha_inter*inter[i]);
  }

  // calculate the likelihood
  Fecundity ~ poisson(F_hat);
}
