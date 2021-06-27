// Beverton-Holt growth model

data{
  int<lower = 1> N;
  int Fecundity[N];
  vector[N] intra;
  real intra_g;
  vector[N] avfa;
  vector[N] brho;
  vector[N] esca;
  vector[N] laca;
  vector[N] vumy;
  vector[N] trhi;
}

parameters{
  //real<lower = 0> lambda;
  real<lower = 0, upper = 10000> lambda;
  real alpha_avfa;
  real alpha_brho;
  real alpha_esca;
  real alpha_laca;
  real alpha_vumy;
  real alpha_trhi;

  //real<lower = 0>  alpha_avfa;
  //real<lower = 0>  alpha_brho;
  //real<lower = 0>  alpha_esca;
  //real<lower = 0>  alpha_laca;
  //real<lower = 0>  alpha_trhi;
  //real<lower = 0>  alpha_vumy;

}

model{
  // create a vector of predictions
  vector[N] F_hat;

  // set priors
  alpha_avfa ~ normal(0, 1);
  alpha_brho ~ normal(0, 1);
  alpha_esca ~ normal(0, 1);
  alpha_laca ~ normal(0, 1);
  alpha_vumy ~ normal(0, 1);
  alpha_trhi ~ normal(0, 1);
  lambda ~ gamma(0.001, 0.001);


  // implement the biological model
  for(i in 1:N){
    F_hat[i] = lambda*intra[i]*intra_g / (1 + alpha_avfa*avfa[i]*0.78 + alpha_brho*brho[i]*0.8 + alpha_esca*esca[i]*0.95 + alpha_laca*laca[i]*0.8 + alpha_vumy*vumy[i]*0.6 + alpha_trhi*trhi[i]*0.2);
  }

  // calculate the likelihood
  Fecundity ~ poisson(F_hat);
}
