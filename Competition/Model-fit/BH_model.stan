// Model for within year fecundity of Bromus

data{
  int<lower = 1> N;
  int Fecundity[N];
  vector[N] intra;
  vector[N] avfa;
  vector[N] brho;
  vector[N] esca;
  vector[N] laca;
  vector[N] trhi;
  vector[N] vumy;
  //  vector[N] dry; 
  // vector[N] wet;
}

parameters{
  real<lower = 0> lambda;
  real alpha_avfa;
  real alpha_brho;
  real alpha_esca;
  real  alpha_laca;
  real alpha_trhi;
  real alpha_vumy;

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
  alpha_avfa ~ normal(0, 1000);
  alpha_brho ~ normal(0, 1000);
  alpha_esca ~ normal(0, 1000);
  alpha_laca ~ normal(0, 1000);
  alpha_trhi ~ normal(0, 1000);
  alpha_vumy ~ normal(0, 1000);
  lambda ~ gamma(0.001, 0.001);


  // implement the biological model
  for(i in 1:N){
    F_hat[i] = lambda*intra[i] / (1 + alpha_avfa*avfa[i] + alpha_brho*brho[i]+ alpha_esca*esca[i]+ alpha_laca*laca[i]+ alpha_trhi*trhi[i]+ alpha_vumy*vumy[i]);
  }

  // calculate the likelihood
  Fecundity ~ poisson(F_hat);
}





