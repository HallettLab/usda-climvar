// Model for within year fecundity of Bromus

data{
  int<lower = 1> N;
  int Fecundity[N];
  vector[N] intra;
  vector[N] avfa;
  vector[N] brho;
  vector[N] laca;
  vector[N] vumy;
  //  vector[N] dry; 
  // vector[N] wet;
}

parameters{
  real<lower = 0> lambda;
  //real<lower = 0>  alpha_avfa;
  //real<lower = 0>  alpha_brho;
  //real<lower = 0>  alpha_laca;
  //real<lower = 0>  alpha_vumy;

  real alpha_avfa;
  real alpha_brho;
  real alpha_laca;
  real alpha_vumy;

}

model{
  // create a vector of predictions
  vector[N] F_hat;

  // set priors
  alpha_avfa ~ normal(0, 1000);
  alpha_brho ~ normal(0, 1000);
  alpha_laca ~ normal(0, 1000);
  alpha_vumy ~ normal(0, 1000);
  lambda ~ gamma(0.001, 0.001);


  // implement the biological model
  for(i in 1:N){
    F_hat[i] = lambda*intra[i] / (1 + exp(alpha_avfa)*avfa[i] + exp(alpha_brho)*brho[i] + exp(alpha_laca)*laca[i] + exp(alpha_vumy)*vumy[i]);
  }

  // calculate the likelihood
  Fecundity ~ poisson(F_hat);
}





