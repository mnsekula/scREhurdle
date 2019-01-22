data {
  int<lower=0> G;               // number of genes
  int<lower=0> N;               // number of cells
  int<lower=1> M;               // number of columns in model matrix X + 1 to include random effects
  int<lower=0> y[G,N];          // expression matrix
  matrix[N, M-1] X;             // model matrix
  vector<lower=0>[M] prior_s;   // scale for cauchy priors (betas/zetas)
}// end data

parameters {
  matrix[M,G] beta_L;           // logistic betas: beta_L[M,] = zeta_L
  matrix[M,G] beta_C;           // t neg binom betas: beta_C[M,] = zeta_C
  vector[N] omega;              // independent random effects
  vector<lower=0>[G] phi;       // overdispersion parameter  
  real lambda1;                 // hyperparameter for phi
  real<lower=0> lambda2;        // hyperparameter for phi         
  real<lower=0> sigma2;         // variance of random effects
}// end parameters

model {
  matrix[N,G] X_beta_L;
  matrix[N,G] X_beta_C;
  matrix[N,M] X_full;
  matrix[N,G] zero_adj;
  X_full = append_col(X,omega);
  
  X_beta_L = X_full*beta_L;
  X_beta_C = X_full*beta_C;
  
  for(g in 1:G){
    for(i in 1:N){
      if (y[g,i] == 0) 
        1 ~ bernoulli_logit(X_beta_L[i,g]);
      else{
        zero_adj[i,g] = neg_binomial_2_log_lpmf(0|X_beta_C[i,g], phi[g]);
        //zero_adj[i,g] = (phi[g]/(exp(X_beta_C[i,g])+phi[g]))^phi[g];
        0 ~ bernoulli_logit(X_beta_L[i,g]);
        y[g,i] ~ neg_binomial_2_log(X_beta_C[i,g], phi[g]);
        target += - log1m_exp(zero_adj[i,g]);
        //target += - log1m(zero_adj[i,g]);
      }
    }// end cell loop
  }// end gene loop
    
  sigma2 ~ inv_gamma(1,1);
  omega ~ normal(0,sqrt(sigma2));

  for(j in 1:M){
    beta_L[j,] ~ cauchy(0,prior_s[j]);
    beta_C[j,] ~ cauchy(0,prior_s[j]);
  }
    
  lambda1 ~ cauchy(0,10);
  lambda2 ~ inv_gamma(0.001,0.001);
  phi ~ lognormal(lambda1,lambda2);
}// end model

generated quantities {
  row_vector[G] beta_L1;        // beta_L corresponding to treatment effect
  row_vector[G] beta_C1;        // beta_L corresponding to treatment effect
  row_vector[G] zeta_L;         // direct access zeta_L (beta_L[M,] = zeta_L)
  row_vector[G] zeta_C;         // direct access zeta_C (beta_C[M,] = zeta_C)
  beta_L1 = beta_L[2,];
  beta_C1 = beta_C[2,];
  zeta_L = beta_L[M,];
  zeta_C = beta_C[M,];
}// end gq

