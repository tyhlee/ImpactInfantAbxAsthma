data {
  int<lower=1> N;                    //number of data points
  int<lower=1> J;                    //number of unique studies
  int<lower=1> K;                   // number of betas
  int<lower=1, upper=J> studyid[N];  //study id
  real log_aOR[N];                   //adjusted aOR
  real log_aOR_se[N];                //adjusted aOR se
  int<lower=0,upper=5> dose[N];      // dose
  real<lower=0> t[N];                // asthma Dx age
  int<lower=1> N_test;
  matrix[N_test,2] x_test;
}

parameters {
  vector[K] beta;
  vector[J] u;               
  real<lower=0>sigma_u;      
  cholesky_factor_corr[K] Lcorr;
  vector<lower=0>[K] sigma;
}

transformed parameters{
  vector[N] log_mu;
  corr_matrix[K] R; // correlation matrix
  cov_matrix[K] Sigma; // V * L * V
  R = multiply_lower_tri_self_transpose(Lcorr); // R = Lcorr * Lcorr'
  Sigma = quad_form_diag(R, sigma);
  for (i in 1:N){
    log_mu[i] = beta[1] + beta[2]*t[i] + beta[3]*dose[i]  +  u[studyid[i]];
  }
}

model {
  //priors
  u ~ normal(0,sigma_u);    
  sigma_u ~ lognormal(0,10);
  sigma ~ cauchy(0, 5); 
  Lcorr ~ lkj_corr_cholesky(1.0);
  beta ~ multi_normal([0,0,0],Sigma);
  //likelihood
  for (i in 1:N){
    log_aOR[i] ~ normal(log_mu[i],log_aOR_se[i]);
  }
}

generated quantities{
  vector[N] y_pred;
  vector[N_test] y_test;
  vector[N_test] y_patrick;
  
  for(i in 1:N){
    y_pred[i] = normal_rng(log_mu[i],log_aOR_se[i]);
  }
  
  for(j in 1:N_test){
    y_test[j] = normal_rng(beta[1] + beta[3]*x_test[j,2]+beta[2]*x_test[j,1],sigma_u);
  }
  
  for(k in 1:N_test){
    y_patrick[k] = beta[1] + beta[3]*x_test[k,2]+beta[2]*x_test[k,1]+u[studyid[1]];
  }
}
