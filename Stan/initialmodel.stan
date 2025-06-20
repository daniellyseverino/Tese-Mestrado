// Stan Model with missing counts

functions {

  real genLog(int t, real a, real b, real c, real f, int logScale) {
    real logV = log(f)+log(a)+log(c)-(c*t)-(f+1)*log( b+exp(-c*t) );
    if (logScale == 1){
      return logV;
    } else {
      return exp(logV);
    }
  }

}


data {
  
  // -----> observed data
  int<lower = 1> T;
  int<lower = 1> D;
  int<lower=0> n[T, D - 1];
  int<lower=0> N[T];
  
}

parameters {

  real<lower = 0> a_theta;
  real<lower = 0> b_theta;
  real<lower = 0> c_theta;
  real<lower = 0> f_theta;
  real<lower = 0> a_alpha;
  real<lower = 0> b_alpha;
  real<lower = 0> c_alpha;
  real<lower = 0> f_alpha;
  vector[D - 1] beta;
  
}

transformed parameters {
  
  matrix<lower = 0>[T, D - 1] lambda;
  vector[T] alpha;
  real<lower = sum(lambda[T, ])> theta[T];

  real b1_theta = log(b_theta);  
  real b1_alpha = log(b_alpha);  

  for(t in 1:T){
    alpha[t] = genLog(t, a_alpha, b_alpha, c_alpha, f_alpha, 1);
  }
  
  for(t in 1:T){
    for(d in 1:(D - 1)){
      lambda[t, d] = exp(alpha[t] + beta[d]);
    }
  }

  for(t in 1:T){
    theta[t] = genLog(t, a_theta, b_theta, c_theta, f_theta, 0);
  }
  
}


model {
  
  // -----> likelihood function
  
  for(t in 1:T){
    for(d in 1:(D - 1)){
      n[t, d] ~ poisson( lambda[t, d] );
    }
  }

  for(t in 1:T){
    N[t] ~ poisson(theta[t]);
  }
  
  // -----> prior distributions
  a_theta ~ gamma(0.1, 0.1);
  b1_theta ~ normal(0, sqrt(20));  // sqrt(1/0.2)
  c_theta ~ gamma(2,9); 
  f_theta ~ gamma(0.01,0.01);
  a_alpha ~ gamma(0.1, 0.1);
  b1_alpha ~ normal(0, sqrt(20));  // sqrt(1/0.2)
  c_alpha ~ gamma(2,9);           //  gamma(2,9)  shape=2, scale=9
  f_alpha ~ gamma(0.01,0.01);
  beta ~ normal(0, 100);
  
}
