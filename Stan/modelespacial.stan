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
  int<lower = 1> S;
  int<lower = 1> T;
  int<lower = 1> D;
  int<lower=0> n[S, T, D - 1];
  int<lower=0> N[S, T];
  
}

parameters {

  vector<lower = 0>[S] a_theta;
  vector<lower = 0>[S] b_theta;
  vector<lower = 0>[S] c_theta;
  vector<lower = 0>[S] f_theta;
  vector<lower = 0>[S] a_alpha;
  vector<lower = 0>[S] b_alpha;
  vector<lower = 0>[S] c_alpha;
  vector<lower = 0>[S] f_alpha;
  matrix[S, D - 1] beta;
  vector<lower = 0>[S] phi_n;
  vector<lower = 0>[S] phi_N;
  
}

transformed parameters {
  
  array[S] matrix<lower = 0>[T, D - 1] lambda;
  matrix[S,T] alpha;
  real<lower = sum(lambda[S,T, ])> theta[S,T];

  real b1_theta = log(b_theta);  
  real b1_alpha = log(b_alpha);  

  for(s in 1:S){
    for(t in 1:T){
      alpha[s,t] = genLog(t, a_alpha[s], b_alpha[s], c_alpha[s], f_alpha[s], 1);
      theta[s,t] = genLog(t, a_theta[s], b_theta[s], c_theta[s], f_theta[s], 0);
    }
  }
  
  for(s in 1:S){
    for(t in 1:T){
      for(d in 1:(D - 1)){
        lambda[s][t, d] = exp(alpha[s,t] + beta[s,d]);
      }
    }
  }

  
}


model {
  
  // -----> likelihood function
  
  for(s in 1:S){
    for(t in 1:T){
      for(d in 1:(D - 1)){
        n[s, t, d] ~ neg_binomial_2( lambda[s][t, d], phi_n[s] );
      }
    }
  }

  for(s in 1:S){
    for(t in 1:T){
      N[s,t] ~ neg_binomial_2(theta[s,t], phi_N[s]);
    }
  }
  
  // -----> prior distributions
  a_theta[s] ~ gamma(0.1, 0.1);
  b1_theta[s] ~ normal(0, sqrt(20));  // sqrt(1/0.2)
  c_theta[s] ~ gamma(2,9); 
  f_theta[s] ~ gamma(0.01,0.01);
  a_alpha[s] ~ gamma(0.1, 0.1);
  b1_alpha[s] ~ normal(0, sqrt(20));  // sqrt(1/0.2)
  c_alpha[s] ~ gamma(2,9);           //  gamma(2,9)  shape=2, scale=9
  f_alpha[s] ~ gamma(0.01,0.01);
  beta ~ normal(0, 100);
  phi_n[s] ~ gamma(4,0.01);
  phi_N[s] ~ gamma(4,0.01);
  
}
