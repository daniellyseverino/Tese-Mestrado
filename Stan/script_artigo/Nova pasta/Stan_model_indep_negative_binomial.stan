functions {
  /**
    * Return the log probability of a proper conditional autoregressive (CAR) prior
  * with a sparse representation for the adjacency matrix
  *
    * @param param Vector containing the parameters with a CAR prior
  * @param tau Precision parameter for the CAR prior (real)
  * @param phi Dependence (usually spatial) parameter for the CAR prior (real)
  * @param W_sparse Sparse representation of adjacency matrix (int array)
  * @param n Length of phi (int)
  * @param W_n Number of adjacent pairs (int)
  * @param D_sparse Number of neighbors for each location (vector)
  * @param lambda Eigenvalues of D^{-1/2}*W*D^{-1/2} (vector)
  *
    * @return Log probability density of CAR prior up to additive constant
  */
    
    real normal_indep_lpdf(vector param, real tau, int n) {
      
      
      return 0.5*(n*log(tau)-tau*dot_product(param,param));
    }
  
}


//
  // Stan model to evaluated the cases of Covid-19 - Poisson model
// model: generalized static logistics


data {
  
  //-----------------------------
    // observed data
  int<lower=1> n; // number of states
  int<lower=1> t; // number of times
  int<lower=0> y[n,t]; // counts of new cases
  vector<lower=0>[n] pop; //one pop per state
  real<lower=0,upper=1> p; // Max proportion of infected
  //-----------------------------
}



parameters {
  real f1;
  vector[n] b1;
  vector[n] c1;
  real <lower=0> tau2_b;
  real<lower=0> tau2_c;
   real<lower=0> phi;
  vector<lower=0, upper=1>[n] a2;
}

transformed parameters{
  
  real f;
  vector[n] b;
  vector<lower=0>[n] c;
  matrix[n, t] mu;
  //  vector[n] a =  exp(-500  + ((log(p)+( exp(f1) * b1) ) + 500 ) .* a2);
  vector[n] a =  exp(-25  + ((log(p)+( exp(f1) * b1) ) + 25 ) .* a2);
  vector[n] pop2;
  vector[n] a1;
  
  
  
  for (i in 1:n) {
    b[i] = exp(b1[i]);
    c[i] = exp(c1[i]);
    
    
    
    pop2[i] = pop[i]/1000000;
    
    a1[i] = a[i]*1000000;
    
  }
  f = exp(f1);
  
  for(i in 1:n){
    for(j in 1:t){
      mu[i,j] = pop2[i]*exp(f1+log(a1[i])+log(c[i])-(c[i]*j)-(f+1)*log( b[i] +exp(-c[i]*j) ) );
    }
  }
  
}


model {
  //----------------------------
    // likelihood function
  for (i in 1:n) {
    row_vector[t] mu_aux;
    mu_aux = mu[i,1:t];
    y[i,1:t] ~ neg_binomial_2(mu_aux,phi); // observed model
  }
  //----------------------
    // prior distributions
  c1 ~ normal(0,1/sqrt(tau2_c));
  b1 ~ normal(0,1/sqrt(tau2_b));
  tau2_b ~ gamma(0.01,0.01);
  tau2_c ~ gamma(0.01,0.01);
  
  f1 ~ normal(log(101)/2 , sqrt(log(101))); // Mean= 1, var =100
  for(i in 1:n){
    //a2[i] ~ beta(0.01,0.01);
    a2[i] ~ uniform(0,1);
  }
  
}



