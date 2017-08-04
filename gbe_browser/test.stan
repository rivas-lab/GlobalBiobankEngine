  data {
    int<lower=1> N;
    real y[N];
  } 
  parameters {
    real mu;
  } 
  model {
    target += normal_lpdf(y | mu, 1);
  } 
