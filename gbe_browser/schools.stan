data {
     int<lower=0> J;          // number of schools
     real betahat[J];               // estimated treatment effects
     real<lower=0> sigma[J];  // standard error of effect estimates
      }
parameters {
  real mu;
    real<lower=0> tau;
      vector[J] eta;
      }
transformed parameters {
  vector[J] beta;
    beta = mu + tau * eta;
    }

model {
  target += normal_lpdf(eta | 0, 1);
    target += normal_lpdf(betahat | beta, sigma);
    }
