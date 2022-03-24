
    data {
      int <lower=1> N;
      vector[N] x;
      vector[N] y;
    }
    parameters {
      real <lower=0> alpha;
      real <lower=0> beta;
      real <lower=0, upper=1> gamma;
      real <lower=0> sigma;
    }
    model {
      vector[N] mu;
      gamma ~ beta(7,7.0/3.0);
      for (i in 1:N){
         mu[i] = alpha - beta * pow(gamma,(x[i]));
         y[i] ~ normal(mu[i], sigma);
      }
      
    }

