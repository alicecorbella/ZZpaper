
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
      mu = alpha + beta * pow(gamma,(x));
      alpha ~ uniform(0,100000);
      beta  ~ uniform(0,100000);
      gamma ~ beta(7,7/3);
      sigma ~ uniform(0,100000);
      y ~ normal(mu, sigma);
    }

