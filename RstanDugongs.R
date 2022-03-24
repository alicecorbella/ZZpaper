library(rstan)
setwd( "/Users/alice/projects/ZZpaper")



Dim=4
data_x =c(1.0, 1.5, 1.5, 1.5, 2.5, 4.0, 5.0, 5.0, 7.0, 8.0, 8.5, 9.0, 9.5,
          9.5, 10, 12, 12, 13, 13, 14.5, 15.5, 15.5, 16.5, 17, 22.5, 29.0, 31.5)
data_y = c(1.80, 1.85, 1.87, 1.77, 2.02, 2.27, 2.15, 2.26, 2.47, 2.19, 2.26,
          2.40, 2.39, 2.41, 2.50, 2.32, 2.32, 2.43, 2.47, 2.56, 2.65, 2.47, 2.64,
          2.56, 2.70, 2.72, 2.57)

write( "
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
", "Dugongs.stan")

stanc("Dugongs.stan")

data = list(y = data_y, x=data_x, N = length(data_x))

model = stan_model("Dugongs.stan")


c=1
init = list(list(alpha = exp(c), beta = exp(c), gamma= exp(c)/(1+exp(c)), sigma = exp(c)))


fit = sampling(model,data=data,iter=20000,chains=1, init=init, warmup=1000, thin=1)

params = extract(fit)

par(mfrow=c(2,2))
plot(params$alpha, type="l")
plot(params$beta, type="l")
plot(params$gamma, type="l")
plot(params$sigma, type="l")



c=40
init = list(list(alpha = exp(c), beta = exp(c), gamma= exp(c)/(1+exp(c)), sigma = exp(c)), 
            list(alpha = exp(c), beta = exp(c), gamma= exp(c)/(1+exp(c)), sigma = exp(-c)),
            list(alpha = exp(sample(c(-1,1), 1)*c), beta = exp(sample(c(-1,1), 1)*c), gamma= exp(-c)/(1+exp(-c)), sigma = exp(sample(c(-1,1), 1)*c)), 
            list(alpha = exp(-c), beta = exp(-c), gamma= exp(-c)/(1+exp(-c)), sigma = exp(-c)))

print(init)

fit = sampling(model,data=data,iter=2000,chains=4, init=init, warmup=1000, thin=1)

params = extract(fit)

par(mfrow=c(2,2))
plot(params$alpha[1:1000], type="l", ylim=c(0, max(params$alpha)))
lines(params$alpha[1001:2000], col=2)
lines(params$alpha[2001:3000], col=3)
lines(params$alpha[3001:4000], col=4)
plot(params$beta[1:1000], type="l", ylim=c(0, max(params$beta)))
lines(params$beta[1001:2000], col=2)
lines(params$beta[2001:3000], col=3)
lines(params$beta[3001:4000], col=4)
plot(params$gamma[1:1000], type="l", ylim=c(0, max(params$gamma)))
lines(params$gamma[1001:2000], col=2)
lines(params$gamma[2001:3000], col=3)
lines(params$gamma[3001:4000], col=4)
plot(params$sigma[1:1000], type="l", ylim=c(0, max(params$sigma)))
lines(params$sigma[1001:2000], col=2)
lines(params$sigma[2001:3000], col=3)
lines(params$sigma[3001:4000], col=4)

