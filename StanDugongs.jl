
using StanSample
using Distributions
using DataFrames
using MonteCarloMeasurements, AxisKeys
using StanSample


CODEwd = "/Users/alice/projects/ZZpaper/"
SAVEwd = "/Users/alice/OneDrive - University of Warwick/Manuscripts/04ADZZ/figures/Sec5.1/"
NameEx = ""
# Loading the packages
include(string(CODEwd, "env.jl"))
include(string(CODEwd, "plotting.jl"))

pt(78)
gr(size = (pt(84), pt(84)), labelfontsize=8, legend = false)

r_insp=10
r_ess=100
budget=200000

Dim=4
data_x = [1.0, 1.5, 1.5, 1.5, 2.5, 4.0, 5.0, 5.0, 7.0, 8.0, 8.5, 9.0, 9.5,
    9.5, 10, 12, 12, 13, 13, 14.5, 15.5, 15.5, 16.5, 17, 22.5, 29.0, 31.5]
data_y = [1.80, 1.85, 1.87, 1.77, 2.02, 2.27, 2.15, 2.26, 2.47, 2.19, 2.26,
    2.40, 2.39, 2.41, 2.50, 2.32, 2.32, 2.43, 2.47, 2.56, 2.65, 2.47, 2.64,
    2.56, 2.70, 2.72, 2.57]


f1=plot(data_x, data_y, seriestype=:scatter, legend=false, color=:black,
    alpha=0.5, xlabel="Zⱼ", ylabel="Yⱼ")
dimnames=["α", "β", "γ", "σ"]



StanDugongs= "
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
      mu = alpha - beta * gamma^(x);
      gamma ~ beta(7,7/3);
      y ~ normal(mu, sigma);
    }
";


m6_1s = SampleModel("m6.1s", StanDugongs);

data = (y = data_y, x=data_x, N = length(data_x))


init = (alpha = exp(1), beta = exp(1), gamma= exp(1)/(1+exp(1)), sigma = exp(1))


rc6_1s = stan_sample(m6_1s; data, init, chains=1, warmup=0, thin=1);


if success(rc6_1s)
    st6_1s = read_samples(m6_1s) # By default a StanTable object is returned

    # Display the schema of the tbl

    st6_1s |> display
    println()

    # Display the draws

    df6_1s = DataFrame(st6_1s)
    df6_1s |> display
    println()

    # Or using a KeyedArray object from AxisKeys.jl

    chns6_1s = read_samples(m6_1s, :keyedarray)
    chns6_1s |> display
end

df6_1s = DataFrame(st6_1s)

plot(plot(df6_1s[:,1]), plot(df6_1s[:,2]), plot(df6_1s[:,3]), plot(df6_1s[:,4]),
    layout=(2,2))


plot(density(df6_1s[100:1000,1]), density(df6_1s[100:1000,2]),
    density(df6_1s[100:1000,3]), density(df6_1s[100:1000,4]),
    layout=(2,2), size=(600,600))

c=8
init = (alpha = exp(c), beta = exp(c), gamma= exp(c)/(1+exp(c)), sigma = exp(c))


rc6_1s = stan_sample(m6_1s; data, init);


if success(rc6_1s)
    st6_1s = read_samples(m6_1s) # By default a StanTable object is returned

    # Display the schema of the tbl

    st6_1s |> display
    println()

    # Display the draws

    df6_1s = DataFrame(st6_1s)
    df6_1s |> display
    println()

    # Or using a KeyedArray object from AxisKeys.jl

    chns6_1s = read_samples(m6_1s, :keyedarray)
    chns6_1s |> display
end

df6_1s = DataFrame(st6_1s)

plot(plot(df6_1s[:,1]),plot(df6_1s[:,2]),
    plot(df6_1s[:,3]),plot(df6_1s[:,4]), layout=(2,2),size=(600,600))
