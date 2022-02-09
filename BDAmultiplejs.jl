# Direcotries
CODEwd = "/Users/alice/projects/ZZpaper/"
SAVEwd = "/Users/alice/OneDrive - University of Warwick/Manuscripts/04ADZZ/figures/sec6.2/"

# Loading the packages
include(string(CODEwd, "env.jl"))

# Include supplementary packages for data IO
using CSV
using DataFrames
using StatsBase
using Serialization
using Dates
using Extremes

Y = open(deserialize, string(CODEwd,"data/BDA_preprocessed.bin"))

# a very small data for the plots
Yvsd = Y[sample(1:(size(Y)[1]), 100), :]
print(Yvsd)

# survival plot for
minimum(Yvsd.DIAGNOSISDATEBEST)

p1=plot([Yvsd.DIAGNOSISDATEBEST[1], Yvsd.VITALSTATUSDATE[1]], [1, 1],
    color=:gray, ylabel="Patients", xlabel="Calendar time")
for i in 1:100
    plot!([Yvsd.DIAGNOSISDATEBEST[i], Yvsd.VITALSTATUSDATE[i]], [i, i], color=:gray)
end
plot!(Yvsd.VITALSTATUSDATE,1:100, seriestype=:scatter,
    markershape=ifelse.(Yvsd.NEWVITALSTATUS.=="D", :star4, :circle),
        markersize=ifelse.(Yvsd.NEWVITALSTATUS.=="D",4, 2.5),
        color=ifelse.(Yvsd.NEWVITALSTATUS.=="D",:red, :green))
plot!([minimum(Yvsd.DIAGNOSISDATEBEST), minimum(Yvsd.DIAGNOSISDATEBEST)], [110,105],
    seriestype=:scatter, markershape=[:star4, :circle], color=[:red, :green])
annotate!([(minimum(Yvsd.DIAGNOSISDATEBEST), 110, ("  Dead", 10, :black, :left)),
        (minimum(Yvsd.DIAGNOSISDATEBEST), 105, ("  Censored", 10, :black, :left))])
savefig(p1, string(SAVEwd, "f1.pdf"))


# a small data to compare the subsampled version with the full-data version

Random.seed!(1234)
Ysd = Y[sample(1:(size(Y)[1]), 5000), :]
print(Ysd[1:20, :])

# model with only age as a covariate
t_vec = Dates.value.(Ysd.TIMEDIAGNOSITOFINALEVENT).+0.01
c_vec = zeros(length(t_vec))
c_vec[findall(Ysd.NEWVITALSTATUS.=="D")] .= 1
z1_vec = (Ysd.AGE.-mean(Ysd.AGE))./std(Ysd.AGE)

histogram(t_vec)
histogram(t_vec[findall(c_vec.==1)])
histogram(t_vec[findall(c_vec.==0)])
histogram(c_vec)


dimnames=["lgα", "β₀", "β₁"]

x = [-0.5, +6.5, -0.4]
t = t_vec
c = c_vec
z1= z1_vec
# transform parameter
α = exp(x[1])

J = length(t)
mll = 0
for j in 1:J
    lgμ_j = x[2] + x[3] *z1[j]
    μ_j   = exp(lgμ_j)
    if c[j]==1
        mll+= -log(α)+(1-α)*log(t[j])+α*log(μ_j)+((t[j]/μ_j)^α)
    elseif c[j]==0
        mll+= ((t[j]/μ_j)^α)
    end
    print(string(j, ", ", mll))
end

function U(x::Vector; t = t_vec, c = c_vec, z1= z1_vec)
    # transform parameter
    α = exp(x[1])
    J = length(t)
    mll = 0
    for j in 1:J
        lgμ_j = x[2] + x[3] *z1[j]
        μ_j   = exp(lgμ_j)
        if c[j]==1
            mll+= -log(α)+(1-α)*log(t[j])+α*log(μ_j)+((t[j]/μ_j)^α)
        elseif c[j]==0
            mll+= ((t[j]/μ_j)^α)
        end
    end
    return mll
end

U([-0.5, +6.5, -0.4])
U([-0.8, 5, -1])
U([0, 35, -2])
U([1, -1, 1])


# Loading the functions for general zig zag
include(string(CODEwd, "functions.jl"))
include(string(CODEwd, "plotting.jl"))

Dim=3
start= [-0.5, +6.5, -0.4]
zzsk =  zz(; NS=1000, x0_0=start, tmax=0.1)
zzsk =  zz(; NS=1000, x0_0=start, tmax=0.1)


zzsk["SK"]
p2a=plot(zzsk["SK"][:, 1], zzsk["SK"][:, 2], title=dimnames[1])
p2b=plot(zzsk["SK"][:, 1], zzsk["SK"][:, 3], title=dimnames[2])
p2c=plot(zzsk["SK"][:, 1], zzsk["SK"][:, 4], title=dimnames[3])

p2d=plot(zzsk["SK"][:, 2], zzsk["SK"][:, 3], xlabel=dimnames[1], ylabel=dimnames[2])
p2e=plot(zzsk["SK"][:, 2], zzsk["SK"][:, 4], xlabel=dimnames[1], ylabel=dimnames[3])
p2f=plot(zzsk["SK"][:, 3], zzsk["SK"][:, 4], xlabel=dimnames[2], ylabel=dimnames[3])


# f2 - tune of t_max
f2=tmplot(R=10, try_tmax=[0.01, 0.025, 0.05, 0.075, 0.1],
    start=[-0.5, 9, -0.5])
savefig(f2, string(SAVEwd, "f2.pdf"))
tmax_tuned=0.05


start1 = [-0.5, +6.5, -0.4]
start2 = [-0.8, 5, -1]
start3 = [0, 35, -2]
start4 = [1, -1, 1]

zzsk1 =  zz(; NS=1000, x0_0=start1, tmax=0.05)
zzsk2 =  zz(; NS=1000, x0_0=start2, tmax=0.05)
zzsk3 =  zz(; NS=1000, x0_0=start3, tmax=0.05)
zzsk4 =  zz(; NS=1000, x0_0=start4, tmax=0.05)


plot(zzsk1["SK"][:, 1], zzsk1["SK"][:, 2], title=dimnames[1])
plot!(zzsk2["SK"][:, 1], zzsk2["SK"][:, 2], title=dimnames[1])
plot!(zzsk3["SK"][:, 1], zzsk3["SK"][:, 2], title=dimnames[1])
plot!(zzsk4["SK"][:, 1], zzsk4["SK"][:, 2], title=dimnames[1])

plot(zzsk1["SK"][:, 1], zzsk1["SK"][:, 3], title=dimnames[2])
plot!(zzsk2["SK"][:, 1], zzsk2["SK"][:, 3], title=dimnames[2])
plot!(zzsk3["SK"][:, 1], zzsk3["SK"][:, 3], title=dimnames[2])
plot!(zzsk4["SK"][:, 1], zzsk4["SK"][:, 3], title=dimnames[2])

plot(zzsk1["SK"][:, 1], zzsk1["SK"][:, 4], title=dimnames[3])
plot!(zzsk2["SK"][:, 1], zzsk2["SK"][:, 4], title=dimnames[3])
plot!(zzsk3["SK"][:, 1], zzsk3["SK"][:, 4], title=dimnames[3])
plot!(zzsk4["SK"][:, 1], zzsk4["SK"][:, 4], title=dimnames[3])

plot(zzsk1["SK"][:, 2], zzsk1["SK"][:, 3], xlabel=dimnames[1], ylabel=dimnames[2])
plot!(zzsk2["SK"][:, 2], zzsk2["SK"][:, 3], xlabel=dimnames[1], ylabel=dimnames[2])
plot!(zzsk3["SK"][:, 2], zzsk3["SK"][:, 3], xlabel=dimnames[1], ylabel=dimnames[2])
plot!(zzsk4["SK"][:, 2], zzsk4["SK"][:, 3], xlabel=dimnames[1], ylabel=dimnames[2])

plot(zzsk1["SK"][:, 2], zzsk1["SK"][:, 4], xlabel=dimnames[1], ylabel=dimnames[3])
plot!(zzsk2["SK"][:, 2], zzsk2["SK"][:, 4], xlabel=dimnames[1], ylabel=dimnames[3])
plot!(zzsk3["SK"][:, 2], zzsk3["SK"][:, 4], xlabel=dimnames[1], ylabel=dimnames[3])
plot!(zzsk4["SK"][:, 2], zzsk4["SK"][:, 4], xlabel=dimnames[1], ylabel=dimnames[3])

plot(zzsk1["SK"][:, 3], zzsk1["SK"][:, 4], xlabel=dimnames[2], ylabel=dimnames[3])
plot!(zzsk2["SK"][:, 3], zzsk2["SK"][:, 4], xlabel=dimnames[2], ylabel=dimnames[3])
plot!(zzsk3["SK"][:, 3], zzsk3["SK"][:, 4], xlabel=dimnames[2], ylabel=dimnames[3])
plot!(zzsk4["SK"][:, 3], zzsk4["SK"][:, 4], xlabel=dimnames[2], ylabel=dimnames[3])

# Lets check how HCM behaves
Lε_tuned=0.03
L_tuned=2
hmc = runHMC(;epsilon=Lε_tuned/L_tuned,L=L_tuned,IT=1000,qs=start)
sum(hmc["accept"])
plot(autocor(hmc["SampleQ"][:, 1]))
plot(autocor(hmc["SampleQ"][:, 2]))
plot(autocor(hmc["SampleQ"][:, 3]))
plot(hmc["SampleQ"][:, 1], title=dimnames[1])
plot(hmc["SampleQ"][:, 2], title=dimnames[2])
plot(hmc["SampleQ"][:, 3], title=dimnames[3])

plot(hmc["SampleQ"][:, 1], hmc["SampleQ"][:, 2], xlabel=dimnames[1], ylabel=dimnames[2])
plot(hmc["SampleQ"][:, 1], hmc["SampleQ"][:, 3], xlabel=dimnames[1], ylabel=dimnames[3])
plot(hmc["SampleQ"][:, 2], hmc["SampleQ"][:, 3], xlabel=dimnames[2], ylabel=dimnames[3])

hmc1 = runHMC(;epsilon=Lε_tuned/L_tuned,L=L_tuned,IT=1000,qs=start1)
hmc2 = runHMC(;epsilon=Lε_tuned/L_tuned,L=L_tuned,IT=1000,qs=start2)
hmc3 = runHMC(;epsilon=Lε_tuned/L_tuned,L=L_tuned,IT=1000,qs=start3)
hmc4 = runHMC(;epsilon=Lε_tuned/L_tuned,L=L_tuned,IT=1000,qs=start4)


plot(hmc1["SampleQ"][:, 1], title=dimnames[1])
plot!(hmc2["SampleQ"][:, 1], title=dimnames[1])
plot!(hmc3["SampleQ"][:, 1], title=dimnames[1])
plot!(hmc4["SampleQ"][:, 1], title=dimnames[1])

plot(hmc1["SampleQ"][:, 2], title=dimnames[2])
plot!(hmc2["SampleQ"][:, 2], title=dimnames[2])
plot!(hmc3["SampleQ"][:, 2], title=dimnames[2])
plot!(hmc4["SampleQ"][:, 2], title=dimnames[2])

plot(hmc1["SampleQ"][:, 3], title=dimnames[3])
plot!(hmc2["SampleQ"][:, 3], title=dimnames[3])
plot!(hmc3["SampleQ"][:, 3], title=dimnames[3])
plot!(hmc4["SampleQ"][:, 3], title=dimnames[3])

plot(hmc1["SampleQ"][:, 1],hmc1["SampleQ"][:, 2], xlabel=dimnames[1], ylabel=dimnames[2])
plot!(hmc2["SampleQ"][:, 1],hmc2["SampleQ"][:,2], xlabel=dimnames[1], ylabel=dimnames[2])
plot!(hmc3["SampleQ"][:, 1],hmc3["SampleQ"][:, 2], xlabel=dimnames[1], ylabel=dimnames[2])
plot!(hmc4["SampleQ"][:, 1],hmc4["SampleQ"][:, 2], xlabel=dimnames[1], ylabel=dimnames[2])

plot(hmc1["SampleQ"][:, 1],hmc1["SampleQ"][:, 3], xlabel=dimnames[1], ylabel=dimnames[3])
plot!(hmc2["SampleQ"][:, 1],hmc2["SampleQ"][:, 3], xlabel=dimnames[1], ylabel=dimnames[3])
plot!(hmc3["SampleQ"][:, 1],hmc3["SampleQ"][:, 3], xlabel=dimnames[1], ylabel=dimnames[3])
plot!(hmc4["SampleQ"][:, 1],hmc4["SampleQ"][:, 3], xlabel=dimnames[1], ylabel=dimnames[3])

plot(hmc1["SampleQ"][:, 2],hmc1["SampleQ"][:, 3], xlabel=dimnames[2], ylabel=dimnames[3])
plot!(hmc2["SampleQ"][:, 2],hmc2["SampleQ"][:, 3], xlabel=dimnames[2], ylabel=dimnames[3])
plot!(hmc3["SampleQ"][:, 2],hmc3["SampleQ"][:, 3], xlabel=dimnames[2], ylabel=dimnames[3])
plot!(hmc4["SampleQ"][:, 2],hmc4["SampleQ"][:, 3], xlabel=dimnames[2], ylabel=dimnames[3])

# SUB SAMPLING LIKELIHOOD AND FUNCTIONS
# function Uj(x::Vector; j::Integer, t = t_vec, c = c_vec, z1= z1_vec)
#     # transform parameter
#     α = exp(x[1])
#     J = length(t)
#     lgμ_j = x[2] + x[3] *z1[j]
#     μ_j   = exp(lgμ_j)
#     if c[j]==1
#         mll_j = -log(α)+(1-α)*log(t[j])+α*log(μ_j)+((t[j]/μ_j)^α)
#     elseif c[j]==0
#         mll_j = ((t[j]/μ_j)^α)
#     end
#     return J*mll_j
# end

function Uj(x::Vector; j::Vector{Int}, t = t_vec, c = c_vec, z1= z1_vec)
    # transform parameter
    α = exp(x[1])
    J = length(t)
    # sub sample size
    sss = length(j)
    mll = 0
    for s in 1:sss
        lgμ_j = x[2] + x[3] *z1[j[s]]
        μ_j   = exp(lgμ_j)
        if c[j[s]]==1
            mll+= -log(α)+(1-α)*log(t[j[s]])+α*log(μ_j)+((t[j[s]]/μ_j)^α)
        elseif c[j[s]]==0
            mll+= ((t[j[s]]/μ_j)^α)
        end
    end
    return (J/sss)*mll
end



typeof([1,2,3])


Uj([-0.5, 9, -0.5], j=[1])
Uj([-0.5, 9, -0.5], j=[1,2])
Uj([-0.5, 9, -0.5], j=[4, 5, 3])
Uj([-0.5, 9, -0.5], j=[4])

# chose a sensible value for the CV x̃
x̃    = [-0.5, 9, -0.5]
∇tot = ForwardDiff.gradient(U, x̃)
# and compute the jth specifc gradient in x̃
∇j   = Array{Float64, 2}(undef, length(t_vec), Dim)
for jth in 1:length(t_vec)
    ∇j[jth, :] = ForwardDiff.gradient(x -> Uj(x; j=[jth]), x̃)
end
∇j

function rateswitchj(t; j, x0, v0, grad=∇tot, gradj=∇j,  Uxj=Uj)
    xt= x0+v0*t
    λ = max.((ForwardDiff.gradient(x -> Uxj(x; j=j), xt)+grad-
                vec(mapslices(mean, gradj[j,:], dims=1))).*v0, 0.0)
    return λ
end




# check if the rate works
plot(plot(t -> rateswitchj(t; j=[2, 3,4,5,6], x0=[-0.5, 9, -0.5], v0=[1,1,1])[1], 0,5),
    plot(t -> rateswitchj(t; j=[2, 3,4,5,6], x0=[-0.5, 9, -0.5], v0=[1,1,1])[2], 0,5),
    plot(t -> rateswitchj(t; j=[2, 3,4,5,6], x0=[-0.5, 9, -0.5], v0=[1,1,1])[3], 0,5), layout=(3,1))

# # now add the function to get an estimate of the maximum
# # parameters


idsamp = reshape(sample(1:J,30*10), 30, 10)

# now add the function to get an estimate of the maximum
# parameters
function getMestPOT(;ss , x, v, tmax, J=length(t_vec), D=Dim,
        sss=5) # this is how many points in each rate to select
    # ss=30*20
    # x=[-0.5, 9, -0.5]
    # v=[1,1,1]
    # tmax=0.1
    # J=length(t_vec)
    # D=Dim
    # get the index of the observations used to estimate the max
    idsamp = reshape(sample(1:J,ss*sss), ss, sss)
    # maximum rate per each dimension (this is unclean might include some nan)
    λ̄uncl     = zeros(ss, 3)
    for s in 1:ss
        for d in 1:D
            Optd_s= optimize(t -> -rateswitchj(t; j=idsamp[s,:], x0=x, v0=v)[d], 0, tmax, iterations=10)
            λ̄uncl[s, d]= -Optim.minimum(Optd_s)
        end
    end
    λ̄ = λ̄uncl[.!isnan.(vec(mapslices(sum, λ̄uncl, dims=2))),:]
    u = mapslices(x -> quantile(x, 0.95), λ̄, dims=1)
    X = Dict()
    Y = Dict()
    k   = zeros(D)
    n   = zeros(D)
    for d in 1:D
        X[d] = filter(x->x>=u[d],λ̄[:,d])
        Y[d] = X[d].-u[d]
        k[d] = length(Y[d])
        n[d] = length(λ̄[:, d])
    end
    # p(x>u)
    ζ̂   = k./n
    # return value
    m  = J
    # parameters
    Mest= zeros(D)
    for d in 1:D
        if all(Y[d].==0)
            Mest[d]= 0
        elseif (u[d]==0)
            Mest[d]= maximum(Y[d])*2
        else
            fitoutd = gpfit(Y[d])
            ϕ̂d  = fitoutd.θ̂[1]
            ξ̂d  = fitoutd.θ̂[2]
            Mest[d] = u[d]+(exp(ϕ̂d))/ξ̂d*((m*ζ̂[d])^ξ̂d -1)
        end
    end
    while any(isnan.(Mest))
        idsamp = sample(1:J,ss)
        # maximum rate per each dimension (this is unclean might include some nan)
        λ̄uncl     = zeros(ss, 3)
        for s in 1:ss
            for d in 1:D
                Optd_s= optimize(t -> -rateswitchj(t; j=idsamp[s], x0=x, v0=v)[d], 0, tmax, iterations=10)
                λ̄uncl[s, d]= -Optim.minimum(Optd_s)
            end
        end
        λ̄ = λ̄uncl[.!isnan.(vec(mapslices(sum, λ̄uncl, dims=2))),:]
        u = mapslices(x -> quantile(x, 0.95), λ̄, dims=1)
        X = Dict()
        Y = Dict()
        k   = zeros(D)
        n   = zeros(D)
        for d in 1:D
            X[d] = filter(x->x>=u[d],λ̄[:,d])
            Y[d] = X[d].-u[d]
            k[d] = length(Y[d])
            n[d] = length(λ̄[:, d])
        end
        # p(x>u)
        ζ̂   = k./n
        # return value
        m  = J
        # parameters
        Mest= zeros(D)
        for d in 1:D
            if all(Y[d].==0)
                Mest[d]= 0
            elseif u[d]==0
                Mest[d]= maximum(Y[d])*2
            else
                fitoutd = gpfit(Y[d])
                ϕ̂d  = fitoutd.θ̂[1]
                ξ̂d  = fitoutd.θ̂[2]
                Mest[d] = u[d]+(exp(ϕ̂d))/ξ̂d*((m*ζ̂[d])^ξ̂d -1)
            end
        end
    end
    return 2*Mest
end


# Random.seed!(1234)
# getMestPOTbis(ss=30*20, x=[-0.5, 9, -0.5], v=[+1,+1, +1], tmax=0.1, J=length(t_vec)) # chosen block size|

Random.seed!(1234)
getMestPOT(ss=30*20, x=[-0.5, 9, -0.5], v=[+1,+1, +1], tmax=0.1, J=length(t_vec)) # chosen block size|

Random.seed!(1234)

NS=1000      # number of skeleton points
x0_0=[-0.5, 9, -0.5]         # initial location
v0_0=missing    # initial velocity
tmax=0.025            # tmax for tuning
B=false        # budget available, if false run until NS
roundopt=true # whether to use our rounded version of the optimization
ε₁ = 1.00e-20   # boundary error : how close to get before turning
ssM = 20*40
NOBS=length(t_vec)


# check if the stopping criterion has been defined and how
if (NS==false && B==false)
    error("Choose one stopping criterion: e.g. either: NS=10000; B=500000")
elseif (NS != false && B !=false)
    error("No multiple stopping criteria: choose either NS=... or  B=... \n
    set the other to false")
end

# get an estimate of the size of the sketon (SS) :
# if that is not given compute it from B and T and
# add at the end some checks to make sure
# we append more lines if there is the need
(NS != false) ? SS = NS : SS = B

# check the dimension of the problem
if (size(x0_0, 2)==1)
    D    = size(x0_0, 1)
else
    error("the intial value x0_0 must be a vertical vector of dimensions Dx1, e.g. [x_01, x_02, ..., x_0D]")
end
# set up the array for the skeleton :
# note: these are horizonatal because it should be more efficient for Julia
x0set = Array{Float64, 2}(undef, D, SS)
v0set = Array{Float64, 2}(undef, D, SS)
t0set = Array{Float64, 2}(undef, 1, SS)
# Array of the number of gradient evaluations with three rows
# 1st: check horizon; 2nd: optimization; 3rd tpp
GradEvals = zeros(Float64, 3, SS)
ErrorOpt  = zeros(Float64, 1, SS)

# setting up they state at the begining of the location
x0set[:, 1] = x0_0

# set velocity to one if missing
if (ismissing(v0_0))
    v0set[:, 1] = ones(D)
else
    v0set[:, 1] = v0_0
end
t0set[1, 1] = 0.0


x0i = x0set[:, 1]
v0i = v0set[:, 1]
ts  = 0
tp  = 0
horizon=tmax
k   = 2

# stopping criterion
NS != false ? keepgoing = (k <= NS) : keepgoing = sum(GradEvals)<=B

# TOADD ↓
# Check if the initial values are in a section were the likelihood is not defined.
# TOADD ↑
while keepgoing         # Run the loop until the chosen condition is met
# check for bounds at the horizon
    while isnan(globalrate(horizon; x0=x0i, v0=v0i)) && horizon>ε₁
        horizon=horizon/2
        # count evaluations [added to the horizon change section]
        GradEvals[1, k] = GradEvals[1, k] + 1
    end
    # if approached the horizon (at the chosen ε₁) switch back
    if horizon <=ε₁
        error("Possible border ahead")
    else # continue with the current horizon
        M   = getMestPOT(ss=ssM, x=x0i, v=v0i, tmax=horizon, J=NOBS)
        # count evaluations [added to the optimization evaluation section]
        GradEvals[2, k] = GradEvals[2, k] + D  #*10 or iterations per grad evals
        if sum(M) == 0   # move dererministically if the rate of switching is 0 in all dims
            ts = ts+horizon
            x0i= x0i+horizon*v0i
        else        # propose a time for the switch in the dimensions
            tp_d   = rand.(Exponential.(1 ./M))
            i0     = argmin(tp_d)
            tp     = minimum(tp_d)

            if tp >= horizon    # move deterministically if reached the horizon
                ts = ts+horizon
                x0i= x0i+horizon*v0i
            else                # evaluate proposal
                accept = false  # start evaluations
                while (tp < horizon) && (accept == false)
                    j_0  = sample(1:NOBS, 5)
                    ar = (rateswitchj(tp; j=j_0, x0=x0i, v0=v0i)[i0])/M[i0]
                    # count evaluations [added to the thinned section]
                    GradEvals[3, k] = GradEvals[3, k] + 1
                    if ar > 1   # if optimization was wrong
                        horizon = tp
                        M   = getMestPOT(ss=ssM, x=x0i, v=v0i, tmax=horizon, J=NOBS)
                        tp_d   = rand.(Exponential.(1 ./M))
                        i0     = argmin(tp_d)
                        tp     = minimum(tp_d)
                        # count evaluations [added to the optimization evaluation section]
                        GradEvals[2, k] = GradEvals[2, k] + D  #*10 or iterations per grad evals
                        ErrorOpt[1, k] = ErrorOpt[1, k] + 1
                    else        # evaluate acceptance
                        if (rand(Bernoulli(ar)))
                            # update location and switch velocity
                            x0i    = x0i + tp * v0i
                            v0i[i0] = -v0i[i0]
                            # save the skeleton point
                            v0set[:, k] = v0i
                            x0set[:, k] = x0i
                            t0set[1, k] = t0set[1, (k-1)]+ts+tp
                            # reset time from skeleton point, horizon,
                            # flag acceptance and increase counter
                            ts = 0.0
                            tp = 0.0
                            horizon = tmax
                            k = k+1
                            accept = true
                            print(string(k ,"\n"))
                        else   # upon rejection increas stochastic time
                            tp_d   = rand.(Exponential.(1 ./M))
                            i0     = argmin(tp_d)
                            tp     = tp+minimum(tp_d)                        end
                    end

                end
                if tp >= horizon && (accept == false) # if exited while loop because horizon reached
                    ts = ts+horizon
                    x0i= x0i+horizon*v0i
                end
            end
        end
    end
    NS != false ? keepgoing = (k <= NS) : keepgoing = sum(GradEvals)<=B

end
if B != false
    x0set=x0set[:,1:(k-1)]
    v0set=v0set[:,1:(k-1)]
    t0set=t0set[:,1:(k-1)]
    GradEvals=GradEvals[:,1:(k-1)]
    ErrorOpt=ErrorOpt[:,1:(k-1)]
end
outsk=hcat(transpose(t0set),transpose(x0set),transpose(v0set),
    transpose(GradEvals),transpose(ErrorOpt))
output=Dict([("SkeletonLocation", x0set), ("SkeletonVelocity", v0set), ("SkeletonTime", t0set)
  , ("GradientEvaluations", GradEvals), ("Errorsoptimization", ErrorOpt), ("SK",outsk)])
return(output)

p1= plot(outsk[:, 1], outsk[:, 2], title=dimnames[1])
plot!(zzsk["SK"][:, 1], zzsk["SK"][:, 2])
p2=plot(outsk[:, 1], outsk[:, 3], title=dimnames[2])
plot!(zzsk["SK"][:, 1], zzsk["SK"][:, 3])
p3=plot(outsk[:, 1], outsk[:, 4], title=dimnames[3])
plot!(zzsk["SK"][:, 1], zzsk["SK"][:, 4])
savefig(p1, string(SAVEwd, "sk1.pdf"))
savefig(p2, string(SAVEwd, "sk2.pdf"))
savefig(p3, string(SAVEwd, "sk3.pdf"))
outss = zzsample(;N=1000, sk=output)
outnss = zzsample(;N=1000, sk=zzsk)

p4=density(outss[100:1000,1], label="with ss", legend=true, linewidth=2)
density!(outnss[100:1000,1], label="without ss",  linewidth=2)


p5=density(outss[100:1000,2], label="with ss", legend=true, linewidth=2)
density!(outnss[100:1000,2], label="without ss", linewidth=2)


p6=density(outss[100:1000,3], label="with ss", legend=true, linewidth=2)
density!(outnss[100:1000,3], label="without ss",  linewidth=2)

savefig(p4, string(SAVEwd, "d1.pdf"))
savefig(p5, string(SAVEwd, "d2.pdf"))
savefig(p6, string(SAVEwd, "d3.pdf"))
