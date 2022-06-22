# Direcotries
CODEwd = "/Users/alice/projects/ZZpaper/"
SAVEwd = "/Users/alice/OneDrive - University of Warwick/Manuscripts/04ADZZ/figures/sec6.2/"

# Loading the packages
include(string(CODEwd, "env.jl"))
# Loading the functions for general zig zag
include(string(CODEwd, "functions.jl"))
include(string(CODEwd, "plotting.jl"))

# Include supplementary packages for data IO
using CSV
using DataFrames
using StatsBase
using Serialization
using Dates
using Extremes

Y = open(deserialize, string(CODEwd,"data/BDA_preprocessed.bin"))

# check why this missing has escaped the filters
Y=Y[.!ismissing.(Y.TIMEDIAGNOSITOFINALEVENT),:]

(print(describe(Y)))

# model with age and STAGE if spreaded
dimnames=["lgα", "β₀", "β₁", "β₂"]
t_vec = Dates.value.(Y.TIMEDIAGNOSITOFINALEVENT).+0.01
c_vec = zeros(length(t_vec))
c_vec[findall(Y.NEWVITALSTATUS.=="D")] .= 1
z1_vec = (Y.AGE.-mean(Y.AGE))./std(Y.AGE)
z2_vec = zeros(length(t_vec))
z2_vec[findall(Y.STAGENUMBER.>=3)] .= 1

function U(x::Vector; t = t_vec, c = c_vec, z1= z1_vec, z2=z2_vec)
    # transform parameter
    α = exp(x[1])
    J = length(t)
    mll = -x[1]
    for j in 1:J
        lgμ_j = x[2] + x[3] *z1[j]+ x[4] *z2[j]
        μ_j   = exp(lgμ_j)
        if c[j]==1
            mll+= -log(α)+(1-α)*log(t[j])+α*log(μ_j)+((t[j]/μ_j)^α)
        elseif c[j]==0
            mll+= ((t[j]/μ_j)^α)
        end
    end
    return mll
end


U([-0.5, +9, -0.4, -0.4])
U([-0.429, 9.421, -0.354, -2.073])
Dim  = 4
start = [-0.429, 9.421, -0.354, -2.073]

opti = optimize(U, start)
x̃=opti.minimizer
start=x̃

# run function for 5 iters just for first compilation
zz(; NS=5, x0_0=start, tmax=0.01)

now()
zzsk_smalltm=zz(; NS=5000, x0_0=start, tmax=0.005)
now()

now()
zzsk=zz(; NS=5000, x0_0=start, tmax=0.01)
now()
zzsk
JLD.save(string(SAVEwd, "/zzALLOBS.jld"),Dict("zz1"=>zzsk))

# zzlong=JLD.load(string(SAVEwd, "/zzALLOBS.jld"))
# zzsk=zzlong["zz1"]

sum(zzsk["GradientEvaluations"])

p1=plot(zzsk["SK"][:, 1], zzsk["SK"][:, 2], title=dimnames[1])
p2=plot(zzsk["SK"][:, 1], zzsk["SK"][:, 3], title=dimnames[2])
p3=plot(zzsk["SK"][:, 1], zzsk["SK"][:, 4], title=dimnames[3])
p4=plot(zzsk["SK"][:, 1], zzsk["SK"][:, 5], title=dimnames[4])
p5 = plot(p1,p2,p3,p4, layout=(2,2), size=(600,600))

# savefig(p5, string(SAVEwd, "4parsk.pdf"))

# obtain an estimate of the mode to
smp = zzsample(N=5000, sk=zzsk)


# SUBSAMPLING


function Uj(x::Vector;j::Vector{Int},  t = t_vec, c = c_vec, z1= z1_vec, z2=z2_vec)
    # transform parameter
    α = exp(x[1])
    J = length(t)
    sss = length(j)  #sub-sample size
    mll = 0
    for s in 1:sss
        lgμ_j = x[2] + x[3] *z1[j[s]]+ x[4] *z2[j[s]]
        μ_j   = exp(lgμ_j)
        if c[j[s]]==1
            mll+= -log(α)+(1-α)*log(t[j[s]])+α*log(μ_j)+((t[j[s]]/μ_j)^α)
        elseif c[j[s]]==0
            mll+= ((t[j[s]]/μ_j)^α)
        end
    end
    return (J/sss)*mll -x[1]
end


# compute the gradient at  x̃  for CV
∇tot = ForwardDiff.gradient(U, x̃)
∇j   = Array{Float64, 2}(undef, length(t_vec), Dim)
for jth in 1:length(t_vec)
    ∇j[jth, :] = ForwardDiff.gradient(x -> Uj(x; j=[jth]), x̃)
end
# done

# run function for 5 iters just for first compilation
zz_w_ss(; NS=5, x0_0=start, tmax=0.005,ssM =1000,
        NOBS=length(t_vec), ssS=50)

now()
print(now())
zzskss_50=zz_w_ss(; NS=5000, x0_0=start, tmax=0.005,ssM =1000,
        NOBS=length(t_vec), ssS=50)
now()
print(now())


now()
print(now())
zzskss_20=zz_w_ss(; NS=5000, x0_0=start, tmax=0.005,ssM =1000,
        NOBS=length(t_vec), ssS=20)
now()
print(now())


sum(zzskss_50["GradientEvaluations"])


px_dict=Dict()
for d in 1:Dim
    px=plot(zzsk["SK"][:, 1], zzsk["SK"][:, d+1], label="All obs",
        title=dimnames[d],legend=true,linewidth=0.9, color=:blue)
    plot!(zzskss_50["SK"][:, 1], zzskss_50["SK"][:, d+1],
        label="50 ss", linewidth=0.9, color=:black)
        px_dict[d]=px
end

pxsk=plot(px_dict[1], px_dict[2], px_dict[3], px_dict[4], layout=(2,2), size=(1000,800))

savefig(pxsk, string(SAVEwd, "/f1.pdf"))

zzsummaries(dms=1, sk=zzskss_50, B=100)["EffectiveSampleSize"]

ess1_b=Vector(undef, 10)
for b in 1:10
    ess1_b[b]=zzsummaries(dms=1, sk=zzsk, B=[2, 5, 10,50, 100,200,300,400, 500, 1000][b])["EffectiveSampleSize"]
end
ess1_b

ess1_b_ss=Vector(undef, 10)
for b in 1:10
    ess1_b_ss[b]=zzsummaries(dms=1, sk=zzskss_50, B=[2, 3, 10, 50,100,200,300,400, 500, 1000][b])["EffectiveSampleSize"]
end
ess1_b_ss
plot(ess1_b)

plot!(ess1_b_ss)


zzsummaries(dms=1, sk=zzskss_20, B=100)["EffectiveSampleSize"]
