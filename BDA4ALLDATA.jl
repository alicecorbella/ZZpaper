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
    mll = 0
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

Dim  =4
start = [-0.5, +9, -0.4, -0.4]
start2= [-0.15, 13, -0.1, -5]
zzsk  =  zz(; NS=10, x0_0=start, tmax=0.02)
zzsk2 =  zz(; NS=10, x0_0=start2, tmax=0.02)

p1=plot(zzsk["SK"][:, 1], zzsk["SK"][:, 2], title=dimnames[1])
# plot!(zzsk2["SK"][:, 1], zzsk2["SK"][:, 2], title=dimnames[1])
p2=plot(zzsk["SK"][:, 1], zzsk["SK"][:, 3], title=dimnames[2])
# plot!(zzsk2["SK"][:, 1], zzsk2["SK"][:, 3], title=dimnames[2])
p3=plot(zzsk["SK"][:, 1], zzsk["SK"][:, 4], title=dimnames[3])
# plot!(zzsk2["SK"][:, 1], zzsk2["SK"][:, 4], title=dimnames[3])
p4=plot(zzsk["SK"][:, 1], zzsk["SK"][:, 5], title=dimnames[4])
# plot!(zzsk2["SK"][:, 1], zzsk2["SK"][:, 5], title=dimnames[4])
p5 = plot(p1,p2,p3,p4, layout=(2,2), size=(600,600))

# savefig(p5, string(SAVEwd, "4parsk.pdf"))

# SUBSAMPLING

function Uj(x::Vector;j::Integer,  t = t_vec, c = c_vec, z1= z1_vec, z2=z2_vec)
    # transform parameter
    α = exp(x[1])
    J = length(t)
    lgμ_j = x[2] + x[3] *z1[j]+ x[4] *z2[j]
    μ_j   = exp(lgμ_j)
    if c[j]==1
        mll_j = -log(α)+(1-α)*log(t[j])+α*log(μ_j)+((t[j]/μ_j)^α)
    elseif c[j]==0
        mll_j = ((t[j]/μ_j)^α)
    end
    return J*mll_j
end


# compute the gradient at "mode" for CV
x̃    = [-0.4, 9.5, -0.3, -2]
∇tot = ForwardDiff.gradient(U, x̃)
∇j   = Array{Float64, 2}(undef, length(t_vec), Dim)
for jth in 1:length(t_vec)
    ∇j[jth, :] = ForwardDiff.gradient(x -> Uj(x; j=jth), x̃)
end
# done

zzskss=zz_w_ss(; NS=1000, x0_0=start, tmax=0.025,ssM =10000, NOBS=length(t_vec))


p1=plot(zzskss["SK"][:, 1], zzskss["SK"][:, 2], title=dimnames[1])
plot!(zzsk["SK"][:, 1], zzsk["SK"][:, 2])
p2=plot(zzskss["SK"][:, 1], zzskss["SK"][:, 3], title=dimnames[2])
plot!(zzsk["SK"][:, 1], zzsk["SK"][:, 3])
p3=plot(zzskss["SK"][:, 1], zzskss["SK"][:, 4], title=dimnames[3])
plot!(zzsk["SK"][:, 1], zzsk["SK"][:, 4])
p4=plot(zzskss["SK"][:, 1], zzskss["SK"][:, 5], title=dimnames[4])
plot!(zzsk["SK"][:, 1], zzsk["SK"][:, 5])
p5 = plot(p1,p2,p3,p4, layout=(2,2), size=(600,600))

savefig(p5, string(SAVEwd, "ALLSScomp.pdf"))
