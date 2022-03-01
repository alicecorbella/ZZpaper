# Direcotries
CODEwd = "/Users/alice/projects/ZZpaper/"
SAVEwd = "/Users/alice/OneDrive - University of Warwick/Manuscripts/04ADZZ/figures/prova/"

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

zzlong=JLD.load(string(SAVEwd, "/zzALLOBS.jld"))
zzsk=zzlong["zz1"]
zzsk2=zzlong["zz2"]
# obtain an estimate of the mode to
smp = zzsample(N=5000, sk=zzsk)
smp2 = zzsample(N=5000, sk=zzsk2)

med1=quantile(vcat(smp[3000:5000, 1],smp[3000:5000, 1]), 0.5)
med2=quantile(vcat(smp[3000:5000, 2],smp[3000:5000, 2]), 0.5)
med3=quantile(vcat(smp[3000:5000, 3],smp[3000:5000, 3]), 0.5)
med4=quantile(vcat(smp[3000:5000, 4],smp[3000:5000, 4]), 0.5)
tmax_tuned=0.0001


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
    return (J/sss)*mll
end


# compute the gradient at "median" for CV
x̃    = [med1, med2, med3, med4]
∇tot = ForwardDiff.gradient(U, x̃)
∇j   = Array{Float64, 2}(undef, length(t_vec), Dim)
for jth in 1:length(t_vec)
    ∇j[jth, :] = ForwardDiff.gradient(x -> Uj(x; j=[jth]), x̃)
end
# done

zzskss_3=zz_w_ss(; NS=5000, x0_0=start, tmax=0.002,ssM =1000,
        NOBS=length(t_vec), ssS=3)
zzskss_5=zz_w_ss(; NS=5000, x0_0=start, tmax=0.002,ssM =1000,
        NOBS=length(t_vec), ssS=5)
zzskss_10=zz_w_ss(; NS=5000, x0_0=start, tmax=0.002,ssM =1000,
        NOBS=length(t_vec), ssS=10)
zzskss_50=zz_w_ss(; NS=5000, x0_0=start, tmax=0.002,ssM =1000,
        NOBS=length(t_vec), ssS=50)
zzskss_100=zz_w_ss(; NS=5000, x0_0=start, tmax=0.002,ssM =1000,
        NOBS=length(t_vec), ssS=100)
zzskss_200=zz_w_ss(; NS=5000, x0_0=start, tmax=0.002,ssM =1000,
        NOBS=length(t_vec), ssS=200)
zzskss_1000=zz_w_ss(; NS=5000, x0_0=start, tmax=0.002,ssM =1000,
        NOBS=length(t_vec), ssS=1000)

# JLD.save(string(SAVEwd, "/zzSSs.jld"),
#     Dict("zzskss_5"=>zzskss_5, "zzskss_10"=>zzskss_10,
#     "zzskss_50"=>zzskss_50, "zzskss_100"=>zzskss_100,
#     "zzskss_200"=>zzskss_200, "zzskss_1000"=>zzskss_1000))


p1=plot(zzskss_5["SK"][:, 1], zzskss_5["SK"][:, 2], title=dimnames[1],
    label="5 ss", legend=true)
plot!(zzskss_10["SK"][:, 1], zzskss_10["SK"][:, 2], label="10 ss")
plot!(zzskss_50["SK"][:, 1], zzskss_50["SK"][:, 2], label="50 ss")
plot!(zzskss_100["SK"][:, 1], zzskss_100["SK"][:, 2], label="100 ss")
plot!(zzskss_200["SK"][:, 1], zzskss_200["SK"][:, 2], label="200 ss")
plot!(zzskss_1000["SK"][:, 1], zzskss_1000["SK"][:, 2], label="1000 ss")
savefig(p1, string(SAVEwd, "SS1_onlyss.pdf"))
plot!(zzsk["SK"][1:2000, 1], zzsk["SK"][1:2000, 2], label="All obs")



p2=plot(zzskss_5["SK"][:, 1], zzskss_5["SK"][:, 3], title=dimnames[2],
    label="5 ss", legend=true)
plot!(zzskss_10["SK"][:, 1], zzskss_10["SK"][:, 3], label="10 ss")
plot!(zzskss_50["SK"][:, 1], zzskss_50["SK"][:, 3], label="50 ss")
plot!(zzskss_100["SK"][:, 1], zzskss_100["SK"][:, 3], label="100 ss")
plot!(zzskss_200["SK"][:, 1], zzskss_200["SK"][:, 3], label="200 ss")
plot!(zzskss_1000["SK"][:, 1], zzskss_1000["SK"][:, 3], label="1000 ss")
savefig(p2, string(SAVEwd, "SS2_onlyss.pdf"))
plot!(zzsk["SK"][1:2000, 1], zzsk["SK"][1:2000, 3], label="All obs")




p3=plot(zzskss_5["SK"][:, 1], zzskss_5["SK"][:, 4], title=dimnames[3],
    label="5 ss", legend=true)
plot!(zzskss_10["SK"][:, 1], zzskss_10["SK"][:, 4], label="10 ss")
plot!(zzskss_50["SK"][:, 1], zzskss_50["SK"][:, 4], label="50 ss")
plot!(zzskss_100["SK"][:, 1], zzskss_100["SK"][:, 4], label="100 ss")
plot!(zzskss_200["SK"][:, 1], zzskss_200["SK"][:, 4], label="200 ss")
plot!(zzskss_1000["SK"][:, 1], zzskss_1000["SK"][:, 4], label="1000 ss")
savefig(p3, string(SAVEwd, "SS3_onlyss.pdf"))
plot!(zzsk["SK"][1:2000, 1], zzsk["SK"][1:2000, 4], label="All obs")




p4=plot(zzskss_5["SK"][:, 1], zzskss_5["SK"][:, 5], title=dimnames[4],
    label="5 ss", legend=true)
plot!(zzskss_10["SK"][:, 1], zzskss_10["SK"][:, 5], label="10 ss")
plot!(zzskss_50["SK"][:, 1], zzskss_50["SK"][:, 5], label="50 ss")
plot!(zzskss_100["SK"][:, 1], zzskss_100["SK"][:, 5], label="100 ss")
plot!(zzskss_200["SK"][:, 1], zzskss_200["SK"][:, 5], label="200 ss")
plot!(zzskss_1000["SK"][:, 1], zzskss_1000["SK"][:, 5], label="1000 ss")
savefig(p4, string(SAVEwd, "SS4_onlyss.pdf"))
plot!(zzsk["SK"][1:2000, 1], zzsk["SK"][1:2000, 5], label="All obs")

savefig(p1, string(SAVEwd, "SS1.pdf"))
savefig(p2, string(SAVEwd, "SS2.pdf"))
savefig(p3, string(SAVEwd, "SS3.pdf"))
savefig(p4, string(SAVEwd, "SS4.pdf"))



zzskss_frommode_10=zz_w_ss(; NS=5000, x0_0=x̃, tmax=0.002,ssM =1000,
        NOBS=length(t_vec), ssS=10)
zzskss_frommode_20=zz_w_ss(; NS=5000, x0_0=x̃, tmax=0.002,ssM =1000,
        NOBS=length(t_vec), ssS=20)
zzskss_frommode_50=zz_w_ss(; NS=5000, x0_0=x̃, tmax=0.002,ssM =1000,
        NOBS=length(t_vec), ssS=50)
px1=plot(zzskss_frommode_50["SK"][:, 1], zzskss_frommode_50["SK"][:, 2],
    title=dimnames[1],label="50 ss", legend=true, linewidth=0.5)
plot!(zzsk["SK"][:, 1], zzsk["SK"][:, 2], label="All obs", linewidth=0.5)
plot!(zzskss_frommode_20["SK"][:, 1], zzskss_frommode_20["SK"][:, 2],
    label="20 ss", linewidth=0.5)

px2=plot(zzskss_frommode_50["SK"][:, 1], zzskss_frommode_50["SK"][:, 3],
    title=dimnames[2],label="50 ss", legend=true, linewidth=0.5)
plot!(zzsk["SK"][:, 1], zzsk["SK"][:, 3], label="All obs", linewidth=0.5)

px3=plot(zzskss_frommode_50["SK"][:, 1], zzskss_frommode_50["SK"][:, 4],
    title=dimnames[3],label="50 ss", legend=true, linewidth=0.5)
plot!(zzsk["SK"][:, 1], zzsk["SK"][:, 4], label="All obs", linewidth=0.5)

px4=plot(zzskss_frommode_50["SK"][:, 1], zzskss_frommode_50["SK"][:, 5],
    title=dimnames[4],label="50 ss", legend=true, linewidth=0.5)
plot!(zzsk["SK"][:, 1], zzsk["SK"][:, 5], label="All obs", linewidth=0.5)

savefig(px1, string(SAVEwd, "SS1_fm.pdf"))
savefig(px2, string(SAVEwd, "SS2_fm.pdf"))
savefig(px3, string(SAVEwd, "SS3_fm.pdf"))
savefig(px4, string(SAVEwd, "SS4_fm.pdf"))


smp = zzsample(N=5000, sk=zzsk)
smp_50ss = zzsample(N=5000, sk=zzskss_frommode_50)
smp_20ss = zzsample(N=5000, sk=zzskss_frommode_20)

density(smp[3000:5000, 1], color=:orange, linewidth=2)
density!(smp_50ss[3000:5000, 1], color=:blue, linewidth=2)

density(smp[2000:5000, 2], color=:orange, linewidth=2)
density!(smp_50ss[2000:5000, 2], color=:blue, linewidth=2)
density!(smp_20ss[2000:5000, 2], color=:red, linewidth=2)
density(smp[2000:5000, 3], color=:orange, linewidth=2)
density!(smp_50ss[2000:5000, 3], color=:blue, linewidth=2)
density(smp[2000:3500, 4], color=:orange, linewidth=2)
density!(smp_50ss[2000:5000, 4], color=:blue, linewidth=2)
density!(smp[3501:5000, 4], color=:red, linewidth=2)





ESS_wholedata = zeros(4)
ESS_withSS  = zeros(4)
for d in 1:4
    ESS_wholedata[d]=ESSbm(smpl=smp[2000:5000, d], nbatches=25)
    ESS_withSS[d]=ESSbm(smpl=smp_50ss[2000:5000, d], nbatches=25)
end
ESS_wholedata

ESS_withSS

print(round.(hcat(ESS_withSS, ESS_wholedata), digits=2))

# # reduced skeleton without burnin
# sk_red=Dict([("SkeletonLocation", zzsk["SkeletonLocation"][:, 2000:5000]),
#         ("SkeletonVelocity", zzsk["SkeletonVelocity"][:, 2000:5000]),
#         ("SkeletonTime", zzsk["SkeletonTime"][:, 2000:5000])])
#
# sk_ss_red=Dict([("SkeletonLocation", zzskss_frommode_50["SkeletonLocation"][:, 2000:5000]),
#         ("SkeletonVelocity", zzskss_frommode_50["SkeletonVelocity"][:, 2000:5000]),
#         ("SkeletonTime", zzskss_frommode_50["SkeletonTime"][:, 2000:5000])])
#







start = [-0.4, +9, -0.4, -0.4]
start2= [-0.45, 10, -0.5, -3]

start3= [-1, 9, -1, -1]



# parameters HMC
Lε_tuned=0.001
L_tuned=2

hmc1 = runHMC(;epsilon=Lε_tuned/L_tuned,L=L_tuned,IT=5000,
    qs=start)

hmc2 = runHMC(;epsilon=Lε_tuned/L_tuned,L=L_tuned,IT=100,
    qs=start3)


plot(hmc2["SampleQ"][:, 1], title=dimnames[1])
plot(hmc2["SampleQ"][:, 2], title=dimnames[2])
plot(hmc2["SampleQ"][:, 1], title=dimnames[3])
plot(hmc2["SampleQ"][:, 1], title=dimnames[4])


plot(autocor(hmc["SampleQ"][:, 1]))
plot!(autocor(hmc["SampleQ"][:, 2]))
plot!(autocor(hmc["SampleQ"][:, 3]))
plot!(autocor(hmc["SampleQ"][:, 4]))
