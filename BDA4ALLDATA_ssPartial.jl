# Direcotries
CODEwd = "/Users/alice/projects/ZZpaper/"
SAVEwd = "/Users/alice/OneDrive - University of Warwick/Manuscripts/04ADZZ/figures/sec6.2bis/"

# Loading the packages
include(string(CODEwd, "env.jl"))
# Loading the functions for general zig zag
include(string(CODEwd, "functions.jl"))
include(string(CODEwd, "plotting.jl"))

gr(size = (pt(84), pt(84)), labelfontsize=8, legend = false)

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

Dim  =4
start = [-0.4, +9, -0.4, -0.4]
start2= [-0.45, 10, -0.5, -3]
# zzsk  =  zz(; NS=5000, x0_0=start, tmax=0.01)
# zzsk2 =  zz(; NS=5000, x0_0=start2, tmax=0.01)
# JLD.save(string(SAVEwd, "/zzALLOBS.jld"),
#     Dict("zz1"=>zzsk, "zz2"=>zzsk2))
zzlong=JLD.load(("/Users/alice/OneDrive - University of Warwick/Manuscripts/04ADZZ/figures/sec6.2/zzALLOBS.jld"))
zzsk=zzlong["zz1"]
zzsk2=zzlong["zz2"]

p1=plot(zzsk["SK"][:, 1], zzsk["SK"][:, 2], ylabel=dimnames[1], color=:gray, linewidth=0.5)
plot!(zzsk2["SK"][:, 1], zzsk2["SK"][:, 2], color=:gray, linewidth=0.5)
p2=plot(zzsk["SK"][:, 1], zzsk["SK"][:, 3], ylabel=dimnames[2], color=:gray, linewidth=0.5)
plot!(zzsk2["SK"][:, 1], zzsk2["SK"][:, 3], color=:gray, linewidth=0.5)
p3=plot(zzsk["SK"][:, 1], zzsk["SK"][:, 4], ylabel=dimnames[3], color=:gray, linewidth=0.5)
plot!(zzsk2["SK"][:, 1], zzsk2["SK"][:, 4], color=:gray, linewidth=0.5)
p4=plot(zzsk["SK"][:, 1], zzsk["SK"][:, 5], ylabel=dimnames[4], color=:gray, linewidth=0.5)
plot!(zzsk2["SK"][:, 1], zzsk2["SK"][:, 5], color=:gray, linewidth=0.5)
p5 = plot(p1,p2,p3,p4, layout=(2,2), size=(pt(84*2),pt(84*2)))

# savefig(p5, string(SAVEwd, "4parsk.pdf"))

# obtain an estimate of the mode to
smp = zzsample(N=5000, sk=zzsk)
smp2 = zzsample(N=5000, sk=zzsk2)

med1=quantile(vcat(smp[3000:5000, 1],smp[3000:5000, 1]), 0.5)
p1=density(vcat(smp[3000:5000, 1], smp2[3000:5000, 1]), linewidth=2,
    title=string(dimnames[1],", Me=", round(med1, digits=3)))
vline!([med1], color=:black, linewidth=2)

med2=quantile(vcat(smp[3000:5000, 2],smp[3000:5000, 2]), 0.5)
p2=density(vcat(smp[3000:5000, 2], smp2[3000:5000, 2]), linewidth=2,
    title=string(dimnames[2],", Me=", round(med2, digits=3)))
vline!([med2], color=:black, linewidth=2)

med3=quantile(vcat(smp[3000:5000, 3],smp[3000:5000, 3]), 0.5)
p3=density(vcat(smp[3000:5000, 3], smp2[3000:5000, 3]), linewidth=2,
    title=string(dimnames[3],", Me=", round(med3, digits=3)))
vline!([med3], color=:black, linewidth=2)

med4=quantile(vcat(smp[3000:5000, 4],smp[3000:5000, 4]), 0.5)
p4=density(vcat(smp[3000:5000, 4], smp2[3000:5000, 4]), linewidth=2,
    title=string(dimnames[1],", Me=", round(med4, digits=3)))
vline!([med4], color=:black, linewidth=2)

p5 = plot(p1,p2,p3,p4, layout=(2,2), size=(pt(84*2),pt(84*2)))
# savefig(p5, string(SAVEwd, "4parsmp.pdf"))



# lets try and find the optimal tmax from the mode
tmplot(;R = 4, try_tmax=[0.0005, 0.001, 0.002, 0.005],
    ns=500, start=[med1, med2, med3, med4])
savefig(f2["NHPP"], string(SAVEwd, "f2a.pdf"))
savefig(f2["opt"], string(SAVEwd, "f2b.pdf"))
savefig(f2["tot"], string(SAVEwd, "f2c.pdf"))
f2tot= plot(f2["NHPP"],f2["opt"],f2["tot"], layout=(1,3),
    size=(pt(252), pt(60)))
savefig(f2tot, string(SAVEwd, "f2.pdf"))
tmax_tuned=0.002

zzsk_short  =  zz(; NS=50, x0_0=[med1, med2, med3, med4], tmax=tmax_tuned)

dplot=Dict()
for d in 1:4
    dplot[d] = plot(zzsk_short["SK"][:, 1], zzsk_short["SK"][:, d+1], title=dimnames[d])
end

plot(dplot[1],dplot[2],dplot[3],dplot[4], layout=(2,2), size=(pt(84*2),pt(84*2)))

# SUBSAMPLING
function Uj(x::Vector;j::Vector{Int},  t = t_vec, c = c_vec, z1= z1_vec, z2=z2_vec)
    # transform parameter
    α = exp(x[1])
    J = length(t)
    sss = length(j)  #sub-sample size
    mll = -x[1]
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
opti = optimize(U, [med1, med2, med3, med4])
x̃    = opti.minimizer
∇tot = ForwardDiff.gradient(U, x̃)
∇j   = Array{Float64, 2}(undef, length(t_vec), Dim)
for jth in 1:length(t_vec)
    ∇j[jth, :] = ForwardDiff.gradient(x -> Uj(x; j=[jth]), x̃)
end
∇j

# RATES AND FUNCTIONS FOR THE SUBSAPLING
function rateswitchj(t; j, x0, v0, grad=∇tot, gradj=∇j,  Uxj=Uj)
    xt= x0+v0*t
    λ = max.((ForwardDiff.gradient(x -> Uxj(x; j=j), xt)+grad-
            vec(mapslices(mean, gradj[j,:], dims=1))).*v0, 0.0)
    return λ
end


function rateswitchj_d1(t; j, x0, v0, ∂=∇tot[1], ∂j=∇j[:,1],  Uxj=Uj)
    xt= x0+v0*t
    λ1 = max((ForwardDiff.derivative(x -> Uxj([x, xt[2], xt[3], xt[4]]; j=j), xt[1])+
            ∂-mean(∂j[j]))*v0[1], 0.0)
    return λ1
end


function rateswitchj_d2(t; j, x0, v0, ∂=∇tot[2], ∂j=∇j[:,2],  Uxj=Uj)
    xt= x0+v0*t
    λ2 = max((ForwardDiff.derivative(x -> Uxj([xt[1], x, xt[3], xt[4]]; j=j), xt[2])+
            ∂-mean(∂j[j]))*v0[2], 0.0)
    return λ2
end


function rateswitchj_d3(t; j, x0, v0, ∂=∇tot[3], ∂j=∇j[:,3],  Uxj=Uj)
    xt= x0+v0*t
    λ3 = max((ForwardDiff.derivative(x -> Uxj([xt[1], xt[2], x, xt[4]]; j=j), xt[3])+
            ∂-mean(∂j[j]))*v0[3], 0.0)
    return λ3
end


function rateswitchj_d4(t; j, x0, v0, ∂=∇tot[4], ∂j=∇j[:,4],  Uxj=Uj)
    xt= x0+v0*t
    λ4 = max((ForwardDiff.derivative(x -> Uxj([xt[1], xt[2], xt[3], x]; j=j), xt[4])+
            ∂-mean(∂j[j]))*v0[4], 0.0)
    return λ4
end

@benchmark rateswitchj(0.01; j=[3, 5, 6], x0=[med1, med2, med3, med4], v0=ones(4))[1]

print(rateswitchj(0.01; j=[2], x0=[med1, med2, med3, med4], v0=ones(4)))
print(rateswitchj_d1(0.01; j=[2], x0=[med1, med2, med3, med4], v0=ones(4)))
print(rateswitchj_d2(0.01; j=[2], x0=[med1, med2, med3, med4], v0=ones(4)))
print(rateswitchj_d3(0.01; j=[2], x0=[med1, med2, med3, med4], v0=ones(4)))
print(rateswitchj_d4(0.01; j=[2], x0=[med1, med2, med3, med4], v0=ones(4)))
print(rateswitchj(0.01; j=[3, 6, 10], x0=[-0.5, +9, -0.4, -0.4], v0=[-1,1,1,-1]))
print(rateswitchj_d1(0.01; j=[3, 6, 10], x0=[-0.5, +9, -0.4, -0.4], v0=[-1,1,1,-1]))
print(rateswitchj_d2(0.01; j=[3, 6, 10], x0=[-0.5, +9, -0.4, -0.4], v0=[-1,1,1,-1]))
print(rateswitchj_d3(0.01; j=[3, 6, 10], x0=[-0.5, +9, -0.4, -0.4], v0=[-1,1,1,-1]))
print(rateswitchj_d4(0.01; j=[3, 6, 10], x0=[-0.5, +9, -0.4, -0.4], v0=[-1,1,1,-1]))

@benchmark rateswitchj(0.01; j=[3, 5, 6], x0=[med1, med2, med3, med4], v0=ones(4))[1]
@benchmark rateswitchj_d1(0.01; j=[3, 5, 6], x0=[med1, med2, med3, med4], v0=ones(4))





zzskss_5=zz_w_ss(; NS=5000, x0_0=x̃, tmax=tmax_tuned,ssM =1000,
        NOBS=length(t_vec), ssS=5)
zzskss_10=zz_w_ss(; NS=5000, x0_0=x̃, tmax=tmax_tuned,ssM =1000,
        NOBS=length(t_vec), ssS=10)
zzskss_10bis=zz_w_ss(; NS=5000, x0_0=x̃, tmax=tmax_tuned,ssM =1000,
        NOBS=length(t_vec), ssS=10)
zzskss_10tris=zz_w_ss(; NS=5000, x0_0=x̃, tmax=tmax_tuned,ssM =1000,
        NOBS=length(t_vec), ssS=10)
zzskss_50=zz_w_ss(; NS=5000, x0_0=x̃, tmax=tmax_tuned,ssM =1000,
        NOBS=length(t_vec), ssS=50)
zzskss_100=zz_w_ss(; NS=5000, x0_0=x̃, tmax=tmax_tuned,ssM =1000,
        NOBS=length(t_vec), ssS=100)
zzskss_200=zz_w_ss(; NS=5000, x0_0=x̃, tmax=tmax_tuned,ssM =1000,
        NOBS=length(t_vec), ssS=200)
zzskss_1000=zz_w_ss(; NS=5000, x0_0=x̃, tmax=tmax_tuned,ssM =1000,
        NOBS=length(t_vec), ssS=1000)

JLD.save(string(SAVEwd, "/zzSSs.jld"),
    Dict("zzskss_5"=>zzskss_5, "zzskss_10"=>zzskss_10,
    "zzskss_50"=>zzskss_50, "zzskss_100"=>zzskss_100,
    "zzskss_200"=>zzskss_200, "zzskss_1000"=>zzskss_1000))


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

smpALLI  = zzsample(N=2000, sk=zzsk)
smpALLII = zzsample(N=2000, sk=zzsk2)
sssmp_5  = zzsample(N=2000, sk=zzskss_5)
sssmp_10  = zzsample(N=2000, sk=zzskss_10)
sssmp_10bis  = zzsample(N=2000, sk=zzskss_10bis)
sssmp_10tris  = zzsample(N=2000, sk=zzskss_10tris)
sssmp_50  = zzsample(N=2000, sk=zzskss_50)
sssmp_100  = zzsample(N=2000, sk=zzskss_100)
sssmp_200  = zzsample(N=2000, sk=zzskss_200)
sssmp_1000  = zzsample(N=2000, sk=zzskss_1000)




pm_dict=Dict()
for d in 1:Dim
    pm=density(smpALLI[1001:2000, d],legend=true,linewidth=1, label="All obs", color=:black)
    density!(smpALLII[1001:2000, d],legend=true,linewidth=1, color=:black)
    density!(sssmp_5[:, d],label="5", linewidth=0.5, color=:red)
    density!(sssmp_10[:, d], label="10 ", linewidth=0.5, color=:blue)
    density!(sssmp_10bis[:, d], label="10bis", linewidth=0.5, color=:lightblue)
    density!(sssmp_10tris[:, d], label="10tris", linewidth=0.5, color=:darkblue)
    density!(sssmp_50[:, d],label="50", linewidth=0.5, color=:green)
    density!(sssmp_100[:, d],label="100", linewidth=0.5, color=:yellow)
    density!(sssmp_200[:, d],label="200 ", linewidth=0.5, color=:orange)
    density!(sssmp_1000[:, d],label="1000", linewidth=0.5, color=:purple)
    pm_dict[d]=pm
end

pmsm=plot(pm_dict[1], pm_dict[2], pm_dict[3], pm_dict[4], layout=(2,2), size=(1000,800))
savefig(pmsm, string(SAVEwd, "pmsm.pdf"))
