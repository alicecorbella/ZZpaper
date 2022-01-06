# Direcotries
CODEwd = "/Users/alice/projects/ZZpaper/"
SAVEwd = "/Users/alice/OneDrive - University of Warwick/Manuscripts/04ADZZ/figures/Sec5.1/"

# Loading the packages
include(string(CODEwd, "env.jl"))

# Loading the functions
include(string(CODEwd, "functions.jl"))
include(string(CODEwd, "plotting.jl"))

r_insp=10
r_ess=100
budget=200000

Dim=4
data_x = [1.0, 1.5, 1.5, 1.5, 2.5, 4.0, 5.0, 5.0, 7.0, 8.0, 8.5, 9.0, 9.5,
    9.5, 10, 12, 12, 13, 13, 14.5, 15.5, 15.5, 16.5, 17, 22.5, 29.0, 31.5]
data_y = [1.80, 1.85, 1.87, 1.77, 2.02, 2.27, 2.15, 2.26, 2.47, 2.19, 2.26,
    2.40, 2.39, 2.41, 2.50, 2.32, 2.32, 2.43, 2.47, 2.56, 2.65, 2.47, 2.64,
    2.56, 2.70, 2.72, 2.57]

dimnames=["α", "β", "γ", "σ"]


function U(x::Vector; xdat=data_x, ydat=data_y)
    α = exp(x[1])
    β = exp(x[2])
    γ = exp(x[3])/(1+exp(x[3]))
    σ = exp(x[4])
    n = size(xdat)[1]
    # minus log prior: all uniform between 0 and 100, only the jacobian matters
    mlp = -(ifelse(α>0, 0, -Inf)+ ifelse(β>0, 0, -Inf)+
        logpdf(Beta(7,7/3), γ) +ifelse(σ>0, 0, -Inf)+
        x[1]+x[2]+x[3]-2*log(1+exp(x[3]))+x[4]
        )
    # minus log posterior is the mlp plus the constant plus the
    # loop over all the variables
    mll = mlp + (n/2)*log(2*pi*σ^2)

    for i in 1:n
        mll+= 0.5*((ydat[i]-α+β*(γ^xdat[i]))/σ)^2
    end
    return mll
end


# f2 - data
f1=plot(data_x, data_y, seriestype=:scatter, legend=false, color=:black,
    alpha=0.5, xlabel="Xⱼ", ylabel="Yⱼ", size = (300, 240))
savefig(f1, string(SAVEwd, "f1.pdf"))

# f2 - tune of t_max
f2=tmplot(R=100, try_tmax=[0.01, 0.02, 0.03,0.04, 0.05, 0.06])
savefig(f2, string(SAVEwd,"f2.pdf"))
tmax_tuned=0.02

# f3 - inspection of the zig zag from mode
start=[0,0,0,0]
# multiple chains from the center
SKinsp= Array{Float64, 3}(undef, 1000, Dim, r_insp)
TMinsp= Array{Float64, 3}(undef, 1000, 1, r_insp)
for r in 1:r_insp
    zzsk =  zz(; NS=1000, x0_0=start, tmax=tmax_tuned)
    SKinsp[:,:,r] = zzsk["SK"][:, 2:(Dim+1)]
    TMinsp[:,:,r] = zzsk["SK"][:, 1]
end

f3a=plot(0,0, alpha=0, xlabel="Time", ylabel="θ₁")
for r in 1:r_insp
    plot!(TMinsp[:, 1, r], SKinsp[:,1,r])
end
f3b=plot(0,0, alpha=0, xlabel="Time", ylabel="θ₂")
for r in 1:r_insp
    plot!(TMinsp[:, 1, r], SKinsp[:,2,r])
end
f3c=plot(0,0, alpha=0, xlabel="Time", ylabel="θ₃")
for r in 1:r_insp
    plot!(TMinsp[:, 1, r], SKinsp[:,3,r])
end
f3d=plot(0,0, alpha=0, xlabel="Time", ylabel="θ₄")
for r in 1:r_insp
    plot!(TMinsp[:, 1, r], SKinsp[:,4,r])
end
f3 = plot(f3a, f3b, f3c, f3d, layout=(2, 2), size=(1000, 800))
savefig(f3, string(SAVEwd, "/f3.pdf"))

# f4 - inspection of the zig zag from tail
startvalues=hcat([5, 5, 5, 5],
                [-5, 5, 5, 5],[5, -5, 5, 5],[5, 5, -5, 5],[5, 5, 5, -5],
                [-5, -5, 5, 5],[-5, 5, -5, 5],[-5, 5, 5, -5],
                [5, -5, -5, 5],[5, -5, 5, -5],[5, 5, -5, -5],
                [-5, -5, -5, 5],[-5, -5, 5, -5],[-5, 5, -5, -5],[5, -5, -5, -5],
                [-5, -5, -5, -5])
SKinsp= Array{Float64, 3}(undef, 1000, Dim, size(startvalues)[2])
TMinsp= Array{Float64, 3}(undef, 1000, 1, size(startvalues)[2])
for r in 1:size(startvalues)[2]
    print(U(startvalues[:, r]))
    zzsk =  zz(; NS=1000, x0_0=startvalues[:, r], tmax=tmax_tuned)
    SKinsp[:,:,r] = zzsk["SK"][:, 2:(Dim+1)]
    TMinsp[:,:,r] = zzsk["SK"][:, 1]
    print(r)
end

f4a=plot(0,0, alpha=0, xlabel="Time", ylabel="θ₁")
for r in 1:r_insp
    plot!(TMinsp[:, 1, r], SKinsp[:,1,r])
end
f4b=plot(0,0, alpha=0, xlabel="Time", ylabel="θ₂")
for r in 1:r_insp
    plot!(TMinsp[:, 1, r], SKinsp[:,2,r])
end
f4c=plot(0,0, alpha=0, xlabel="Time", ylabel="θ₃")
for r in 1:r_insp
    plot!(TMinsp[:, 1, r], SKinsp[:,3,r])
end
f4d=plot(0,0, alpha=0, xlabel="Time", ylabel="θ₄")
for r in 1:r_insp
    plot!(TMinsp[:, 1, r], SKinsp[:,4,r])
end
f4 = plot(f4a, f4b, f4c, f4d, layout=(2, 2), size=(1000, 800))
savefig(f4, string(SAVEwd, "/f4.pdf"))


# f5 - density inspection
zzsk = zz(; NS=100000, x0_0=start, tmax=tmax_tuned)
zzsm = zzsample(;N=10000, sk=zzsk)
f5a=density(zzsm[1000:10000, 1], xlabel="θ₁", color=:black, linewidth=3, alpha=0.7)
f5b=density(zzsm[1000:10000, 2], xlabel="θ₂", color=:black, linewidth=3, alpha=0.7)
f5c=density(zzsm[1000:10000, 3], xlabel="θ₃", color=:black, linewidth=3, alpha=0.7)
f5d=density(zzsm[1000:10000, 4], xlabel="θ₄", color=:black, linewidth=3, alpha=0.7)
f5 = plot(f5a, f5b, f5c, f5d, layout=(2, 2), size=(1000, 800))
savefig(f5, string(SAVEwd, "/f5.pdf"))


# parameters HMC
Lε_tuned=0.1
L_tuned=10



# f6 - inspection of the HMC from mode
start=[0,0,0,0]
# multiple chains from the center
HMCinsp= Array{Float64, 3}(undef, 1000, Dim, r_insp)
for r in 1:r_insp
    hmc = runHMC(;epsilon=Lε_tuned/L_tuned,L=L_tuned,IT=1000,qs=start)
    HMCinsp[:,:,r] = hmc["SampleQ"][1:1000, 1:(Dim)]
end

f6a=plot(0,0, alpha=0, xlabel="Iterations", ylabel="θ₁")
for r in 1:r_insp
    plot!(HMCinsp[:,1,r])
end
f6b=plot(0,0, alpha=0, xlabel="Iterations", ylabel="θ₂")
for r in 1:r_insp
    plot!(HMCinsp[:,2,r])
end
f6c=plot(0,0, alpha=0, xlabel="Iterations", ylabel="θ₃")
for r in 1:r_insp
    plot!(HMCinsp[:,3,r])
end
f6d=plot(0,0, alpha=0, xlabel="Iterations", ylabel="θ₄")
for r in 1:r_insp
    plot!(HMCinsp[:,4,r])
end
f6 = plot(f6a, f6b, f6c, f6d, layout=(2, 2), size=(1000, 800))
savefig(f6, string(SAVEwd, "/f6.pdf"))


# f7 - inspection of the HMC from far

startvalues=hcat([5, 5, 5, 5],
                [-5, 5, 5, 5],[5, -5, 5, 5],[5, 5, -5, 5],[5, 5, 5, -5],
                [-5, -5, 5, 5],[-5, 5, -5, 5],[-5, 5, 5, -5],
                [5, -5, -5, 5],[5, -5, 5, -5],[5, 5, -5, -5],
                [-5, -5, -5, 5],[-5, -5, 5, -5],[-5, 5, -5, -5],[5, -5, -5, -5],
                [-5, -5, -5, -5])

HMCinsp= Array{Float64, 3}(undef, 1000, Dim, size(startvalues)[2])
for r in 1: size(startvalues)[2]
    hmc = runHMC(;epsilon=Lε_tuned/L_tuned,L=L_tuned,IT=1000,qs=startvalues[:,r])
    HMCinsp[:,:,r] = hmc["SampleQ"][1:1000, 1:(Dim)]
end

HMCinsp

f7a=plot(0,0, alpha=0, xlabel="Time", ylabel="θ₁")
for r in 1:r_insp
    plot!(HMCinsp[:,1,r])
end
f7b=plot(0,0, alpha=0, xlabel="Time", ylabel="θ₂")
for r in 1:r_insp
    plot!(HMCinsp[:,2,r])
end
f7c=plot(0,0, alpha=0, xlabel="Time", ylabel="θ₃")
for r in 1:r_insp
    plot!(HMCinsp[:,3,r])
end
f7d=plot(0,0, alpha=0, xlabel="Time", ylabel="θ₄")
for r in 1:r_insp
    plot!(HMCinsp[:,4,r])
end
f7 = plot(f7a, f7b, f7c, f7d, layout=(2, 2), size=(1000, 800))
savefig(f7, string(SAVEwd, "/f7.pdf"))

# f8 - density hmc
hmc = runHMC(;epsilon=Lε_tuned/L_tuned,L=L_tuned,IT=100000,qs=start)
f8a=density(hmc["SampleQ"][1000:10:100000, 1], xlabel="θ₁", color=:black, linewidth=3, alpha=0.7)
f8b=density(hmc["SampleQ"][1000:10:100000, 2], xlabel="θ₂", color=:black, linewidth=3, alpha=0.7)
f8c=density(hmc["SampleQ"][1000:10:100000, 3], xlabel="θ₃", color=:black, linewidth=3, alpha=0.7)
f8d=density(hmc["SampleQ"][1000:10:100000, 4], xlabel="θ₄", color=:black, linewidth=3, alpha=0.7)
f8 = plot(f8a, f8b, f8c, f8d, layout=(2, 2), size=(1000, 800))

savefig(f8, string(SAVEwd, "/f8.pdf"))



# f8 - density hmc
f58a=density(hmc["SampleQ"][1000:10:100000, 1], xlabel="θ₁", color=:red, linewidth=3, alpha=0.7)
density!(zzsm[1000:10000, 1], xlabel="θ₁", color=:blue, linewidth=3, alpha=0.7)
f58b=density(hmc["SampleQ"][1000:10:100000, 2], xlabel="θ₂", color=:red, linewidth=3, alpha=0.7)
density!(zzsm[1000:10000, 2], xlabel="θ₂", color=:blue, linewidth=3, alpha=0.7)
f58c=density(hmc["SampleQ"][1000:10:100000, 3], xlabel="θ₃", color=:red, linewidth=3, alpha=0.7)
density!(zzsm[1000:10000, 3], xlabel="θ₃", color=:blue, linewidth=3, alpha=0.7)
f58d=density(hmc["SampleQ"][1000:10:100000, 4], xlabel="θ₄", color=:red, linewidth=3, alpha=0.7)
density!(zzsm[1000:10000, 4], xlabel="θ₄", color=:blue, linewidth=3, alpha=0.7)
f58 = plot(f58a, f58b, f58c, f58d, layout=(2, 2), size=(1000, 800))
savefig(f58, string(SAVEwd, "/f58.pdf"))


zzs = Dict()
hmcs= Dict()
for r in 1:r_ess
    start=rand(Dim)
    zzs[[r]]=zz(NS=false, B=budget, x0_0=start, tmax=tmax_tuned)
    hmcs[[r]]=runHMC(;epsilon=Lε_tuned/L_tuned,L=L_tuned,IT=Int(round(budget/(L_tuned+1))),qs=start)
    print(r)
end

# f9 - chose nbatches for zz
f9= bs_zz(;try_nb=[5, 10,25, 50, 100, 500])
savefig(f9, string(SAVEwd, "/f9.pdf"))
nbZZ_tuned = 25

# f10 - chose nbatches for hmc
f10= bs_hmc(;try_nb=[5, 10,25, 50, 100, 500])
savefig(f10, string(SAVEwd, "/f10.pdf"))
nbHMC_tuned = 50

# summaries and f11 - plot of the ESS
outESS = ESSsummaries()
outESS


touch(string(SAVEwd,"/summaries.txt" ))
io = open(string(SAVEwd, "/summaries.txt"), "w")
write(io, outESS["string"])
close(io)
f11=violin(string.(transpose(1:Dim)), (outESS["essZZ"]), side=:right, linewidth=0,
    label="ZZ", color=:blue, alpha=0.7, xlabel="Dim", ylabel="ESS", legend=false,
    ylims=(0, maximum(hcat(outESS["essZZ"], outESS["essHMC"]))*1.1))
violin!(string.(transpose(1:Dim)), (outESS["essHMC"]), side=:left, linewidth=0,
    label="HMC", color=:red)
annotate!([(Dim+0.2, maximum(hcat(outESS["essZZ"], outESS["essHMC"]))*1.05, ("ZZ", 14, :blue, :center)),
    (Dim+0.2, maximum(hcat(outESS["essZZ"], outESS["essHMC"]))*0.95, ("HMC", 14, :red, :center))])
savefig(f11, string(SAVEwd, "/f11.pdf"))

# save ESS
JLD.save(string(SAVEwd, "/out.jld"),
    Dict("outESS"=>outESS, "zzs"=>zzs, "hmcs"=>hmcs))
