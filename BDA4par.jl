# Direcotries
CODEwd = "/Users/alice/projects/ZZpaper/"
SAVEwd = "/Users/alice/OneDrive - University of Warwick/Manuscripts/04ADZZ/figures/sec5.2/"

# Loading the packages
include(string(CODEwd, "env.jl"))
# Loading the functions for general zig zag
include(string(CODEwd, "functions.jl"))
include(string(CODEwd, "plotting.jl"))


pt(78)
gr(size = (pt(84), pt(84)), labelfontsize=8, legend = false)

r_tmax=100
r_insp=10
r_ess=100
budget=200000


# Include supplementary packages for data IO
using CSV
using DataFrames
using StatsBase
using Serialization
using Dates
using Extremes

Y = open(deserialize, string(CODEwd,"data/BDA_preprocessed.bin"))
Random.seed!(1234)
Ysd = Y[sample(1:(size(Y)[1]), 500), :]
print(Ysd[1:20, :])

# model with age and STAGE if spreaded
dimnames=["lgα", "β₀", "β₁", "β₂"]
t_vec = Dates.value.(Ysd.TIMEDIAGNOSITOFINALEVENT).+0.01
c_vec = zeros(length(t_vec))
c_vec[findall(Ysd.NEWVITALSTATUS.=="D")] .= 1
z1_vec = (Ysd.AGE.-mean(Ysd.AGE))./std(Ysd.AGE)
z2_vec = zeros(length(t_vec))
z2_vec[findall(Ysd.STAGENUMBER.>=3)] .= 1

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

opti = optimize(U, [-0.3, 10.0, 0, -1.5])
x̃=opti.minimizer
# profile likelihoods
r_halfwidth= [2,10, 3,3]
plot(x -> U(vcat(x, x̃[2:4])), xlims=(x̃[1]-r_halfwidth[1], x̃[1]+r_halfwidth[1]))
plot(x -> U(vcat(x̃[1],x, x̃[3:4])), xlims=(x̃[2]-r_halfwidth[2], x̃[2]+r_halfwidth[2]))
plot(x -> U(vcat(x̃[1:2],x, x̃[4])), xlims=(x̃[3]-r_halfwidth[3], x̃[3]+r_halfwidth[3]))
plot(x -> U(vcat(x̃[1:3], x)), xlims=(x̃[4]-r_halfwidth[4], x̃[4]+r_halfwidth[4]))

# profile likelihoods on larger ranges
r_halfwidth= [2,10, 3,3]*5
plot(x -> U(vcat(x, x̃[2:4])), xlims=(x̃[1]-r_halfwidth[1], x̃[1]+r_halfwidth[1]))
plot(x -> U(vcat(x̃[1],x, x̃[3:4])), xlims=(x̃[2]-r_halfwidth[2], x̃[2]+r_halfwidth[2]))
plot(x -> U(vcat(x̃[1:2],x, x̃[4])), xlims=(x̃[3]-r_halfwidth[3], x̃[3]+r_halfwidth[3]))
plot(x -> U(vcat(x̃[1:3], x)), xlims=(x̃[4]-r_halfwidth[4], x̃[4]+r_halfwidth[4]))


contour(range(-.5, -.45, length=100),range(9.3,9.8, length=100),
    [(U([x,y, x̃[3], x̃[4]])) for x in range(-.5, -.45, length=100),
        y in range(9.3,9.8, length=100)], width=1,
        legend=true, xlabel=dimnames[1], ylabel=dimnames[2])


contour(range(-.5, -.45, length=100),range(-.75,-0.7, length=100),
    [(U([x,x̃[2],y,x̃[4]])) for x in range(-.5, -.45, length=100),
        y in range(-.75,-0.7, length=100)], width=1,
        legend=true, xlabel=dimnames[1], ylabel=dimnames[3])

contour(range(-.5, -.45, length=100),range(-3,-2, length=100),
    [(U([x,x̃[2], x̃[3],y])) for x in range(-.5, -.45, length=100),
        y in range(-3,-2, length=100)], width=1,
        legend=true, xlabel=dimnames[1], ylabel=dimnames[4])




U([-0.5, +9, -0.4, -0.4])
U([4, 10, 5, 5])
x=[4, 10, 5, 5]
Dim  =4
start = [-0.5, +9, -0.4, -0.4]
start2= [-0.15, 13, -0.1, -5]
zzsk  =  zz(; NS=1000, x0_0=start, tmax=0.02)
zzsk2 =  zz(; NS=1000, x0_0=start2, tmax=0.02)

p1=plot(zzsk["SK"][:, 1], zzsk["SK"][:, 2], ylabel=dimnames[1])
plot!(zzsk2["SK"][:, 1], zzsk2["SK"][:, 2])
p2=plot(zzsk["SK"][:, 1], zzsk["SK"][:, 3], ylabel=dimnames[2])
plot!(zzsk2["SK"][:, 1], zzsk2["SK"][:, 3])
p3=plot(zzsk["SK"][:, 1], zzsk["SK"][:, 4], ylabel=dimnames[3])
plot!(zzsk2["SK"][:, 1], zzsk2["SK"][:, 4])
p4=plot(zzsk["SK"][:, 1], zzsk["SK"][:, 5], ylabel=dimnames[4])
plot!(zzsk2["SK"][:, 1], zzsk2["SK"][:, 5])
p5 = plot(p1,p2,p3,p4, layout=(2,2), size=(600,600))



f2=tmplot(R=r_tmax, try_tmax=[0.05, 0.075, 0.1, 0.15, 0.2],
    start=[-0.41, 9.5, -0.35, -2.1])
savefig(f2["NHPP"], string(SAVEwd, "f2a.pdf"))
savefig(f2["opt"], string(SAVEwd, "f2b.pdf"))
savefig(f2["tot"], string(SAVEwd, "f2c.pdf"))
f2tot= plot(f2["NHPP"],f2["opt"],f2["tot"], layout=(1,3),
    size=(pt(252), pt(60)))
savefig(f2tot, string(SAVEwd, "f2.pdf"))
tmax_tuned=0.1


start2= [-0.15, 13, -0.1, -5]
zzsk_plot  =  zz(; NS=2000, x0_0=start2, tmax=tmax_tuned)

p1=plot(zzsk_plot["SK"][:, 1], zzsk_plot["SK"][:, 2], ylabel=dimnames[1], color=:gray, linewidth=0.5)
p2=plot(zzsk_plot["SK"][:, 1], zzsk_plot["SK"][:, 3], ylabel=dimnames[2], color=:gray, linewidth=0.5)
p3=plot(zzsk_plot["SK"][:, 1], zzsk_plot["SK"][:, 4], ylabel=dimnames[3], color=:gray, linewidth=0.5)
p4=plot(zzsk_plot["SK"][:, 1], zzsk_plot["SK"][:, 5], ylabel=dimnames[4], color=:gray, linewidth=0.5)
p5 = plot(p1,p2,p3,p4, layout=(2,2), size=(pt(84*2),pt(84*2)))
savefig(p5, string(SAVEwd, "out4par1skel.pdf"))

# inspection of the zz from mode

# f3 - inspection of the zig zag from mode
start_μ=[-0.45,9.5,-0.8,-2]
# multiple chains from the center
SKinsp= Array{Float64, 3}(undef, 1000, Dim, r_insp)
TMinsp= Array{Float64, 3}(undef, 1000, 1, r_insp)
for r in 1:r_insp
    start= rand(MvNormal(start_μ,Diagonal(ones(Dim)).*0.01+zeros( Dim,Dim)))
    zzsk =  zz(; NS=1000, x0_0=start, tmax=tmax_tuned)
    SKinsp[:,:,r] = zzsk["SK"][:, 2:(Dim+1)]
    TMinsp[:,:,r] = zzsk["SK"][:, 1]
end

f3a=plot(0,0, alpha=0, xlabel="Time", ylabel=dimnames[1])
for r in 1:r_insp
    plot!(TMinsp[:, 1, r], SKinsp[:,1,r], linewidth=0.5, color=:gray)
end
f3b=plot(0,0, alpha=0, xlabel="Time", ylabel=dimnames[2])
for r in 1:r_insp
    plot!(TMinsp[:, 1, r], SKinsp[:,2,r], linewidth=0.5, color=:gray)
end
f3c=plot(0,0, alpha=0, xlabel="Time", ylabel=dimnames[3])
for r in 1:r_insp
    plot!(TMinsp[:, 1, r], SKinsp[:,3,r], linewidth=0.5, color=:gray)
end
f3d=plot(0,0, alpha=0, xlabel="Time", ylabel=dimnames[4])
for r in 1:r_insp
    plot!(TMinsp[:, 1, r], SKinsp[:,4,r], linewidth=0.5, color=:gray)
end
savefig(f3a, string(SAVEwd, "/f3a.pdf"))
savefig(f3b, string(SAVEwd, "/f3b.pdf"))
savefig(f3c, string(SAVEwd, "/f3c.pdf"))
savefig(f3d, string(SAVEwd, "/f3d.pdf"))
f3 = plot(f3a, f3b, f3c, f3d, layout=(2, 2), size=(pt(84*2), pt(84*2)))
savefig(f3, string(SAVEwd, "/f3.pdf"))

Random.seed!(1234)
SKsize=1000
SKinsp= Array{Float64, 3}(undef, SKsize, Dim, r_insp)
TMinsp= Array{Float64, 3}(undef, SKsize, 1, r_insp)
for r in 1:r_insp
    start=[sample([-5,2]), sample([-1,40]),
            sample([-4,3]), sample([-5,1])]
    print(start)
    zzsk =  zz(; NS=SKsize, x0_0=start, tmax=tmax_tuned)
    SKinsp[:,:,r] = zzsk["SK"][:, 2:(Dim+1)]
    TMinsp[:,:,r] = zzsk["SK"][:, 1]
end


f4a=plot(0,0, alpha=0, xlabel="Time", ylabel=dimnames[1])
for r in 1:r_insp
    plot!(TMinsp[:, 1, r], SKinsp[:,1,r], linewidth=0.5, color=:gray)
end
f4b=plot(0,0, alpha=0, xlabel="Time", ylabel=dimnames[2])
for r in 1:r_insp
    plot!(TMinsp[:, 1, r], SKinsp[:,2,r], linewidth=0.5, color=:gray)
end
f4c=plot(0,0, alpha=0, xlabel="Time", ylabel=dimnames[3])
for r in 1:r_insp
    plot!(TMinsp[:, 1, r], SKinsp[:,3,r], linewidth=0.5, color=:gray)
end
f4d=plot(0,0, alpha=0, xlabel="Time", ylabel=dimnames[4])
for r in 1:r_insp
    plot!(TMinsp[:, 1, r], SKinsp[:,4,r], linewidth=0.5, color=:gray)
end
f4 = plot(f4a, f4b, f4c, f4d, layout=(2, 2), size=(pt(84*2), pt(84*2)))
savefig(f4, string(SAVEwd, "/f4.pdf"))

# f5 - density inspection
zzsk = zz(; NS=100000, x0_0=start_μ, tmax=tmax_tuned)
zzsm = zzsample(;N=10000, sk=zzsk)
f5a=density(zzsm[1000:10000, 1], xlabel=dimnames[1], color=:black, linewidth=3, alpha=0.7)
f5b=density(zzsm[1000:10000, 2], xlabel=dimnames[2], color=:black, linewidth=3, alpha=0.7)
f5c=density(zzsm[1000:10000, 3], xlabel=dimnames[3], color=:black, linewidth=3, alpha=0.7)
f5d=density(zzsm[1000:10000, 4], xlabel=dimnames[4], color=:black, linewidth=3, alpha=0.7)
f5 = plot(f5a, f5b, f5c, f5d, layout=(2, 2), size=(pt(84*2), pt(84*2)))
savefig(f5, string(SAVEwd, "/f5.pdf"))

# tuning the HMC


# parameters HMC
Lε_tuned=0.08
L_tuned=2

hmc = runHMC(;epsilon=Lε_tuned/L_tuned,L=L_tuned,IT=1000,qs=start_μ)
plot(hmc["SampleQ"][:, 1], title=dimnames[1])
plot(hmc["SampleQ"][:, 2], title=dimnames[2])
plot(hmc["SampleQ"][:, 1], title=dimnames[3])
plot(hmc["SampleQ"][:, 1], title=dimnames[4])


plot(autocor(hmc["SampleQ"][:, 1]))
plot!(autocor(hmc["SampleQ"][:, 2]))
plot!(autocor(hmc["SampleQ"][:, 3]))
plot!(autocor(hmc["SampleQ"][:, 4]))


sum(hmc["accept"])/1000


# f6 - inspection of the HMC from mode
# multiple chains from the center
# f6 - inspection of the HMC from mode
start_μ=[-0.45,9.5,-0.8,-2]
# multiple chains from the center
HMCinsp= Array{Float64, 3}(undef, 1000, Dim, r_insp)
for r in 1:r_insp
    start= rand(MvNormal(start_μ,Diagonal(ones(Dim)).*0.01+zeros( Dim,Dim)))
    hmc = runHMC(;epsilon=Lε_tuned/L_tuned,L=L_tuned,IT=1000,qs=start)
    HMCinsp[:,:,r] = hmc["SampleQ"][1:1000, 1:(Dim)]
end

f6a=plot(0,0, alpha=0, xlabel="Iterations", ylabel=dimnames[1])
for r in 1:r_insp
    plot!(HMCinsp[:,1,r], linewidth=0.5, color=:gray)
end
f6b=plot(0,0, alpha=0, xlabel="Iterations", ylabel=dimnames[2])
for r in 1:r_insp
    plot!(HMCinsp[:,2,r], linewidth=0.5, color=:gray)
end
f6c=plot(0,0, alpha=0, xlabel="Iterations", ylabel=dimnames[3])
for r in 1:r_insp
    plot!(HMCinsp[:,3,r], linewidth=0.5, color=:gray)
end
f6d=plot(0,0, alpha=0, xlabel="Iterations", ylabel=dimnames[4])
for r in 1:r_insp
    plot!(HMCinsp[:,4,r], linewidth=0.5, color=:gray)
end
f6 = plot(f6a, f6b, f6c, f6d, layout=(2, 2), size=(pt(84*2), pt(84*2)))
savefig(f6, string(SAVEwd, "/f6.pdf"))

# inspection of the hmc from far (same values as before)
Random.seed!(1234)
HMCinsp= Array{Float64, 3}(undef, 10000, Dim, r_insp)
for r in 1: r_insp
    start=[sample([-5,2]), sample([-1,40]),
            sample([-4,3]), sample([-5,1])]
    hmc = runHMC(;epsilon=Lε_tuned/L_tuned,L=L_tuned,IT=10000,qs=start)
    HMCinsp[:,:,r] = hmc["SampleQ"][1:10000, 1:(Dim)]
end

HMCinsp

f7a=plot(0,0, alpha=0, xlabel="Iterations", ylabel=dimnames[1])
for r in 1:r_insp
    plot!(HMCinsp[:,1,r], linewidth=0.5, color=:gray)
end
f7b=plot(0,0, alpha=0, xlabel="Iterations", ylabel=dimnames[2])
for r in 1:r_insp
    plot!(HMCinsp[:,2,r], linewidth=0.5, color=:gray)
end
f7c=plot(0,0, alpha=0, xlabel="Iterations", ylabel=dimnames[3])
for r in 1:r_insp
    plot!(HMCinsp[:,3,r], linewidth=0.5, color=:gray)
end
f7d=plot(0,0, alpha=0, xlabel="Iterations", ylabel=dimnames[4])
for r in 1:r_insp
    plot!(HMCinsp[:,4,r], linewidth=0.5, color=:gray)
end
f7 = plot(f7a, f7b, f7c, f7d, layout=(2, 2), size=(pt(84*2), pt(84*2)))
savefig(f7, string(SAVEwd, "/f7.pdf"))


# f8 - density hmc
hmc = runHMC(;epsilon=Lε_tuned/L_tuned,L=L_tuned,IT=100000,qs=start_μ)
f8a=density(hmc["SampleQ"][1000:10:100000, 1], xlabel=dimnames[1], color=:black, linewidth=3, alpha=0.7)
f8b=density(hmc["SampleQ"][1000:10:100000, 2], xlabel=dimnames[2], color=:black, linewidth=3, alpha=0.7)
f8c=density(hmc["SampleQ"][1000:10:100000, 3], xlabel=dimnames[3], color=:black, linewidth=3, alpha=0.7)
f8d=density(hmc["SampleQ"][1000:10:100000, 4], xlabel=dimnames[4], color=:black, linewidth=3, alpha=0.7)
f8 = plot(f8a, f8b, f8c, f8d, layout=(2, 2), size=(pt(84*2), pt(84*2)))

savefig(f8, string(SAVEwd, "/f8.pdf"))


# f8 - density hmc
f58a=density(hmc["SampleQ"][1000:10:100000, 1], xlabel=dimnames[1], color=:red, linewidth=3, alpha=0.7)
density!(zzsm[1000:10000, 1],  color=:blue, linewidth=3, alpha=0.7)
f58b=density(hmc["SampleQ"][1000:10:100000, 2], xlabel=dimnames[2], color=:red, linewidth=3, alpha=0.7)
density!(zzsm[1000:10000, 2], color=:blue, linewidth=3, alpha=0.7)
f58c=density(hmc["SampleQ"][1000:10:100000, 3], xlabel=dimnames[3], color=:red, linewidth=3, alpha=0.7)
density!(zzsm[1000:10000, 3],  color=:blue, linewidth=3, alpha=0.7)
f58d=density(hmc["SampleQ"][1000:10:100000, 4], xlabel=dimnames[4], color=:red, linewidth=3, alpha=0.7)
density!(zzsm[1000:10000, 4], color=:blue, linewidth=3, alpha=0.7)
f58 = plot(f58a, f58b, f58c, f58d, layout=(2, 2), size=(pt(84*2), pt(84*2)))
savefig(f58, string(SAVEwd, "/f58.pdf"))

r_ess=100
zzs = Dict()
hmcs= Dict()
for r in 1:r_ess
    start=rand(MvNormal(start_μ,Diagonal(ones(Dim)).*0.01+zeros( Dim,Dim)))
    zzs[[r]]=zz(NS=false, B=budget, x0_0=start, tmax=tmax_tuned)
    hmcs[[r]]=runHMC(;epsilon=Lε_tuned/L_tuned,L=L_tuned,IT=Int(round(budget/(L_tuned+1))),qs=start)
    print(r)
end

# f9 - chose nbatches for zz
f9= bs_zz(;try_nb=[5, 10,25, 50, 100, 500])
savefig(f9, string(SAVEwd, "/f9.pdf"))
nbZZ_tuned = 100

# f10 - chose nbatches for hmc
f10= bs_hmc(;try_nb=[5, 10,25, 50, 100, 500])
savefig(f10, string(SAVEwd, "/f10.pdf"))
nbHMC_tuned = 100

# summaries and f11 - plot of the ESS
outESS = ESSsummaries()
outESS


touch(string(SAVEwd,"/summaries.txt" ))
io = open(string(SAVEwd, "/summaries.txt"), "w")
write(io, outESS["string"])
close(io)
f11=violin(string.(transpose(1:Dim)), (outESS["essZZ"]), side=:right, linewidth=0,
    label="ZZ", color=:black, alpha=0.7, xlabel="Dim", ylabel="ESS", legend=false,
    ylims=(0, maximum(hcat(outESS["essZZ"], outESS["essHMC"]))*1.1))
violin!(string.(transpose(1:Dim)), (outESS["essHMC"]), side=:left, linewidth=0,
    label="HMC", color=:lightgray)
annotate!([(Dim+0.2, maximum(hcat(outESS["essZZ"], outESS["essHMC"]))*1.05, ("ZZ", 14, :black, :center)),
    (Dim+0.2, maximum(hcat(outESS["essZZ"], outESS["essHMC"]))*0.95, ("HMC", 14, :lightgray, :center))])
savefig(f11, string(SAVEwd, "/f11.pdf"))

# save ESS
JLD.save(string(SAVEwd, "/out.jld"),
    Dict("outESS"=>outESS, "zzs"=>zzs, "hmcs"=>hmcs))
