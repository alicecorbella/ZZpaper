# Direcotries
CODEwd = "/Users/alice/projects/ZZpaper/"
SAVEwd = "/Users/alice/OneDrive - University of Warwick/Manuscripts/04ADZZ/figures/performance/"

# Loading the packages
include(string(CODEwd, "env.jl"))

# Loading the functions
include(string(CODEwd, "functions.jl"))
include(string(CODEwd, "plotting.jl"))


r_tmax=100
r_insp=10
r_ess=100
budget=20000



# For large-sized journals the figures should be 84 mm
# (for double-column text areas), or 174 mm
# (for single-column text areas) wide and not higher than 234 mm.

pt(78)


gr(size = (pt(84), pt(84)), labelfontsize=8, legend = false,
    xtickfontsize=8,ytickfontsize=8)

# -----------------------------------------
# Example on Bivariate Isotropic normal
# -----------------------------------------
NameEx = "IsoNorm2d"
mkdir(string(SAVEwd, NameEx))
Dim = 2
μ = zeros(Dim)
Σ = Diagonal(ones(Dim))+zeros( Dim,Dim)
function U(x::Vector, m=μ, s=Σ)
    -(-2/2*log(2π)-0.5*log(det(s))-
    0.5*(transpose(x-m)*inv(s)*(x-m))
     )
end

# f1 - plot the target
f1=densplot()
savefig(f1, string(SAVEwd, NameEx, "/f1.pdf"))

# f2 - tune of t_max
f2=tmplot(R=r_tmax, try_tmax=[2, 2.25, 2.5, 2.75, 3, 3.5])
savefig(f2["NHPP"], string(SAVEwd, NameEx, "/f2a.pdf"))
savefig(f2["opt"], string(SAVEwd, NameEx, "/f2b.pdf"))
savefig(f2["tot"], string(SAVEwd, NameEx, "/f2c.pdf"))
f2tot= plot(f2["NHPP"],f2["opt"],f2["tot"], layout=(1,3),
    size=(pt(252), pt(60)))
savefig(f2tot, string(SAVEwd, NameEx, "/f2.pdf"))
tmax_tuned=2.5

# f3 - inspection of the zig zag from mode
# multiple chains from the center
f3=plot(0,0, alpha=0, xlabel="X₁", ylabel="X₂", size = (pt(84), pt(60)))
startpoints = Array{Float64, 2}(undef, r_insp, 2)
endpoints = Array{Float64, 2}(undef, r_insp, 2)
for r in 1:r_insp
    start= rand(MvNormal(μ,Σ))
    startpoints[r, :] = start
    zzsk =  zz(; NS=1000, x0_0=start, tmax=tmax_tuned)
    plot!(zzsk["SK"][:, 2], zzsk["SK"][:, 3], color=:gray, linewidth=0.3)
    endpoints[r, :] = zzsk["SK"][1000, 2:3]
end
plot!(startpoints[:,1], startpoints[:,2], seriestype=:scatter,
    markersize=2, color=:black)
plot!(endpoints[:,1], endpoints[:,2], seriestype=:scatter,
        markersize=2, markershape=:star5, color=:black)
display(f3)
savefig(f3, string(SAVEwd, NameEx, "/f3.pdf"))

# f4 - inspection of the zig zag from tail
qsize=0.000001
lim  = quantile(Normal(0,1), qsize)
f4=tailstart(xrange=[-lim, lim], yrange=[-lim, lim], nside=10)
savefig(f4, string(SAVEwd, NameEx, "/f4.pdf"))


# f5 - density inspection
start= rand(MvNormal(μ,Σ))
zzsk = zz(; NS=100000, x0_0=start, tmax=tmax_tuned)
zzsm = zzsample(;N=10000, sk=zzsk)
f5=plot(zzsm[:, 1], zzsm[:, 2], seriestype=:scatter,
        color=:black, alpha=0.5, markersize=2,
        xlabel="X₁", ylabel="X₂")
contour!(range(-4, 4, length=100),range(-4,4, length=100),[(pdf(MvNormal(μ, Σ), [x,y])) for y in range(-4, 4, length=100), x in range(-4,4, length=100)],
    fill=(false,cgrad(:grays,[0,0.1,1.0])), color=:lightgray, linewidth=.5)
savefig(f5, string(SAVEwd, NameEx, "/f5.pdf"))

# parameters HMC
Lε_tuned=1
L_tuned=2

# f6 - inspection of the HMC from mode
# multiple chains from the center
startpoints = Array{Float64, 2}(undef, r_insp, 2)
endpoints = Array{Float64, 2}(undef, r_insp, 2)

f6=plot(0,0, alpha=0, xlabel="X₁", ylabel="X₂", size = (pt(84), pt(60)))
for r in 1:r_insp
    start= rand(MvNormal(μ,Σ))
    startpoints[r, :] = start
    hmc = runHMC(;epsilon=Lε_tuned/L_tuned,L=L_tuned,IT=1000,qs=start)
    plot!(hmc["SampleQ"][:, 1], hmc["SampleQ"][:, 2], color=:gray, linewidth=0.3)
    endpoints[r, :] = hmc["SampleQ"][1000, 1:2]
end
plot!(startpoints[:,1], startpoints[:,2], seriestype=:scatter,
    markersize=2, color=:black)
plot!(endpoints[:,1], endpoints[:,2], seriestype=:scatter,
        markersize=2, markershape=:star5, color=:black)
savefig(f6, string(SAVEwd, NameEx, "/f6.pdf"))


# f7 - inspection of the HMC from far
f7=tailstart_hmc(xrange=[-lim, lim], yrange=[-lim, lim], nside=10)
savefig(f7, string(SAVEwd, NameEx, "/f7.pdf"))

# f8 - density hmc
hmc = runHMC(;epsilon=Lε_tuned/L_tuned,L=L_tuned,IT=100000,qs=start)
f8=plot(hmc["SampleQ"][1:10:100000, 1], hmc["SampleQ"][1:10:100000, 2], seriestype=:scatter,
        color=:black, alpha=0.5, markersize=2, xlabel="X₁", ylabel="X₂")
contour!(range(-4, 4, length=100),range(-4,4, length=100),
    [(pdf(MvNormal(μ, Σ), [x,y])) for y in range(-4, 4, length=100), x in range(-4,4, length=100)],
    fill=(false,cgrad(:grays,[0,0.1,1.0])), color=:lightgray, linewidth=0.5)
savefig(f8, string(SAVEwd, NameEx, "/f8.pdf"))

# -------- PERFORMANCE COMPARISON -------- #
# generate the sample
zzs = Dict()
hmcs= Dict()
for r in 1:r_ess
    start=rand(Dim)
    zzs[[r]]=zz(NS=false, B=budget, x0_0=start, tmax=tmax_tuned)
    hmcs[[r]]=runHMC(;epsilon=Lε_tuned/L_tuned,L=L_tuned,IT=Int(round(budget/(L_tuned+1))),qs=start)
    print(r)
end

# f9 - chose nbatches for zz
f9= bs_zz(;try_nb=[10, 100, 500, 1000, 2000])
savefig(f9, string(SAVEwd, NameEx, "/f9.pdf"))
nbZZ_tuned = 1000

# f10 - chose nbatches for hmc
f10= bs_hmc(;try_nb=[10, 100, 500, 1000, 2000])
savefig(f10, string(SAVEwd, NameEx, "/f10.pdf"))
nbHMC_tuned = 100

# summaries and f11 - plot of the ESS
outESS = ESSsummaries()
outESS


touch(string(SAVEwd, NameEx,"/summaries.txt" ))
io = open(string(SAVEwd, NameEx, "/summaries.txt"), "w")
write(io, outESS["string"])
close(io)
f11=violin(string.(transpose(1:Dim)), (outESS["essZZ"]), side=:right, linewidth=0,
    label="ZZ", color=:black, xlabel="Dim", ylabel="ESS", legend=false,
    ylims=(0, maximum(hcat(outESS["essZZ"], outESS["essHMC"]))*1.1))
violin!(string.(transpose(1:Dim)), (outESS["essHMC"]), side=:left, linewidth=0,
    label="HMC", color=:gray)
annotate!([(Dim+0.2, maximum(hcat(outESS["essZZ"], outESS["essHMC"]))*1.05, ("ZZ", 14, :black, :center)),
    (Dim+0.2, maximum(hcat(outESS["essZZ"], outESS["essHMC"]))*0.95, ("HMC", 14, :gray, :center))])
savefig(f11, string(SAVEwd, NameEx, "/f11.pdf"))

# save ESS
JLD.save(string(SAVEwd, NameEx, "/out.jld"),
    Dict("outESS"=>outESS, "zzs"=>zzs, "hmcs"=>hmcs))

# delete everything
f1=f2=f3=f4=f5=f6=f7=f8=f9=f10=f11=nothing
NameEx=io=nothing
tmax_tuned=L_tuned=Lε_tuned=nbZZ_tuned=nbHMC_tuned=start=nothing
zzs=zzsk=zzsm=hmc=hmcs=nothing
Σ=μ=nothing



# -----------------------------------------
# Example on Bivariate Correlated normal
# -----------------------------------------
NameEx = "CorrNorm2d"
mkdir(string(SAVEwd, NameEx))
Dim = 2
μ = zeros(Dim)
Σ = Diagonal(ones(Dim))+zeros( Dim,Dim)
Σ[1,2]=0.90
Σ[2,1]=0.90
function U(x::Vector, m=μ, s=Σ)
    -(-2/2*log(2π)-0.5*log(det(s))-
    0.5*(transpose(x-m)*inv(s)*(x-m))
     )
end

# f1 - plot the target
f1=densplot()
savefig(f1, string(SAVEwd, NameEx, "/f1.pdf"))

# f2 - tune of t_max
f2=tmplot(R=r_tmax, try_tmax=[0.1,0.5,1,1.5,2,2.5])
savefig(f2["NHPP"], string(SAVEwd, NameEx, "/f2a.pdf"))
savefig(f2["opt"], string(SAVEwd, NameEx, "/f2b.pdf"))
savefig(f2["tot"], string(SAVEwd, NameEx, "/f2c.pdf"))
f2tot= plot(f2["NHPP"],f2["opt"],f2["tot"], layout=(1,3),
    size=(pt(252), pt(60)))
savefig(f2tot, string(SAVEwd, NameEx, "/f2.pdf"))
tmax_tuned=1

# f3 - inspection of the zig zag from mode
# multiple chains from the center
f3=plot(0,0, alpha=0, xlabel="X₁", ylabel="X₂", size = (pt(84), pt(60)))
startpoints = Array{Float64, 2}(undef, r_insp, 2)
endpoints = Array{Float64, 2}(undef, r_insp, 2)
for r in 1:r_insp
    start= rand(MvNormal(μ,Σ))
    startpoints[r, :] = start
    zzsk =  zz(; NS=1000, x0_0=start, tmax=tmax_tuned)
    plot!(zzsk["SK"][:, 2], zzsk["SK"][:, 3], color=:gray, linewidth=0.3)
    endpoints[r, :] = zzsk["SK"][1000, 2:3]
end
plot!(startpoints[:,1], startpoints[:,2], seriestype=:scatter,
    markersize=2, color=:black)
plot!(endpoints[:,1], endpoints[:,2], seriestype=:scatter,
        markersize=2, markershape=:star5, color=:black)
savefig(f3, string(SAVEwd, NameEx, "/f3.pdf"))

# f4 - inspection of the zig zag from tail
qsize=0.000001
lim  = quantile(Normal(0,1), qsize)
f4=tailstart(xrange=[-lim, lim], yrange=[-lim, lim], nside=10)
savefig(f4, string(SAVEwd, NameEx, "/f4.pdf"))


# f5 - density inspection
zzsk = zz(; NS=100000, x0_0=rand(MvNormal(μ,Σ)), tmax=tmax_tuned)
zzsm = zzsample(;N=10000, sk=zzsk)
f5=plot(zzsm[:, 1], zzsm[:, 2], seriestype=:scatter,
        color=:black, alpha=0.5, markersize=2,
        xlabel="X₁", ylabel="X₂")
contour!(range(-4, 4, length=100),range(-4,4, length=100),[(pdf(MvNormal(μ, Σ), [x,y])) for y in range(-4, 4, length=100), x in range(-4,4, length=100)],
    fill=(false,cgrad(:grays,[0,0.1,1.0])), color=:lightgray, linewidth=.5)
savefig(f5, string(SAVEwd, NameEx, "/f5.pdf"))

# parameters HMC
Lε_tuned=1.5
L_tuned=3

# f6 - inspection of the HMC from mode
# multiple chains from the center
startpoints = Array{Float64, 2}(undef, r_insp, 2)
endpoints = Array{Float64, 2}(undef, r_insp, 2)
f6=plot(0,0, alpha=0, xlabel="X₁", ylabel="X₂", size = (pt(84), pt(60)))
for r in 1:r_insp
    start= rand(MvNormal(μ,Σ))
    startpoints[r, :] = start
    hmc = runHMC(;epsilon=Lε_tuned/L_tuned,L=L_tuned,IT=1000,qs=start)
    plot!(hmc["SampleQ"][:, 1], hmc["SampleQ"][:, 2], color=:gray, linewidth=0.3)
    endpoints[r, :] = hmc["SampleQ"][1000, 1:2]
end
plot!(startpoints[:,1], startpoints[:,2], seriestype=:scatter,
    markersize=2, color=:black)
plot!(endpoints[:,1], endpoints[:,2], seriestype=:scatter,
        markersize=2, markershape=:star5, color=:black)
savefig(f6, string(SAVEwd, NameEx, "/f6.pdf"))

# f7 - inspection of the HMC from far
f7=tailstart_hmc(xrange=[-lim, lim], yrange=[-lim, lim], nside=10)
savefig(f7, string(SAVEwd, NameEx, "/f7.pdf"))

# f8 - density hmc
hmc = runHMC(;epsilon=Lε_tuned/L_tuned,L=L_tuned,IT=100000,qs=rand(MvNormal(μ,Σ)))
f8=plot(hmc["SampleQ"][1:10:100000, 1], hmc["SampleQ"][1:10:100000, 2], seriestype=:scatter,
        color=:black, alpha=0.5, markersize=2, xlabel="X₁", ylabel="X₂")
contour!(range(-4, 4, length=100),range(-4,4, length=100),[(pdf(MvNormal(μ, Σ), [x,y])) for y in range(-4, 4, length=100), x in range(-4,4, length=100)], fill=(false,cgrad(:grays,[0,0.1,1.0]))
    , color=:lightgray, linewidth=0.5)
savefig(f8, string(SAVEwd, NameEx, "/f8.pdf"))

# -------- PERFORMANCE COMPARISON -------- #
# generate the sample
zzs = Dict()
hmcs= Dict()
for r in 1:r_ess
    start=rand(Dim)
    zzs[[r]]=zz(NS=false, B=budget, x0_0=start, tmax=tmax_tuned)
    hmcs[[r]]=runHMC(;epsilon=Lε_tuned/L_tuned,L=L_tuned,IT=Int(round(budget/(L_tuned+1))),qs=start)
    print(r)
end

# f9 - chose nbatches for zz
f9= bs_zz(;try_nb=[10, 100, 500, 1000])
savefig(f9, string(SAVEwd, NameEx, "/f9.pdf"))
nbZZ_tuned = 100

# f10 - chose nbatches for hmc
f10= bs_hmc(;try_nb=[10, 100, 500, 1000, 2000])
savefig(f10, string(SAVEwd, NameEx, "/f10.pdf"))
nbHMC_tuned = 100

# summaries and f11 - plot of the ESS
outESS = ESSsummaries()
touch(string(SAVEwd, NameEx,"/summaries.txt" ))
io = open(string(SAVEwd, NameEx, "/summaries.txt"), "w")
write(io, outESS["string"])
close(io)
f11=violin(string.(transpose(1:Dim)), (outESS["essZZ"]), side=:right, linewidth=0,
    label="ZZ", color=:black, xlabel="Dim", ylabel="ESS", legend=false,
    ylims=(0, maximum(hcat(outESS["essZZ"], outESS["essHMC"]))*1.1))
violin!(string.(transpose(1:Dim)), (outESS["essHMC"]), side=:left, linewidth=0,
    label="HMC", color=:gray)
annotate!([(Dim+0.2, maximum(hcat(outESS["essZZ"], outESS["essHMC"]))*1.05, ("ZZ", 14, :black, :center)),
    (Dim+0.2, maximum(hcat(outESS["essZZ"], outESS["essHMC"]))*0.95, ("HMC", 14, :gray, :center))])
savefig(f11, string(SAVEwd, NameEx, "/f11.pdf"))

# save ESS
JLD.save(string(SAVEwd, NameEx, "/out.jld"),
    Dict("outESS"=>outESS, "zzs"=>zzs, "hmcs"=>hmcs))

# delete everything
f1=f2=f3=f4=f5=f6=f7=f8=f9=f10=f11=nothing
NameEx=io=nothing
tmax_tuned=L_tuned=Lε_tuned=nbZZ_tuned=nbHMC_tuned=start=nothing
zzs=zzsk=zzsm=hmc=hmcs=nothing
Σ=μ=nothing




# -----------------------------------------
# Example on Normal with different scales
# -----------------------------------------
NameEx = "ScaledNorm2d"
mkdir(string(SAVEwd, NameEx))
Dim = 2
μ = zeros(Dim)
Σ = Diagonal([1.0, 100.0])+zeros( Dim,Dim)
function U(x::Vector, m=μ, s=Σ)
    -(-2/2*log(2π)-0.5*log(det(s))-
    0.5*(transpose(x-m)*inv(s)*(x-m))
     )
end

# f1 - plot the target
f1=densplot()
savefig(f1, string(SAVEwd, NameEx, "/f1.pdf"))

# f2 - tune of t_max
f2=tmplot(R=r_tmax, try_tmax=[3,3.5, 4, 4.5, 5, 6])
savefig(f2["NHPP"], string(SAVEwd, NameEx, "/f2a.pdf"))
savefig(f2["opt"], string(SAVEwd, NameEx, "/f2b.pdf"))
savefig(f2["tot"], string(SAVEwd, NameEx, "/f2c.pdf"))
f2tot= plot(f2["NHPP"],f2["opt"],f2["tot"], layout=(1,3),
    size=(pt(252), pt(60)))
savefig(f2tot, string(SAVEwd, NameEx, "/f2.pdf"))
tmax_tuned=3.5

# f3 - inspection of the zig zag from mode
# multiple chains from the center
f3=plot(0,0, alpha=0, xlabel="X₁", ylabel="X₂", size = (pt(84), pt(60)))
startpoints = Array{Float64, 2}(undef, r_insp, 2)
endpoints = Array{Float64, 2}(undef, r_insp, 2)
for r in 1:r_insp
    start= rand(MvNormal(μ,Σ))
    startpoints[r, :] = start
    zzsk =  zz(; NS=1000, x0_0=start, tmax=tmax_tuned)
    plot!(zzsk["SK"][:, 2], zzsk["SK"][:, 3], color=:gray, linewidth=0.3)
    endpoints[r, :] = zzsk["SK"][1000, 2:3]
end
plot!(startpoints[:,1], startpoints[:,2], seriestype=:scatter,
    markersize=2, color=:black)
plot!(endpoints[:,1], endpoints[:,2], seriestype=:scatter,
        markersize=2, markershape=:star5, color=:black)
savefig(f3, string(SAVEwd, NameEx, "/f3.pdf"))

# f4 - inspection of the zig zag from tail
qsize=0.000001
lim1  = quantile(Normal(0,1), qsize)
lim2  = quantile(Normal(0,10), qsize)
f4=tailstart(xrange=[-lim1, lim1], yrange=[-lim2, lim2], nside=10)
savefig(f4, string(SAVEwd, NameEx, "/f4.pdf"))


# f5 - density inspection
zzsk = zz(; NS=100000, x0_0=rand(MvNormal(μ,Σ)), tmax=tmax_tuned)
zzsm = zzsample(;N=10000, sk=zzsk)
f5=plot(zzsm[:, 1], zzsm[:, 2], seriestype=:scatter,
        color=:black, alpha=0.5, markersize=2,
        xlabel="X₁", ylabel="X₂")
contour!(range(-4, 4, length=100),range(-40,40, length=100),[(pdf(MvNormal(μ, Σ), [x,y])) for y in range(-40, 40, length=100), x in range(-4,4, length=100)],
    fill=(false,cgrad(:grays,[0,0.1,1.0])), color=:lightgray, linewidth=.5)
savefig(f5, string(SAVEwd, NameEx, "/f5.pdf"))

# parameters HMC
Lε_tuned=1.5
L_tuned=2

# f6 - inspection of the HMC from mode
# multiple chains from the center
startpoints = Array{Float64, 2}(undef, r_insp, 2)
endpoints = Array{Float64, 2}(undef, r_insp, 2)

f6=plot(0,0, alpha=0, xlabel="X₁", ylabel="X₂", size = (pt(84), pt(60)))
for r in 1:r_insp
    start= rand(MvNormal(μ,Σ))
    startpoints[r, :] = start
    hmc = runHMC(;epsilon=Lε_tuned/L_tuned,L=L_tuned,IT=1000,qs=start)
    plot!(hmc["SampleQ"][:, 1], hmc["SampleQ"][:, 2], color=:gray, linewidth=0.3)
    endpoints[r, :] = hmc["SampleQ"][1000, 1:2]
end
plot!(startpoints[:,1], startpoints[:,2], seriestype=:scatter,
    markersize=2, color=:black)
plot!(endpoints[:,1], endpoints[:,2], seriestype=:scatter,
        markersize=2, markershape=:star5, color=:black)
savefig(f6, string(SAVEwd, NameEx, "/f6.pdf"))

# f7 - inspection of the HMC from far
f7=tailstart_hmc(xrange=[-lim1, lim1], yrange=[-lim2, lim2], nside=10)
savefig(f7, string(SAVEwd, NameEx, "/f7.pdf"))

# f8 - density hmc
hmc = runHMC(;epsilon=Lε_tuned/L_tuned,L=L_tuned,IT=100000,qs=rand(MvNormal(μ,Σ)))
f8=plot(hmc["SampleQ"][1:10:100000, 1], hmc["SampleQ"][1:10:100000, 2], seriestype=:scatter,
        color=:black, alpha=0.5, markersize=2, xlabel="X₁", ylabel="X₂")
contour!(range(-4, 4, length=100),range(-40,40, length=100),[(pdf(MvNormal(μ, Σ), [x,y])) for y in range(-40, 40, length=100), x in range(-4,4, length=100)], fill=(false,cgrad(:grays,[0,0.1,1.0]))
    , color=:lightgray, linewidth=0.5)
savefig(f8, string(SAVEwd, NameEx, "/f8.pdf"))

# -------- PERFORMANCE COMPARISON -------- #
# generate the sample
zzs = Dict()
hmcs= Dict()
for r in 1:r_ess
    start=rand(Dim)
    zzs[[r]]=zz(NS=false, B=budget, x0_0=start, tmax=tmax_tuned)
    hmcs[[r]]=runHMC(;epsilon=Lε_tuned/L_tuned,L=L_tuned,IT=Int(round(budget/(L_tuned+1))),qs=start)
    print(r)
end

# f9 - chose nbatches for zz
f9= bs_zz(;try_nb=[10, 100, 500, 1000, 2000])
savefig(f9, string(SAVEwd, NameEx, "/f9.pdf"))
nbZZ_tuned = 100

# f10 - chose nbatches for hmc
f10= bs_hmc(;try_nb=[ 5, 10, 100, 500, 1000, 2000])
savefig(f10, string(SAVEwd, NameEx, "/f10.pdf"))
nbHMC_tuned = 10

# summaries and f11 - plot of the ESS
outESS = ESSsummaries()
touch(string(SAVEwd, NameEx,"/summaries.txt" ))
io = open(string(SAVEwd, NameEx, "/summaries.txt"), "w")
write(io, outESS["string"])
close(io)
f11=violin(string.(transpose(1:Dim)), (outESS["essZZ"]), side=:right, linewidth=0,
    label="ZZ", color=:black, xlabel="Dim", ylabel="ESS", legend=false,
    ylims=(0, maximum(hcat(outESS["essZZ"], outESS["essHMC"]))*1.1))
violin!(string.(transpose(1:Dim)), (outESS["essHMC"]), side=:left, linewidth=0,
    label="HMC", color=:gray)
annotate!([(Dim+0.2, maximum(hcat(outESS["essZZ"], outESS["essHMC"]))*1.05, ("ZZ", 14, :black, :center)),
    (Dim+0.2, maximum(hcat(outESS["essZZ"], outESS["essHMC"]))*0.95, ("HMC", 14, :gray, :center))])
savefig(f11, string(SAVEwd, NameEx, "/f11.pdf"))

# save ESS
JLD.save(string(SAVEwd, NameEx, "/out.jld"),
    Dict("outESS"=>outESS, "zzs"=>zzs, "hmcs"=>hmcs))

# delete everything
f1=f2=f3=f4=f5=f6=f7=f8=f9=f10=f11=nothing
NameEx=io=nothing
tmax_tuned=L_tuned=Lε_tuned=nbZZ_tuned=nbHMC_tuned=start=nothing
zzs=zzsk=zzsm=hmc=hmcs=nothing
Σ=μ=nothing




# -----------------------------------------
# Example on Bimodal  normal
# -----------------------------------------
NameEx = "BimodNorm2d"
mkdir(string(SAVEwd, NameEx))
Dim = 2
μ₁ = zeros(Dim).-2
μ₂ = zeros(Dim).+2
Σ₁ = Diagonal(ones(Dim)*2)+zeros( Dim,Dim)
Σ₂ = Diagonal(ones(Dim)*2)+zeros( Dim,Dim)
θ = [0.5, 0.5]

function g(x::Vector, m1=μ₁, m2=μ₂, s1=Σ₁, s2=Σ₂, p=θ)
    gx = p[1]*(det(2π.*s1))^(-.5) * exp(-0.5*(transpose(x.-m1)*inv(s1)*(x-m1)))+
        p[2]*(det(2π.*s2))^(-.5) * exp(-0.5*(transpose(x.-m2)*inv(s2)*(x-m2)))
    return gx
end
 g([-4,-4])

function U(x::Vector, m1=μ₁, m2=μ₂, s1=Σ₁, s2=Σ₂, p=θ)
    gx = p[1]*(det(2π.*s1))^(-.5) * exp(-0.5*(transpose(x.-m1)*inv(s1)*(x-m1)))+
        p[2]*(det(2π.*s2))^(-.5) * exp(-0.5*(transpose(x.-m2)*inv(s2)*(x-m2)))
    return -log(gx)
end

f1=densplot(xrange=[-7,7], yrange=[-7,7])
savefig(f1, string(SAVEwd, NameEx, "/f1.pdf"))

# f2 - tune of t_max
# note tmax should be tested starting from the difference between the modes

f2=tmplot(R=r_tmax, try_tmax=[2, 3, 4, 5, 6])
savefig(f2["NHPP"], string(SAVEwd, NameEx, "/f2a.pdf"))
savefig(f2["opt"], string(SAVEwd, NameEx, "/f2b.pdf"))
savefig(f2["tot"], string(SAVEwd, NameEx, "/f2c.pdf"))
f2tot= plot(f2["NHPP"],f2["opt"],f2["tot"], layout=(1,3),
    size=(pt(252), pt(60)))
savefig(f2tot, string(SAVEwd, NameEx, "/f2.pdf"))
tmax_tuned=4

# f3 - inspection of the zig zag from mode
# multiple chains from the center
f3=plot(0,0, alpha=0, xlabel="X₁", ylabel="X₂", size = (pt(84), pt(60)))
startpoints = Array{Float64, 2}(undef, r_insp, 2)
endpoints = Array{Float64, 2}(undef, r_insp, 2)
for r in 1:r_insp
    mixt = sample([1,2])
    start= ifelse(mixt==1, rand(MvNormal(μ₁,Σ₁)), rand(MvNormal(μ₂,Σ₂)))
    startpoints[r, :] = start
    zzsk =  zz(; NS=1000, x0_0=start, tmax=tmax_tuned)
    plot!(zzsk["SK"][:, 2], zzsk["SK"][:, 3], color=:gray, linewidth=0.3)
    endpoints[r, :] = zzsk["SK"][1000, 2:3]
end
plot!(startpoints[:,1], startpoints[:,2], seriestype=:scatter,
    markersize=2, color=:black)
plot!(endpoints[:,1], endpoints[:,2], seriestype=:scatter,
        markersize=2, markershape=:star5, color=:black)
savefig(f3, string(SAVEwd, NameEx, "/f3.pdf"))

# f4 - inspection of the zig zag from tail
qsize=0.000001
lim  = quantile(Normal(-2,2), qsize*2)
f4=tailstart(xrange=[-lim, lim], yrange=[-lim, lim], nside=10)
savefig(f4, string(SAVEwd, NameEx, "/f4.pdf"))



# f5 - density inspection
zzsk = zz(; NS=100000, x0_0=ifelse(rand([1,2])==1,
        rand(MvNormal(μ₁,Σ₁)), rand(MvNormal(μ₂,Σ₂))), tmax=tmax_tuned)
zzsm = zzsample(;N=10000, sk=zzsk)
f5=plot(zzsm[:, 1], zzsm[:, 2], seriestype=:scatter,
        color=:black, alpha=0.5, markersize=2,
        xlabel="X₁", ylabel="X₂")
contour!(range(-6, 6, length=100),range(-6,6, length=100),
    [(exp(-U([x,y]))) for y in range(-6, 6, length=100),
        x in range(-6,6, length=100)],
        fill=(false,cgrad(:grays,[0,0.1,1.0])), color=:lightgray,linewidth=.5)
savefig(f5, string(SAVEwd, NameEx, "/f5.pdf"))

# parameters HMC
Lε_tuned=3
L_tuned=2
hmc = runHMC(;epsilon=Lε_tuned/L_tuned,L=L_tuned,IT=1000,qs=start)
plot(plot(hmc["SampleQ"][:, 1]), plot(hmc["SampleQ"][:, 2]),
        layout=(2,1), size=(1000,800))
ac=plot(autocor(hmc["SampleQ"][:, 1]))
plot!(autocor(hmc["SampleQ"][:, 2]))
sum(hmc["accept"])/1000


# f6 - inspection of the HMC from mode
# multiple chains from the center
startpoints = Array{Float64, 2}(undef, r_insp, 2)
endpoints = Array{Float64, 2}(undef, r_insp, 2)
f6=plot(0,0, alpha=0, xlabel="X₁", ylabel="X₂", size = (pt(84), pt(60)))
for r in 1:r_insp
    mixt = sample([1,2])
    start= ifelse(mixt==1, rand(MvNormal(μ₁,Σ₁)), rand(MvNormal(μ₂,Σ₂)))
    startpoints[r, :] = start
    hmc = runHMC(;epsilon=Lε_tuned/L_tuned,L=L_tuned,IT=1000,qs=start)
    plot!(hmc["SampleQ"][:, 1], hmc["SampleQ"][:, 2], color=:gray, linewidth=0.3)
    endpoints[r, :] = hmc["SampleQ"][1000, 1:2]
end
plot!(startpoints[:,1], startpoints[:,2], seriestype=:scatter,
    markersize=2, color=:black)
plot!(endpoints[:,1], endpoints[:,2], seriestype=:scatter,
        markersize=2, markershape=:star5, color=:black)
savefig(f6, string(SAVEwd, NameEx, "/f6.pdf"))

# f7 - inspection of the HMC from far
f7=tailstart_hmc(xrange=[-lim, lim], yrange=[-lim, lim], nside=10)
savefig(f7, string(SAVEwd, NameEx, "/f7.pdf"))

# f8 - density hmc
hmc = runHMC(;epsilon=Lε_tuned/L_tuned,L=L_tuned,IT=100000,qs=ifelse(rand([1,2])==1,
        rand(MvNormal(μ₁,Σ₁)), rand(MvNormal(μ₂,Σ₂))))
f8=plot(hmc["SampleQ"][1:10:100000, 1], hmc["SampleQ"][1:10:100000, 2], seriestype=:scatter,
        color=:black, alpha=0.5, markersize=2, xlabel="X₁", ylabel="X₂")
contour!(range(-6, 6, length=100),range(-6,6, length=100),
    [(exp(-U([x,y]))) for y in range(-6, 6, length=100),
        x in range(-6,6, length=100)],
        fill=(false,cgrad(:grays,[0,0.1,1.0]))
            , color=:lightgray, linewidth=0.5)
savefig(f8, string(SAVEwd, NameEx, "/f8.pdf"))

# -------- PERFORMANCE COMPARISON -------- #
# generate the sample
zzs = Dict()
hmcs= Dict()
for r in 1:r_ess
    start=rand(Dim)
    zzs[[r]]=zz(NS=false, B=budget, x0_0=start, tmax=tmax_tuned)
    hmcs[[r]]=runHMC(;epsilon=Lε_tuned/L_tuned,L=L_tuned,IT=Int(round(budget/(L_tuned+1))),qs=start)
    print(r)
end

# f9 - chose nbatches for zz
f9= bs_zz(;try_nb=[3, 5, 10, 50, 100, 200, 500])
savefig(f9, string(SAVEwd, NameEx, "/f9.pdf"))
nbZZ_tuned = 10

# f10 - chose nbatches for hmc
f10= bs_hmc(;try_nb=[5, 10, 100, 500, 1000, 2000])
savefig(f10, string(SAVEwd, NameEx, "/f10.pdf"))
nbHMC_tuned = 100

# summaries and f11 - plot of the ESS
outESS = ESSsummaries()
outESS


touch(string(SAVEwd, NameEx,"/summaries.txt" ))
io = open(string(SAVEwd, NameEx, "/summaries.txt"), "w")
write(io, outESS["string"])
close(io)
f11=violin(string.(transpose(1:Dim)), (outESS["essZZ"]), side=:right, linewidth=0,
    label="ZZ", color=:black, xlabel="Dim", ylabel="ESS", legend=false,
    ylims=(0, maximum(hcat(outESS["essZZ"], outESS["essHMC"]))*1.1))
violin!(string.(transpose(1:Dim)), (outESS["essHMC"]), side=:left, linewidth=0,
    label="HMC", color=:gray)
annotate!([(Dim+0.2, maximum(hcat(outESS["essZZ"], outESS["essHMC"]))*1.05, ("ZZ", 14, :black, :center)),
    (Dim+0.2, maximum(hcat(outESS["essZZ"], outESS["essHMC"]))*0.95, ("HMC", 14, :gray, :center))])
savefig(f11, string(SAVEwd, NameEx, "/f11.pdf"))

# save ESS
JLD.save(string(SAVEwd, NameEx, "/out.jld"),
    Dict("outESS"=>outESS, "zzs"=>zzs, "hmcs"=>hmcs))

# delete everything
f1=f2=f3=f4=f5=f6=f7=f8=f9=f10=f11=nothing
NameEx=io=nothing
tmax_tuned=L_tuned=Lε_tuned=nbZZ_tuned=nbHMC_tuned=start=nothing
zzs=zzsk=zzsm=hmc=hmcs=nothing
Σ=μ=nothing





# -----------------------------------------
# Example on light tailed doistributions
# -----------------------------------------
NameEx = "LightTail2d"
mkdir(string(SAVEwd, NameEx))
Dim = 2
function gx(x)
    exp(-(sum(x.^4))/4)/((gamma(0.25)^2)/(2))
end
function U(x::Vector)
    sum((x).^4)/4
end


# f1 - plot the target
f1=densplot()
savefig(f1, string(SAVEwd, NameEx, "/f1.pdf"))

# f2 - tune of t_max
f2=tmplot(R=r_tmax, try_tmax=[0.5, 1, 1.5, 2,2.5, 3])
savefig(f2["NHPP"], string(SAVEwd, NameEx, "/f2a.pdf"))
savefig(f2["opt"], string(SAVEwd, NameEx, "/f2b.pdf"))
savefig(f2["tot"], string(SAVEwd, NameEx, "/f2c.pdf"))
f2tot= plot(f2["NHPP"],f2["opt"],f2["tot"], layout=(1,3),
    size=(pt(252), pt(60)))
savefig(f2tot, string(SAVEwd, NameEx, "/f2.pdf"))
tmax_tuned=1.5

# f3 - inspection of the zig zag from mode
# multiple chains from the center
f3=plot(0,0, alpha=0, xlabel="X₁", ylabel="X₂", size = (pt(84), pt(60)))
startpoints = Array{Float64, 2}(undef, r_insp, 2)
endpoints = Array{Float64, 2}(undef, r_insp, 2)
for r in 1:r_insp
    start= rand(MvNormal([0,0],Diagonal(ones(Dim)).*0.1+zeros( Dim,Dim)))
    startpoints[r, :] = start
    zzsk =  zz(; NS=1000, x0_0=start, tmax=tmax_tuned)
    plot!(zzsk["SK"][:, 2], zzsk["SK"][:, 3], color=:gray, linewidth=0.3)
    endpoints[r, :] = zzsk["SK"][1000, 2:3]
end
plot!(startpoints[:,1], startpoints[:,2], seriestype=:scatter,
    markersize=2, color=:black)
plot!(endpoints[:,1], endpoints[:,2], seriestype=:scatter,
        markersize=2, markershape=:star5, color=:black)
savefig(f3, string(SAVEwd, NameEx, "/f3.pdf"))

# f4 - inspection of the zig zag from tail
qsize=0.000001
lim  = quantile(Normal(0,1), qsize)
f4=tailstart(xrange=[-lim, lim], yrange=[-lim, lim], nside=10)
savefig(f4, string(SAVEwd, NameEx, "/f4.pdf"))


# f5 - density inspection
zzsk = zz(; NS=100000, x0_0=rand(MvNormal([0,0],Diagonal(ones(Dim)).*0.1+zeros( Dim,Dim))), tmax=tmax_tuned)
zzsm = zzsample(;N=10000, sk=zzsk)
f5=plot(zzsm[:, 1], zzsm[:, 2], seriestype=:scatter,
        color=:black, alpha=0.5, markersize=2,
        xlabel="X₁", ylabel="X₂")
contour!(range(-4, 4, length=100),range(-4,4, length=100),[(gx([x,y])) for y in range(-4, 4, length=100), x in range(-4,4, length=100)], fill=(false,cgrad(:grays,[0,0.1,1.0]))
    , color=:lightgray,linewidth=.5)
savefig(f5, string(SAVEwd, NameEx, "/f5.pdf"))

# parameters HMC
Lε_tuned=1
L_tuned=2

# f6 - inspection of the HMC from mode
# multiple chains from the center
startpoints = Array{Float64, 2}(undef, r_insp, 2)
endpoints = Array{Float64, 2}(undef, r_insp, 2)
f6=plot(0,0, alpha=0, xlabel="X₁", ylabel="X₂", size = (pt(84), pt(60)))
for r in 1:r_insp
    start= rand(MvNormal([0,0],Diagonal(ones(Dim)).*0.1+zeros( Dim,Dim)))
    startpoints[r, :] = start
    hmc = runHMC(;epsilon=Lε_tuned/L_tuned,L=L_tuned,IT=1000,qs=start)
    plot!(hmc["SampleQ"][:, 1], hmc["SampleQ"][:, 2], color=:gray, linewidth=0.3)
    endpoints[r, :] = hmc["SampleQ"][1000, 1:2]
end
plot!(startpoints[:,1], startpoints[:,2], seriestype=:scatter,
    markersize=2, color=:black)
plot!(endpoints[:,1], endpoints[:,2], seriestype=:scatter,
        markersize=2, markershape=:star5, color=:black)
savefig(f6, string(SAVEwd, NameEx, "/f6.pdf"))

# f7 - inspection of the HMC from far
f7=tailstart_hmc(xrange=[-lim, lim], yrange=[-lim, lim], nside=10)
savefig(f7, string(SAVEwd, NameEx, "/f7.pdf"))

# f8 - density hmc
hmc = runHMC(;epsilon=Lε_tuned/L_tuned,L=L_tuned,IT=100000,qs=rand(MvNormal([0,0],Diagonal(ones(Dim)).*0.1+zeros( Dim,Dim))))
f8=plot(hmc["SampleQ"][1:10:100000, 1], hmc["SampleQ"][1:10:100000, 2], seriestype=:scatter,
        color=:black, alpha=0.5, markersize=2, xlabel="X₁", ylabel="X₂")
contour!(range(-4, 4, length=100),range(-4,4, length=100),[(gx([x,y])) for y in range(-4, 4, length=100), x in range(-4,4, length=100)], fill=(false,cgrad(:grays,[0,0.1,1.0]))
    , color=:lightgray, linewidth=0.5)
savefig(f8, string(SAVEwd, NameEx, "/f8.pdf"))

# -------- PERFORMANCE COMPARISON -------- #
# generate the sample
zzs = Dict()
hmcs= Dict()
for r in 1:r_ess
    start=rand(Dim)
    zzs[[r]]=zz(NS=false, B=budget, x0_0=start, tmax=tmax_tuned)
    hmcs[[r]]=runHMC(;epsilon=Lε_tuned/L_tuned,L=L_tuned,IT=Int(round(budget/(L_tuned+1))),qs=start)
    print(r)
end

# f9 - chose nbatches for zz
f9= bs_zz(;try_nb=[10, 100, 500, 1000, 2000])
savefig(f9, string(SAVEwd, NameEx, "/f9.pdf"))
nbZZ_tuned = 1000

# f10 - chose nbatches for hmc
f10= bs_hmc(;try_nb=[10, 100, 500, 1000, 2000])
savefig(f10, string(SAVEwd, NameEx, "/f10.pdf"))
nbHMC_tuned = 100

# summaries and f11 - plot of the ESS
outESS = ESSsummaries()
touch(string(SAVEwd, NameEx,"/summaries.txt" ))
io = open(string(SAVEwd, NameEx, "/summaries.txt"), "w")
write(io, outESS["string"])
close(io)
f11=violin(string.(transpose(1:Dim)), (outESS["essZZ"]), side=:right, linewidth=0,
    label="ZZ", color=:black,  xlabel="Dim", ylabel="ESS", legend=false,
    ylims=(0, maximum(hcat(outESS["essZZ"], outESS["essHMC"]))*1.1))
violin!(string.(transpose(1:Dim)), (outESS["essHMC"]), side=:left, linewidth=0,
    label="HMC", color=:gray)
annotate!([(Dim+0.2, maximum(hcat(outESS["essZZ"], outESS["essHMC"]))*1.05, ("ZZ", 14, :black, :center)),
    (Dim+0.2, maximum(hcat(outESS["essZZ"], outESS["essHMC"]))*0.95, ("HMC", 14, :gray, :center))])
savefig(f11, string(SAVEwd, NameEx, "/f11.pdf"))

# save ESS
JLD.save(string(SAVEwd, NameEx, "/out.jld"),
    Dict("outESS"=>outESS, "zzs"=>zzs, "hmcs"=>hmcs))

# delete everything
f1=f2=f3=f4=f5=f6=f7=f8=f9=f10=f11=nothing
NameEx=io=nothing
tmax_tuned=L_tuned=Lε_tuned=nbZZ_tuned=nbHMC_tuned=start=nothing
zzs=zzsk=zzsm=hmc=hmcs=nothing
Σ=μ=nothing



# -----------------------------------------
# Example on Heavy Tailed normal
# -----------------------------------------
NameEx = "HeavyTail2d"
mkdir(string(SAVEwd, NameEx))
Dim=2
μ = zeros(Dim)
Σ = Diagonal(ones(Dim))+zeros( Dim,Dim)
# heavy-tail distribution
function gx(x; ν=1, m=μ, s=Σ)
    p=size(x)[1]
    (gamma((ν+p)/2)/(gamma(ν/2)*((ν*π)^(p/2))*det(s)))*((1+((1/ν)*(transpose(x-m)*inv(s)*(x-m))))^(-(ν+p)/2))
end
function U(x::Vector; ν=1, m=μ, s=Σ)
    p=size(x)[1]
    -loggamma((ν+p)/2)+loggamma(ν/2)+(p/2)*log(ν*π)+1/2*log(det(Σ))+((ν+p)/2)*log(1+((1/ν)*(transpose(x-m)*inv(s)*(x-m))))
end


# f1 - plot the target
f1=densplot()
savefig(f1, string(SAVEwd, NameEx, "/f1.pdf"))

# f2 - tune of t_max
f2=tmplot(R=r_tmax, try_tmax=[2, 3.5, 5, 6.5, 8, 9.5])
savefig(f2["NHPP"], string(SAVEwd, NameEx, "/f2a.pdf"))
savefig(f2["opt"], string(SAVEwd, NameEx, "/f2b.pdf"))
savefig(f2["tot"], string(SAVEwd, NameEx, "/f2c.pdf"))
f2tot= plot(f2["NHPP"],f2["opt"],f2["tot"], layout=(1,3),
    size=(pt(252), pt(60)))
savefig(f2tot, string(SAVEwd, NameEx, "/f2.pdf"))
tmax_tuned=5

# f3 - inspection of the zig zag from mode
# multiple chains from the center
f3=plot(0,0, alpha=0, xlabel="X₁", ylabel="X₂", size = (pt(84), pt(60)))
startpoints = Array{Float64, 2}(undef, r_insp, 2)
endpoints = Array{Float64, 2}(undef, r_insp, 2)
for r in 1:r_insp
    start= rand(MvNormal([0,0],Diagonal(ones(Dim)).*2+zeros( Dim,Dim)))
    startpoints[r, :] = start
    zzsk =  zz(; NS=1000, x0_0=start, tmax=tmax_tuned)
    plot!(zzsk["SK"][:, 2], zzsk["SK"][:, 3], color=:gray, linewidth=0.3)
    endpoints[r, :] = zzsk["SK"][1000, 2:3]
end
plot!(startpoints[:,1], startpoints[:,2], seriestype=:scatter,
    markersize=2, color=:black)
plot!(endpoints[:,1], endpoints[:,2], seriestype=:scatter,
        markersize=2, markershape=:star5, color=:black)
savefig(f3, string(SAVEwd, NameEx, "/f3.pdf"))

# f4 - inspection of the zig zag from tail
qsize=0.000001
lim  = quantile(TDist(2), qsize)
f4=tailstart(xrange=[-lim, lim], yrange=[-lim, lim], nside=10)
savefig(f4, string(SAVEwd, NameEx, "/f4.pdf"))


# f5 - density inspection
zzsk = zz(; NS=100000, x0_0=rand(MvNormal([0,0],Diagonal(ones(Dim)).*2+zeros( Dim,Dim))), tmax=tmax_tuned)
zzsm = zzsample(;N=10000, sk=zzsk)
f5=plot(zzsm[:, 1], zzsm[:, 2], seriestype=:scatter,
        color=:black, alpha=0.5, markersize=2,
        xlabel="X₁", ylabel="X₂")
contour!(range(-50, 50, length=100),range(-50,50, length=100),[(gx([x,y])) for y in range(-50, 50, length=100), x in range(-50,50, length=100)], fill=(false,cgrad(:grays,[0,0.1,1.0]))
    , color=:lightgray,linewidth=.5)
savefig(f5, string(SAVEwd, NameEx, "/f5.pdf"))

# parameters HMC
Lε_tuned=10
L_tuned=10

# f6 - inspection of the HMC from mode
# multiple chains from the center
startpoints = Array{Float64, 2}(undef, r_insp, 2)
endpoints = Array{Float64, 2}(undef, r_insp, 2)
f6=plot(0,0, alpha=0, xlabel="X₁", ylabel="X₂", size = (pt(84), pt(60)))
for r in 1:r_insp
    start= rand(MvNormal([0,0],Diagonal(ones(Dim)).*2+zeros( Dim,Dim)))
    startpoints[r, :] = start
    hmc = runHMC(;epsilon=Lε_tuned/L_tuned,L=L_tuned,IT=1000,qs=start)
    plot!(hmc["SampleQ"][:, 1], hmc["SampleQ"][:, 2], color=:gray, linewidth=0.3)
    endpoints[r, :] = hmc["SampleQ"][1000, 1:2]
end
plot!(startpoints[:,1], startpoints[:,2], seriestype=:scatter,
    markersize=2, color=:black)
plot!(endpoints[:,1], endpoints[:,2], seriestype=:scatter,
        markersize=2, markershape=:star5, color=:black)
savefig(f6, string(SAVEwd, NameEx, "/f6.pdf"))

# f7 - inspection of the HMC from far
f7=tailstart_hmc(xrange=[-lim, lim], yrange=[-lim, lim], nside=10)
savefig(f7, string(SAVEwd, NameEx, "/f7.pdf"))

# f8 - density hmc
hmc = runHMC(;epsilon=Lε_tuned/L_tuned,L=L_tuned,IT=100000,qs=rand(MvNormal([0,0],Diagonal(ones(Dim)).*2+zeros( Dim,Dim))))
f8=plot(hmc["SampleQ"][1:10:100000, 1], hmc["SampleQ"][1:10:100000, 2], seriestype=:scatter,
        color=:black, alpha=0.5, markersize=2, xlabel="X₁", ylabel="X₂")
contour!(range(-4, 4, length=100),range(-4,4, length=100),[(gx([x,y])) for y in range(-4, 4, length=100), x in range(-4,4, length=100)], fill=(false,cgrad(:grays,[0,0.1,1.0]))
    , color=:lightgray, linewidth=0.5)
savefig(f8, string(SAVEwd, NameEx, "/f8.pdf"))

# -------- PERFORMANCE COMPARISON -------- #
# generate the sample
zzs = Dict()
hmcs= Dict()
for r in 1:r_ess
    start=rand(Dim)
    zzs[[r]]=zz(NS=false, B=budget, x0_0=start, tmax=tmax_tuned)
    hmcs[[r]]=runHMC(;epsilon=Lε_tuned/L_tuned,L=L_tuned,IT=Int(round(budget/(L_tuned+1))),qs=start)
    print(r)
end

# f9 - chose nbatches for zz
f9= bs_tp_zz(;try_nb=[5, 10, 100, 500, 1000])
savefig(f9, string(SAVEwd, NameEx, "/f9.pdf"))
nbZZ_tuned = 10

# f10 - chose nbatches for hmc
f10= bs_hmc(;try_nb=[5, 10, 100, 500, 1000])
savefig(f10, string(SAVEwd, NameEx, "/f10.pdf"))
nbHMC_tuned = 10

# summaries and f11 - plot of the ESS
outESS = ESSsummaries()
touch(string(SAVEwd, NameEx,"/summaries.txt" ))
io = open(string(SAVEwd, NameEx, "/summaries.txt"), "w")
write(io, outESS["string"])
close(io)
f11=violin(string.(transpose(1:Dim)), (outESS["essZZ"]), side=:right, linewidth=0,
    label="ZZ", color=:black, xlabel="Dim", ylabel="ESS", legend=false,
    ylims=(0, maximum(hcat(outESS["essZZ"], outESS["essHMC"]))*1.1))
violin!(string.(transpose(1:Dim)), (outESS["essHMC"]), side=:left, linewidth=0,
    label="HMC", color=:gray)
annotate!([(Dim+0.2, maximum(hcat(outESS["essZZ"], outESS["essHMC"]))*1.05, ("ZZ", 14, :black, :center)),
    (Dim+0.2, maximum(hcat(outESS["essZZ"], outESS["essHMC"]))*0.95, ("HMC", 14, :gray, :center))])
savefig(f11, string(SAVEwd, NameEx, "/f11.pdf"))

# save ESS
JLD.save(string(SAVEwd, NameEx, "/out.jld"),
    Dict("outESS"=>outESS, "zzs"=>zzs, "hmcs"=>hmcs))

# delete everything
f1=f2=f3=f4=f5=f6=f7=f8=f9=f10=f11=nothing
NameEx=io=nothing
tmax_tuned=L_tuned=Lε_tuned=nbZZ_tuned=nbHMC_tuned=start=nothing
zzs=zzsk=zzsm=hmc=hmcs=nothing
Σ=μ=nothing

# HERE STILL TO DO!!!

# -----------------------------------------
# Example on 10variate Isotropic normal
# -----------------------------------------
NameEx = "IsoNorm10d"
mkdir(string(SAVEwd, NameEx))
Dim = 10
μ = zeros(Dim)
Σ = Diagonal(ones(Dim))+zeros( Dim,Dim)
function U(x::Vector, m=μ, s=Σ)
    -(-2/2*log(2π)-0.5*log(det(s))-
    0.5*(transpose(x-m)*inv(s)*(x-m))
     )
end

# f2 - tune of t_max
f2=tmplot(R=r_tmax, try_tmax=[0.2, 0.4, 0.6, 0.8, 1])
savefig(f2, string(SAVEwd, NameEx, "/f2.pdf"))
tmax_tuned=0.6

# f3 - inspection of the zig zag from mode
# multiple chains from the center
zzsk_Dict=Dict()
for r in 1:r_insp
    start=zeros(Dim)
    zzsk_Dict[r] =  zz(; NS=1000, x0_0=start, tmax=tmax_tuned)
end

for d in 1:Dim
    fx = plot(zzsk_Dict[1]["SK"][:, 1], zzsk_Dict[1]["SK"][:, d+1],
        xlabel="Time", ylabel=string("X",d))
    for r in 1:r_insp
        plot!(zzsk_Dict[r]["SK"][:, 1], zzsk_Dict[r]["SK"][:, d+1])
    end
    savefig(fx, string(SAVEwd, NameEx, "/f3d", d,".pdf"))
end


# f4 - inspection of the zig zag from tail

eachdim_range = [-7,7]

# multiple chains from some sampled extremes
zzsk_Dict=Dict()
for r in 1:r_insp
    start = sample(eachdim_range, 10)
    zzsk_Dict[r] =  zz(; NS=1000, x0_0=start, tmax=tmax_tuned)
end

for d in 1:Dim
    fx = plot(zzsk_Dict[1]["SK"][:, 1], zzsk_Dict[1]["SK"][:, d+1],
        xlabel="Time", ylabel=string("X",d))
    for r in 1:r_insp
        plot!(zzsk_Dict[r]["SK"][:, 1], zzsk_Dict[r]["SK"][:, d+1])
    end
    savefig(fx, string(SAVEwd, NameEx, "/f4d", d,".pdf"))
end

# CONTINUE FROM HERE


# f5 - density inspection
zzsk = zz(; NS=100000, x0_0=start, tmax=tmax_tuned)
zzsm = zzsample(;N=10000, sk=zzsk)
f5=plot(zzsm[:, 1], zzsm[:, 2], seriestype=:scatter,
        color=:black, alpha=0.5, markersize=3,
        xlabel="X₁", ylabel="X₂")
contour!(range(-4, 4, length=100),range(-4,4, length=100),
    [(pdf(MvNormal([0,0], [[1.0,0]  [0,1]]), [x,y])) for y in range(-4, 4, length=100), x in range(-4,4, length=100)],
        fill=(false,cgrad(:grays,[0,0.1,1.0])))
savefig(f5, string(SAVEwd, NameEx, "/f5a.pdf"))
f5=plot(zzsm[:, 9], zzsm[:, 10], seriestype=:scatter,
        color=:black, alpha=0.5, markersize=3,
        xlabel="X₉", ylabel="X₁₀")
contour!(range(-4, 4, length=100),range(-4,4, length=100),
    [(pdf(MvNormal([0,0], [[1.0,0]  [0,1]]), [x,y])) for y in range(-4, 4, length=100), x in range(-4,4, length=100)],
        fill=(false,cgrad(:grays,[0,0.1,1.0])))
savefig(f5, string(SAVEwd, NameEx, "/f5b.pdf"))


# parameters HMC
Lε_tuned=1
L_tuned=2

hmc = runHMC(;epsilon=Lε_tuned/L_tuned,L=L_tuned,IT=1000,qs=start)


plot(plot(hmc["SampleQ"][:, 1]), plot(hmc["SampleQ"][:, 2]),
    plot(hmc["SampleQ"][:, 3]), plot(hmc["SampleQ"][:, 4]),
    plot(hmc["SampleQ"][:, 5]), plot(hmc["SampleQ"][:, 6]),
    plot(hmc["SampleQ"][:, 7]), plot(hmc["SampleQ"][:, 8]),
    plot(hmc["SampleQ"][:, 9]), plot(hmc["SampleQ"][:, 10]),
    layout=(5,2), size=(1000,800))


ac=plot(autocor(hmc["SampleQ"][:, 1]))
for d in 1:Dim
    plot!(autocor(hmc["SampleQ"][:, d]))
end
display(ac)



sum(hmc["accept"])/1000


# f6 - inspection of the HMC from mode
# multiple chains from the center
hmc_Dict=Dict()
for r in 1:r_insp
    start=zeros(Dim)
    hmc_Dict[r] =  runHMC(;epsilon=Lε_tuned/L_tuned,L=L_tuned,IT=1000,qs=start)
end

for d in 1:Dim
    fx = plot(hmc_Dict[1]["SampleQ"][:, 1],
        xlabel="Iterations", ylabel=string("X",d))
    for r in 1:r_insp
        plot!(hmc_Dict[r]["SampleQ"][:, d])
    end
    savefig(fx, string(SAVEwd, NameEx, "/f6d", d,".pdf"))
end

# f7 - inspection of the HMC from far
hmc_Dict=Dict()
for r in 1:r_insp
    start = sample(eachdim_range, 10)
    hmc_Dict[r] =  runHMC(;epsilon=Lε_tuned/L_tuned,L=L_tuned,IT=1000,qs=start)
end
for d in 1:Dim
    fx = plot(hmc_Dict[1]["SampleQ"][:, 1],
        xlabel="Iterations", ylabel=string("X",d))
    for r in 1:r_insp
        plot!(hmc_Dict[r]["SampleQ"][:, d])
    end
    savefig(fx, string(SAVEwd, NameEx, "/f7d", d,".pdf"))
end



# f8 - density hmc
hmc = runHMC(;epsilon=Lε_tuned/L_tuned,L=L_tuned,IT=100000,qs=start)
f8=plot(hmc["SampleQ"][1:10:100000, 1], hmc["SampleQ"][1:10:100000, 2], seriestype=:scatter,
        color=:black, alpha=0.5, markersize=3, xlabel="X₁", ylabel="X₂")
contour!(range(-4, 4, length=100),range(-4,4, length=100),
    [(pdf(MvNormal([0,0], [[1.0,0]  [0,1]]), [x,y])) for y in range(-4, 4, length=100), x in range(-4,4, length=100)],
        fill=(false,cgrad(:grays,[0,0.1,1.0])))
savefig(f8, string(SAVEwd, NameEx, "/f8a.pdf"))
f8=plot(hmc["SampleQ"][1:10:100000, 9], hmc["SampleQ"][1:10:100000, 10], seriestype=:scatter,
        color=:black, alpha=0.5, markersize=3, xlabel="X₉", ylabel="X₁₀")
contour!(range(-4, 4, length=100),range(-4,4, length=100),
    [(pdf(MvNormal([0,0], [[1.0,0]  [0,1]]), [x,y])) for y in range(-4, 4, length=100), x in range(-4,4, length=100)],
        fill=(false,cgrad(:grays,[0,0.1,1.0])))
savefig(f8, string(SAVEwd, NameEx, "/f8b.pdf"))

# -------- PERFORMANCE COMPARISON -------- #
# generate the sample
zzs = Dict()
hmcs= Dict()
for r in 1:r_ess
    start=rand(Dim)
    zzs[[r]]=zz(NS=false, B=budget, x0_0=start, tmax=tmax_tuned)
    hmcs[[r]]=runHMC(;epsilon=Lε_tuned/L_tuned,L=L_tuned,IT=Int(round(budget/(L_tuned+1))),qs=start)
    print(r)
end

# f9 - chose nbatches for zz
f9= bs_zz(;try_nb=[100,200,300,400,500])
savefig(f9, string(SAVEwd, NameEx, "/f9.pdf"))
nbZZ_tuned = 300

# f10 - chose nbatches for hmc
f10= bs_hmc(;try_nb=[10,50, 100,200,300,400,500])
savefig(f10, string(SAVEwd, NameEx, "/f10.pdf"))
nbHMC_tuned = 100

# summaries and f11 - plot of the ESS
outESS = ESSsummaries()
outESS


touch(string(SAVEwd, NameEx,"/summaries.txt" ))
io = open(string(SAVEwd, NameEx, "/summaries.txt"), "w")
write(io, outESS["string"])
close(io)
f11=violin(string.(transpose(1:Dim)), (outESS["essZZ"]), side=:right, linewidth=0,
    label="ZZ", color=:black, alpha=0.7, xlabel="Dim", ylabel="ESS", legend=false,
    ylims=(0, maximum(hcat(outESS["essZZ"], outESS["essHMC"]))*1.1))
violin!(string.(transpose(1:Dim)), (outESS["essHMC"]), side=:left, linewidth=0,
    label="HMC", color=:gray)
annotate!([(Dim+0.2, maximum(hcat(outESS["essZZ"], outESS["essHMC"]))*1.05, ("ZZ", 14, :black, :center)),
    (Dim+0.2, maximum(hcat(outESS["essZZ"], outESS["essHMC"]))*0.95, ("HMC", 14, :gray, :center))])
savefig(f11, string(SAVEwd, NameEx, "/f11.pdf"))

# save ESS
JLD.save(string(SAVEwd, NameEx, "/out.jld"),
    Dict("outESS"=>outESS, "zzs"=>zzs, "hmcs"=>hmcs))

# delete everything
f1=f2=f3=f4=f5=f6=f7=f8=f9=f10=f11=nothing
NameEx=io=nothing
tmax_tuned=L_tuned=Lε_tuned=nbZZ_tuned=nbHMC_tuned=start=nothing
zzs=zzsk=zzsm=hmc=hmcs=nothing
Σ=μ=nothing

det(Diagonal(ones(Dim))*0.2+zeros( Dim,Dim).+0.8)
