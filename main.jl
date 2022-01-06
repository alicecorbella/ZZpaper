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
savefig(f2, string(SAVEwd, NameEx, "/f2.pdf"))
tmax_tuned=2.5

# f3 - inspection of the zig zag from mode
start=[0,0]
# multiple chains from the center
f3=plot(0,0, alpha=0, xlabel="X₁", ylabel="X₂")
for r in 1:r_insp
    zzsk =  zz(; NS=1000, x0_0=start, tmax=tmax_tuned)
    plot!(zzsk["SK"][:, 2], zzsk["SK"][:, 3])
end
savefig(f3, string(SAVEwd, NameEx, "/f3.pdf"))

# f4 - inspection of the zig zag from tail
f4=tailstart(xrange=[-7, 7], yrange=[-7, 7], nside=10)
savefig(f4, string(SAVEwd, NameEx, "/f4.pdf"))


# f5 - density inspection
zzsk = zz(; NS=100000, x0_0=start, tmax=tmax_tuned)
zzsm = zzsample(;N=10000, sk=zzsk)
f5=plot(zzsm[:, 1], zzsm[:, 2], seriestype=:scatter,
        color=:black, alpha=0.5, markersize=3,
        xlabel="X₁", ylabel="X₂")
contour!(range(-4, 4, length=100),range(-4,4, length=100),[(pdf(MvNormal(μ, Σ), [x,y])) for y in range(-4, 4, length=100), x in range(-4,4, length=100)], fill=(false,cgrad(:grays,[0,0.1,1.0])))
savefig(f5, string(SAVEwd, NameEx, "/f5.pdf"))

# parameters HMC
Lε_tuned=1
L_tuned=2

# f6 - inspection of the HMC from mode
start=[0,0]
# multiple chains from the center
f6=plot(0,0, alpha=0, xlabel="X₁", ylabel="X₂")
for r in 1:r_insp
    hmc = runHMC(;epsilon=Lε_tuned/L_tuned,L=L_tuned,IT=1000,qs=start)
    plot!(hmc["SampleQ"][:, 1], hmc["SampleQ"][:, 2])
end
savefig(f6, string(SAVEwd, NameEx, "/f6.pdf"))

# f7 - inspection of the HMC from far
f7=tailstart_hmc(xrange=[-7, 7], yrange=[-7, 7], nside=10)
savefig(f7, string(SAVEwd, NameEx, "/f7.pdf"))

# f8 - density hmc
hmc = runHMC(;epsilon=Lε_tuned/L_tuned,L=L_tuned,IT=100000,qs=start)
f8=plot(hmc["SampleQ"][1:10:100000, 1], hmc["SampleQ"][1:10:100000, 2], seriestype=:scatter,
        color=:black, alpha=0.5, markersize=3, xlabel="X₁", ylabel="X₂")
contour!(range(-4, 4, length=100),range(-4,4, length=100),[(pdf(MvNormal(μ, Σ), [x,y])) for y in range(-4, 4, length=100), x in range(-4,4, length=100)], fill=(false,cgrad(:grays,[0,0.1,1.0])))
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
    label="ZZ", color=:blue, alpha=0.7, xlabel="Dim", ylabel="ESS", legend=false,
    ylims=(0, maximum(hcat(outESS["essZZ"], outESS["essHMC"]))*1.1))
violin!(string.(transpose(1:Dim)), (outESS["essHMC"]), side=:left, linewidth=0,
    label="HMC", color=:red)
annotate!([(Dim+0.2, maximum(hcat(outESS["essZZ"], outESS["essHMC"]))*1.05, ("ZZ", 14, :blue, :center)),
    (Dim+0.2, maximum(hcat(outESS["essZZ"], outESS["essHMC"]))*0.95, ("HMC", 14, :red, :center))])
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
savefig(f2, string(SAVEwd, NameEx, "/f2.pdf"))
tmax_tuned=1

# f3 - inspection of the zig zag from mode
start=[0,0]
# multiple chains from the center
f3=plot(0,0, alpha=0, xlabel="X₁", ylabel="X₂")
for r in 1:r_insp
    zzsk =  zz(; NS=1000, x0_0=start, tmax=tmax_tuned)
    plot!(zzsk["SK"][:, 2], zzsk["SK"][:, 3])
end
savefig(f3, string(SAVEwd, NameEx, "/f3.pdf"))

# f4 - inspection of the zig zag from tail
f4=tailstart(xrange=[-7, 7], yrange=[-7, 7], nside=10)
savefig(f4, string(SAVEwd, NameEx, "/f4.pdf"))


# f5 - density inspection
zzsk = zz(; NS=100000, x0_0=start, tmax=tmax_tuned)
zzsm = zzsample(;N=10000, sk=zzsk)
f5=plot(zzsm[:, 1], zzsm[:, 2], seriestype=:scatter,
        color=:black, alpha=0.5, markersize=3,
        xlabel="X₁", ylabel="X₂")
contour!(range(-4, 4, length=100),range(-4,4, length=100),[(pdf(MvNormal(μ, Σ), [x,y])) for y in range(-4, 4, length=100), x in range(-4,4, length=100)], fill=(false,cgrad(:grays,[0,0.1,1.0])))
savefig(f5, string(SAVEwd, NameEx, "/f5.pdf"))

# parameters HMC
Lε_tuned=1.5
L_tuned=3

# f6 - inspection of the HMC from mode
start=[0,0]
# multiple chains from the center
f6=plot(0,0, alpha=0, xlabel="X₁", ylabel="X₂")
for r in 1:r_insp
    hmc = runHMC(;epsilon=Lε_tuned/L_tuned,L=L_tuned,IT=1000,qs=start)
    plot!(hmc["SampleQ"][:, 1], hmc["SampleQ"][:, 2])
end
savefig(f6, string(SAVEwd, NameEx, "/f6.pdf"))

# f7 - inspection of the HMC from far
f7=tailstart_hmc(xrange=[-7, 7], yrange=[-7, 7], nside=10)
savefig(f7, string(SAVEwd, NameEx, "/f7.pdf"))

# f8 - density hmc
hmc = runHMC(;epsilon=Lε_tuned/L_tuned,L=L_tuned,IT=100000,qs=start)
f8=plot(hmc["SampleQ"][1:10:100000, 1], hmc["SampleQ"][1:10:100000, 2], seriestype=:scatter,
        color=:black, alpha=0.5, markersize=3, xlabel="X₁", ylabel="X₂")
contour!(range(-4, 4, length=100),range(-4,4, length=100),[(pdf(MvNormal(μ, Σ), [x,y])) for y in range(-4, 4, length=100), x in range(-4,4, length=100)], fill=(false,cgrad(:grays,[0,0.1,1.0])))
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
    label="ZZ", color=:blue, alpha=0.7, xlabel="Dim", ylabel="ESS", legend=false,
    ylims=(0, maximum(hcat(outESS["essZZ"], outESS["essHMC"]))*1.1))
violin!(string.(transpose(1:Dim)), (outESS["essHMC"]), side=:left, linewidth=0,
    label="HMC", color=:red)
annotate!([(Dim+0.2, maximum(hcat(outESS["essZZ"], outESS["essHMC"]))*1.05, ("ZZ", 14, :blue, :center)),
    (Dim+0.2, maximum(hcat(outESS["essZZ"], outESS["essHMC"]))*0.95, ("HMC", 14, :red, :center))])
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
savefig(f2, string(SAVEwd, NameEx, "/f2.pdf"))
tmax_tuned=4

# f3 - inspection of the zig zag from mode
start=[0,0]
# multiple chains from the center
f3=plot(0,0, alpha=0, xlabel="X₁", ylabel="X₂")
for r in 1:r_insp
    zzsk =  zz(; NS=1000, x0_0=start, tmax=tmax_tuned)
    plot!(zzsk["SK"][:, 2], zzsk["SK"][:, 3])
end
savefig(f3, string(SAVEwd, NameEx, "/f3.pdf"))

# f4 - inspection of the zig zag from tail
f4=tailstart(xrange=[-7, 7], yrange=[-70, 70], nside=10)
savefig(f4, string(SAVEwd, NameEx, "/f4.pdf"))


# f5 - density inspection
zzsk = zz(; NS=100000, x0_0=start, tmax=tmax_tuned)
zzsm = zzsample(;N=10000, sk=zzsk)
f5=plot(zzsm[:, 1], zzsm[:, 2], seriestype=:scatter,
        color=:black, alpha=0.5, markersize=3,
        xlabel="X₁", ylabel="X₂")
contour!(range(-4, 4, length=100),range(-40,40, length=100),[(pdf(MvNormal(μ, Σ), [x,y])) for y in range(-40, 40, length=100), x in range(-4,4, length=100)], fill=(false,cgrad(:grays,[0,0.1,1.0])))
savefig(f5, string(SAVEwd, NameEx, "/f5.pdf"))

# parameters HMC
Lε_tuned=1.5
L_tuned=2

# f6 - inspection of the HMC from mode
start=[0,0]
# multiple chains from the center
f6=plot(0,0, alpha=0, xlabel="X₁", ylabel="X₂")
for r in 1:r_insp
    hmc = runHMC(;epsilon=Lε_tuned/L_tuned,L=L_tuned,IT=1000,qs=start)
    plot!(hmc["SampleQ"][:, 1], hmc["SampleQ"][:, 2])
end
savefig(f6, string(SAVEwd, NameEx, "/f6.pdf"))

# f7 - inspection of the HMC from far
f7=tailstart_hmc(xrange=[-7, 7], yrange=[-70, 70], nside=10)
savefig(f7, string(SAVEwd, NameEx, "/f7.pdf"))

# f8 - density hmc
hmc = runHMC(;epsilon=Lε_tuned/L_tuned,L=L_tuned,IT=100000,qs=start)
f8=plot(hmc["SampleQ"][1:10:100000, 1], hmc["SampleQ"][1:10:100000, 2], seriestype=:scatter,
        color=:black, alpha=0.5, markersize=3, xlabel="X₁", ylabel="X₂")
contour!(range(-4, 4, length=100),range(-40,40, length=100),[(pdf(MvNormal(μ, Σ), [x,y])) for y in range(-40, 40, length=100), x in range(-4,4, length=100)], fill=(false,cgrad(:grays,[0,0.1,1.0])))
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
    label="ZZ", color=:blue, alpha=0.7, xlabel="Dim", ylabel="ESS", legend=false,
    ylims=(0, maximum(hcat(outESS["essZZ"], outESS["essHMC"]))*1.1))
violin!(string.(transpose(1:Dim)), (outESS["essHMC"]), side=:left, linewidth=0,
    label="HMC", color=:red)
annotate!([(Dim+0.2, maximum(hcat(outESS["essZZ"], outESS["essHMC"]))*1.05, ("ZZ", 14, :blue, :center)),
    (Dim+0.2, maximum(hcat(outESS["essZZ"], outESS["essHMC"]))*0.95, ("HMC", 14, :red, :center))])
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
savefig(f2, string(SAVEwd, NameEx, "/f2.pdf"))
tmax_tuned=1.5

# f3 - inspection of the zig zag from mode
start=[0,0]
# multiple chains from the center
f3=plot(0,0, alpha=0, xlabel="X₁", ylabel="X₂")
for r in 1:r_insp
    zzsk =  zz(; NS=1000, x0_0=start, tmax=tmax_tuned)
    plot!(zzsk["SK"][:, 2], zzsk["SK"][:, 3])
end
savefig(f3, string(SAVEwd, NameEx, "/f3.pdf"))

# f4 - inspection of the zig zag from tail
f4=tailstart(xrange=[-7, 7], yrange=[-7, 7], nside=10)
savefig(f4, string(SAVEwd, NameEx, "/f4.pdf"))


# f5 - density inspection
zzsk = zz(; NS=100000, x0_0=start, tmax=tmax_tuned)
zzsm = zzsample(;N=10000, sk=zzsk)
f5=plot(zzsm[:, 1], zzsm[:, 2], seriestype=:scatter,
        color=:black, alpha=0.5, markersize=3,
        xlabel="X₁", ylabel="X₂")
contour!(range(-4, 4, length=100),range(-4,4, length=100),[(gx([x,y])) for y in range(-4, 4, length=100), x in range(-4,4, length=100)], fill=(false,cgrad(:grays,[0,0.1,1.0])))
savefig(f5, string(SAVEwd, NameEx, "/f5.pdf"))

# parameters HMC
Lε_tuned=1
L_tuned=2

# f6 - inspection of the HMC from mode
start=[0,0]
# multiple chains from the center
f6=plot(0,0, alpha=0, xlabel="X₁", ylabel="X₂")
for r in 1:r_insp
    hmc = runHMC(;epsilon=Lε_tuned/L_tuned,L=L_tuned,IT=1000,qs=start)
    plot!(hmc["SampleQ"][:, 1], hmc["SampleQ"][:, 2])
end
savefig(f6, string(SAVEwd, NameEx, "/f6.pdf"))

# f7 - inspection of the HMC from far
f7=tailstart_hmc(xrange=[-7, 7], yrange=[-7, 7], nside=10)
f7=tailend_hmc(xrange=[-7, 7], yrange=[-7, 7], nside=10)
savefig(f7, string(SAVEwd, NameEx, "/f7.pdf"))

# f8 - density hmc
hmc = runHMC(;epsilon=Lε_tuned/L_tuned,L=L_tuned,IT=100000,qs=start)
f8=plot(hmc["SampleQ"][1:10:100000, 1], hmc["SampleQ"][1:10:100000, 2], seriestype=:scatter,
        color=:black, alpha=0.5, markersize=3, xlabel="X₁", ylabel="X₂")
contour!(range(-4, 4, length=100),range(-4,4, length=100),[(gx([x,y])) for y in range(-4, 4, length=100), x in range(-4,4, length=100)], fill=(false,cgrad(:grays,[0,0.1,1.0])))
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
nbZZ_tuned = 500

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
    label="ZZ", color=:blue, alpha=0.7, xlabel="Dim", ylabel="ESS", legend=false,
    ylims=(0, maximum(hcat(outESS["essZZ"], outESS["essHMC"]))*1.1))
violin!(string.(transpose(1:Dim)), (outESS["essHMC"]), side=:left, linewidth=0,
    label="HMC", color=:red)
annotate!([(Dim+0.2, maximum(hcat(outESS["essZZ"], outESS["essHMC"]))*1.05, ("ZZ", 14, :blue, :center)),
    (Dim+0.2, maximum(hcat(outESS["essZZ"], outESS["essHMC"]))*0.95, ("HMC", 14, :red, :center))])
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
# mkdir(string(SAVEwd, NameEx))
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
savefig(f2, string(SAVEwd, NameEx, "/f2.pdf"))
tmax_tuned=5

# f3 - inspection of the zig zag from mode
start=[0,0]
# multiple chains from the center
f3=plot(0,0, alpha=0, xlabel="X₁", ylabel="X₂")
for r in 1:r_insp
    zzsk =  zz(; NS=1000, x0_0=start, tmax=tmax_tuned)
    plot!(zzsk["SK"][:, 2], zzsk["SK"][:, 3])
end
savefig(f3, string(SAVEwd, NameEx, "/f3.pdf"))

# f4 - inspection of the zig zag from tail
Random.seed!(050620)
f4=tailstart(xrange=[-500, 500], yrange=[-500, 500], nside=10)
savefig(f4, string(SAVEwd, NameEx, "/f4.pdf"))


# f5 - density inspection
zzsk = zz(; NS=100000, x0_0=start, tmax=tmax_tuned)
zzsm = zzsample(;N=10000, sk=zzsk)
f5=plot(zzsm[:, 1], zzsm[:, 2], seriestype=:scatter,
        color=:black, alpha=0.5, markersize=3,
        xlabel="X₁", ylabel="X₂")
contour!(range(-50, 50, length=100),range(-50,50, length=100),[(gx([x,y])) for y in range(-50, 50, length=100), x in range(-50,50, length=100)], fill=(false,cgrad(:grays,[0,0.1,1.0])))
savefig(f5, string(SAVEwd, NameEx, "/f5.pdf"))

# parameters HMC
Lε_tuned=10
L_tuned=10

# f6 - inspection of the HMC from mode
start=[0,0]
# multiple chains from the center
f6=plot(0,0, alpha=0, xlabel="X₁", ylabel="X₂")
for r in 1:r_insp
    hmc = runHMC(;epsilon=Lε_tuned/L_tuned,L=L_tuned,IT=1000,qs=start)
    plot!(hmc["SampleQ"][:, 1], hmc["SampleQ"][:, 2])
end
savefig(f6, string(SAVEwd, NameEx, "/f6.pdf"))

# f7 - inspection of the HMC from far
f7=tailstart_hmc(xrange=[-500, 500], yrange=[-500, 500], nside=10)
savefig(f7, string(SAVEwd, NameEx, "/f7.pdf"))

# f8 - density hmc
hmc = runHMC(;epsilon=Lε_tuned/L_tuned,L=L_tuned,IT=100000,qs=start)
f8=plot(hmc["SampleQ"][1:10:100000, 1], hmc["SampleQ"][1:10:100000, 2], seriestype=:scatter,
        color=:black, alpha=0.5, markersize=3, xlabel="X₁", ylabel="X₂")
contour!(range(-4, 4, length=100),range(-4,4, length=100),[(gx([x,y])) for y in range(-4, 4, length=100), x in range(-4,4, length=100)], fill=(false,cgrad(:grays,[0,0.1,1.0])))
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
    label="ZZ", color=:blue, alpha=0.7, xlabel="Dim", ylabel="ESS", legend=false,
    ylims=(0, maximum(hcat(outESS["essZZ"], outESS["essHMC"]))*1.1))
violin!(string.(transpose(1:Dim)), (outESS["essHMC"]), side=:left, linewidth=0,
    label="HMC", color=:red)
annotate!([(Dim+0.2, maximum(hcat(outESS["essZZ"], outESS["essHMC"]))*1.05, ("ZZ", 14, :blue, :center)),
    (Dim+0.2, maximum(hcat(outESS["essZZ"], outESS["essHMC"]))*0.95, ("HMC", 14, :red, :center))])
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
