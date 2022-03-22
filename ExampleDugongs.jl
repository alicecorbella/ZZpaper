# Direcotries
CODEwd = "/Users/alice/projects/ZZpaper/"
SAVEwd = "/Users/alice/OneDrive - University of Warwick/Manuscripts/04ADZZ/figures/Sec5.1/"
NameEx = ""
# Loading the packages
include(string(CODEwd, "env.jl"))

# Loading the functions
include(string(CODEwd, "functions.jl"))
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


U([12, 12, 12, 12])

U([12, -12, -12, -12])


# f2 - data
f1=plot(data_x, data_y, seriestype=:scatter, legend=false, color=:black,
    alpha=0.5, xlabel="Zⱼ", ylabel="Yⱼ")
savefig(f1, string(SAVEwd, "f1.pdf"))

# f2 - tune of t_max
f2=tmplot(R=100, try_tmax=[0.01, 0.02, 0.03,0.04, 0.05, 0.06])
savefig(f2["NHPP"], string(SAVEwd, NameEx, "/f2a.pdf"))
savefig(f2["opt"], string(SAVEwd, NameEx, "/f2b.pdf"))
savefig(f2["tot"], string(SAVEwd, NameEx, "/f2c.pdf"))
f2tot= plot(f2["NHPP"],f2["opt"],f2["tot"], layout=(1,3),
    size=(pt(252), pt(60)))
savefig(f2tot, string(SAVEwd, NameEx, "/f2.pdf"))
tmax_tuned=0.02

# f3 - inspection of the zig zag from mode

# multiple chains from the center
SKinsp= Array{Float64, 3}(undef, 1000, Dim, r_insp)
TMinsp= Array{Float64, 3}(undef, 1000, 1, r_insp)
for r in 1:r_insp
    start= rand(MvNormal([1,0,2,-2],Diagonal(ones(Dim)).*0.01+zeros( Dim,Dim)))
    zzsk =  zz(; NS=1000, x0_0=start, tmax=tmax_tuned)
    SKinsp[:,:,r] = zzsk["SK"][:, 2:(Dim+1)]
    TMinsp[:,:,r] = zzsk["SK"][:, 1]
end

f3a=plot(0,0, alpha=0, xlabel="Time", ylabel="X₁")
for r in 1:r_insp
    plot!(TMinsp[:, 1, r], exp.(SKinsp[:,1,r]), linewidth=0.5, color=:gray)
end
f3b=plot(0,0, alpha=0, xlabel="Time", ylabel="X₂")
for r in 1:r_insp
    plot!(TMinsp[:, 1, r], exp.(SKinsp[:,2,r]), linewidth=0.5, color=:gray)
end
f3c=plot(0,0, alpha=0, xlabel="Time", ylabel="X₃")
for r in 1:r_insp
    plot!(TMinsp[:, 1, r], exp.(SKinsp[:,3,r])./(1 .+exp.(SKinsp[:,3,r])), linewidth=0.5, color=:gray)
end
f3d=plot(0,0, alpha=0, xlabel="Time", ylabel="X₄")
for r in 1:r_insp
    plot!(TMinsp[:, 1, r], exp.(SKinsp[:,4,r]), linewidth=0.5, color=:gray)
end
# savefig(f3a, string(SAVEwd, "/f3a.pdf"))
# savefig(f3b, string(SAVEwd, "/f3b.pdf"))
# savefig(f3c, string(SAVEwd, "/f3c.pdf"))
# savefig(f3d, string(SAVEwd, "/f3d.pdf"))
f3 = plot(f3a, f3b, f3c, f3d, layout=(2, 2), size=(pt(84*2), pt(84*2)))
# savefig(f3, string(SAVEwd, "/f3.pdf"))



f3a=plot(0,0, alpha=0, xlabel="Time", ylabel="X₁")
for r in 1:r_insp
    plot!(TMinsp[:, 1, r], SKinsp[:,1,r], linewidth=0.5, color=:gray)
end
f3b=plot(0,0, alpha=0, xlabel="Time", ylabel="X₂")
for r in 1:r_insp
    plot!(TMinsp[:, 1, r], SKinsp[:,2,r], linewidth=0.5, color=:gray)
end
f3c=plot(0,0, alpha=0, xlabel="Time", ylabel="X₃")
for r in 1:r_insp
    plot!(TMinsp[:, 1, r], SKinsp[:,3,r], linewidth=0.5, color=:gray)
end
f3d=plot(0,0, alpha=0, xlabel="Time", ylabel="X₄")
for r in 1:r_insp
    plot!(TMinsp[:, 1, r], SKinsp[:,4,r], linewidth=0.5, color=:gray)
end
# savefig(f3a, string(SAVEwd, "/f3a.pdf"))
# savefig(f3b, string(SAVEwd, "/f3b.pdf"))
# savefig(f3c, string(SAVEwd, "/f3c.pdf"))
# savefig(f3d, string(SAVEwd, "/f3d.pdf"))
f3 = plot(f3a, f3b, f3c, f3d, layout=(2, 2), size=(pt(84*2), pt(84*2)))
# savefig(f3, string(SAVEwd, "/f3.pdf"))





function zz(;
    NS,             # number of skeleton points
    x0_0,           # initial location
    v0_0=missing,   # initial velocity
    tmax,           # tmax for tuning
    B=false,        # budget available, if false run until NS
    roundopt=true,  # whether to use our rounded version of the optimization
    ε₁ = 1.00e-20   # boundary error : how close to get before turning
)
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

            # # compute the rate at the current state for weights
            # λ0 = rateswitch(0; x0=x0i, v0=v0i)
            # Λ0 = sum(λ0)
            # # count evaluations [added to the horizon change section]
            # GradEvals[1, k] = GradEvals[1, k] + 1
            # # compute the weights: if global rate is 0 then  sample discrete uniform
            # w  = (Λ0==0) ? ones(length(λ0)) : (λ0/Λ0)
            # # chose dimension to switch
            # m  = sample(1:D, Weights(w))
            # # switch velocity and the retain location (border ahead)
            # v0i[m] = -v0i[m]
            # # save the skeleton point
            # v0set[:, k] = v0i
            # x0set[:, k] = x0i
            # t0set[1, k] = t0set[1, k-1]+ts
            # # reset time from skeleton point, horizon, and increase counter
            # ts = 0.0
            # horizon = tmax
            # k = k+1
        else # continue with the current horizon
            opt = myopt(upper_bound=horizon, x0val=x0i, v0val=v0i, round=roundopt)
            Λ̄   = opt["Λbar"]
            # count evaluations [added to the optimization evaluation section]
            GradEvals[2, k] = GradEvals[2, k] + opt["evals"]
            if Λ̄ == 0   # move dererministically if the rate of switching is 0
                ts = ts+horizon
                x0i= x0i+horizon*v0i
            elseif Λ̄>=(10.0^10)
                # chose dimension to switch
                m  = argmax(rateswitch(0; x0=x0i, v0=v0i))
                # update location and switch velocity
                x0i    = x0i + tp * v0i
                v0i[m] = -v0i[m]
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
            else        # propose a time for the switch
                tp = rand(Exponential(1/Λ̄))
                if tp >= horizon    # move deterministically if reached the horizon
                    ts = ts+horizon
                    x0i= x0i+horizon*v0i
                else                # evaluate proposal
                    accept = false  # start evaluations
                    while (tp < horizon) && (accept == false)
                        λt = rateswitch(tp; x0=x0i, v0=v0i)
                        Λt = sum(λt)
                        ar = Λt/Λ̄
                        # count evaluations [added to the thinned section]
                        GradEvals[3, k] = GradEvals[3, k] + 1
                        if ar > 1   # if optimization was wrong
                            horizon = tp
                            opt = myopt(upper_bound=horizon, x0val=x0i, v0val=v0i, round=roundopt)
                            Λ̄   = opt["Λbar"]
                            # restart with the new horizon/optimum
                            tp  = rand(Exponential(1/Λ̄))
                            # count evaluations [added to the optimization evaluation section]
                            GradEvals[2, k] = GradEvals[2, k] + opt["evals"]
                            ErrorOpt[1, k] = ErrorOpt[1, k] + 1
                        else        # evaluate acceptance
                            if (rand(Bernoulli(ar)))
                                # chose dimension to switch
                                m  = sample(1:D, Weights(λt/Λt))
                                # update location and switch velocity
                                x0i    = x0i + tp * v0i
                                v0i[m] = -v0i[m]
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
                            else   # upon rejection increas stochastic time
                                tp=tp+rand(Exponential(1/Λ̄))
                            end
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
end



# f4 - inspection of the zig zag from tail
Random.seed!(1234)
SKinsp= Array{Float64, 3}(undef, 1000, Dim, 10)
TMinsp= Array{Float64, 3}(undef, 1000, 1, 10)
for r in 1:10
    start = rand([-12,12], 4)
    print(start)
    print(U(start))
    # print(U(startvalues[:, r]))
    zzsk =  zz(; NS=1000,x0_0=start, v0_0=[-1,-1,-1,-1],tmax=tmax_tuned)
    SKinsp[:,:,r] = zzsk["SK"][:, 2:(Dim+1)]
    TMinsp[:,:,r] = zzsk["SK"][:, 1]
    print(r)
end

#
# start = [12, 12, -12,-12]
# print(start)
# print(U(start))
# plot()
# # print(U(startvalues[:, r]))
# zzsk =  zz(; NS=10, x0_0=start, v0_0=[-1,-1,-1,-1],tmax=0.0001)
# plot(zzsk["SK"][:, 1], zzsk["SK"][:, 5])
#
# plot(t -> globalrate(t; x0=[12, 12, -12,-12], v0=[-1,-1,1,1]), xlims=(0,0.00000005))
# plot(t -> rateswitch(t; x0=[12, 12, -12,-12], v0=[-1,-1,1,1])[1], xlims=(0,0.00000005))
# plot(t -> rateswitch(t; x0=[12, 12, -12,-12], v0=[-1,-1,1,1])[2], xlims=(0,0.00000005))
# plot(t -> rateswitch(t; x0=[12, 12, -12,-12], v0=[-1,-1,1,1])[3], xlims=(0,0.00000005))
# plot(t -> rateswitch(t; x0=[12, 12, -12,-12], v0=[-1,-1,1,1])[4], xlims=(0,0.00000005))
#
# plot(t -> globalrate(t; x0=[12, 12, -12,-12], v0=[1,1,1,1]), xlims=(0,0.00000005))
# plot(t -> rateswitch(t; x0=[12, 12, -12,-12], v0=[1,1,1,1])[1], xlims=(0,0.00000005))
# plot(t -> rateswitch(t; x0=[12, 12, -12,-12], v0=[1,1,1,1])[2], xlims=(0,0.00000005))
# plot(t -> rateswitch(t; x0=[12, 12, -12,-12], v0=[1,1,1,1])[3], xlims=(0,0.00000005))
# plot(t -> rateswitch(t; x0=[12, 12, -12,-12], v0=[1,1,1,1])[4], xlims=(0,0.00000005))
#
# plot(t -> globalrate(t; x0=[12, 12, 12, 12], v0=[1,-1,-1,-1]), xlims=(0,.05))
# plot(t -> rateswitch(t; x0=[12, 12, 12, 12], v0=[1,-1,-1,-1])[1], xlims=(0,5))
# plot(t -> rateswitch(t; x0=[12, 12, 12, 12], v0=[1,-1,-1,-1])[2], xlims=(0,5))
# plot(t -> rateswitch(t; x0=[12, 12, 12, 12], v0=[1,-1,-1,-1])[3], xlims=(0,5))
# plot(t -> rateswitch(t; x0=[12, 12, 12, 12], v0=[1,-1,-1,-1])[4], xlims=(0,5))
#
#
# plot(t -> globalrate(t; x0=[12, 12, 12, 12], v0=[-1,-1,-1,-1]), xlims=(0,.05))
# plot(t -> rateswitch(t; x0=[12, 12, 12, 12], v0=[-1,-1,-1,-1])[1], xlims=(0,5))
# plot(t -> rateswitch(t; x0=[12, 12, 12, 12], v0=[-1,-1,-1,-1])[2], xlims=(0,5))
# plot(t -> rateswitch(t; x0=[12, 12, 12, 12], v0=[-1,-1,-1,-1])[3], xlims=(0,5))
# plot(t -> rateswitch(t; x0=[12, 12, 12, 12], v0=[-1,-1,-1,-1])[4], xlims=(0,5))
#
#
# plot(t -> globalrate(t; x0=[12, 12, -12,-12], v0=[-1,-1,-1,-1]), xlims=(0,0.00000005))
# plot(t -> globalrate(t; x0=[12, 12, -12,-12], v0=[1,1,-1,-1]), xlims=(0,0.00000005))
# plot(t -> globalrate(t; x0=[12, 12, -12,-12], v0=[-1,1,-1,1]), xlims=(0,0.00000005))
# plot(t -> globalrate(t; x0=[12, 12, -12,-12], v0=[-1,1,1,-1]), xlims=(0,0.00000005))
#
#
# rateswitch(0.002; x0=[12, 12, -12,-12], v0=[-1,-1,-1,-1])./ globalrate(0.002; x0=[12, 12, -12,-12], v0=[-1,-1,-1,-1])
#


f4a=plot(0,0, alpha=0, xlabel="Time", ylabel="X₁")
for r in 1:r_insp
    plot!(TMinsp[:, 1, r], SKinsp[:,1,r], linewidth=0.5)
end
f4b=plot(0,0, alpha=0, xlabel="Time", ylabel="X₂")
for r in 1:r_insp
    plot!(TMinsp[:, 1, r], SKinsp[:,2,r], linewidth=0.5)
end
f4c=plot(0,0, alpha=0, xlabel="Time", ylabel="X₃")
for r in 1:r_insp
    plot!(TMinsp[:, 1, r], SKinsp[:,3,r], linewidth=0.5)
end
f4d=plot(0,0, alpha=0, xlabel="Time", ylabel="X₄")
for r in 1:r_insp
    plot!(TMinsp[:, 1, r], SKinsp[:,4,r], linewidth=0.5)
end
# savefig(f4a, string(SAVEwd, "/f4aWSC.pdf"))
# savefig(f4b, string(SAVEwd, "/f4bWSC.pdf"))
# savefig(f4c, string(SAVEwd, "/f4cWSC.pdf"))
# savefig(f4d, string(SAVEwd, "/f4dWSC.pdf"))
f4 = plot(f4a, f4b, f4c, f4d, layout=(2, 2), size=(pt(84*2), pt(84*2)))
# savefig(f4, string(SAVEwd, "/f4WSC.pdf"))


# f5 - density inspection
zzsk = zz(; NS=100000, x0_0=rand(MvNormal([1,0,2,-2],Diagonal(ones(Dim)).*0.01+zeros( Dim,Dim))),
    tmax=tmax_tuned)
zzsm = zzsample(;N=10000, sk=zzsk)
f5a=density(zzsm[1000:10000, 1], xlabel="X₁", color=:black, linewidth=3, alpha=0.7)
f5b=density(zzsm[1000:10000, 2], xlabel="X₂", color=:black, linewidth=3, alpha=0.7)
f5c=density(zzsm[1000:10000, 3], xlabel="X₃", color=:black, linewidth=3, alpha=0.7)
f5d=density(zzsm[1000:10000, 4], xlabel="X₄", color=:black, linewidth=3, alpha=0.7)
f5 = plot(f5a, f5b, f5c, f5d, layout=(2, 2), size=(pt(84*2), pt(84*2)))
savefig(f5, string(SAVEwd, "/f5.pdf"))


# parameters HMC
Lε_tuned=0.1
L_tuned=10

hmc = runHMC(;epsilon=Lε_tuned/L_tuned,L=L_tuned,IT=1000,qs=[1,0,2,-2])
plot(plot(hmc["SampleQ"][:, 1]), plot(hmc["SampleQ"][:, 2]),
     plot(hmc["SampleQ"][:, 3]), plot(hmc["SampleQ"][:, 4]),
        layout=(4,1), size=(1000,1500))
ac=plot(autocor(hmc["SampleQ"][:, 1]))
plot!(autocor(hmc["SampleQ"][:, 2]))
plot!(autocor(hmc["SampleQ"][:, 3]))
plot!(autocor(hmc["SampleQ"][:, 4]))
sum(hmc["accept"])/1000



# f6 - inspection of the HMC from mode
# multiple chains from the center
HMCinsp= Array{Float64, 3}(undef, 1000, Dim, r_insp)
for r in 1:r_insp
    start= rand(MvNormal([1,0,2,-2],Diagonal(ones(Dim)).*0.01+zeros( Dim,Dim)))
    hmc = runHMC(;epsilon=Lε_tuned/L_tuned,L=L_tuned,IT=1000,qs=start)
    HMCinsp[:,:,r] = hmc["SampleQ"][1:1000, 1:(Dim)]
end

f6a=plot(0,0, alpha=0, xlabel="Iterations", ylabel="X₁")
for r in 1:r_insp
    plot!(HMCinsp[:,1,r], linewidth=0.5, color=:gray)
end
f6b=plot(0,0, alpha=0, xlabel="Iterations", ylabel="X₂")
for r in 1:r_insp
    plot!(HMCinsp[:,2,r], linewidth=0.5, color=:gray)
end
f6c=plot(0,0, alpha=0, xlabel="Iterations", ylabel="X₃")
for r in 1:r_insp
    plot!(HMCinsp[:,3,r], linewidth=0.5, color=:gray)
end
f6d=plot(0,0, alpha=0, xlabel="Iterations", ylabel="X₄")
for r in 1:r_insp
    plot!(HMCinsp[:,4,r], linewidth=0.5, color=:gray)
end
f6 = plot(f6a, f6b, f6c, f6d, layout=(2, 2), size=(pt(84*2), pt(84*2)))
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
    start = rand([-4,4], 4)
    print(U(start))
    hmc = runHMC(;epsilon=Lε_tuned/L_tuned,L=L_tuned,IT=1000,qs=start)
    HMCinsp[:,:,r] = hmc["SampleQ"][1:1000, 1:(Dim)]
end

HMCinsp

f7a=plot(0,0, alpha=0, xlabel="Time", ylabel="X₁")
for r in 1:r_insp
    plot!(HMCinsp[:,1,r], linewidth=0.5, color=:gray)
end
f7b=plot(0,0, alpha=0, xlabel="Time", ylabel="X₂")
for r in 1:r_insp
    plot!(HMCinsp[:,2,r], linewidth=0.5, color=:gray)
end
f7c=plot(0,0, alpha=0, xlabel="Time", ylabel="X₃")
for r in 1:r_insp
    plot!(HMCinsp[:,3,r], linewidth=0.5, color=:gray)
end
f7d=plot(0,0, alpha=0, xlabel="Time", ylabel="X₄")
for r in 1:r_insp
    plot!(HMCinsp[:,4,r], linewidth=0.5, color=:gray)
end
savefig(f7a, string(SAVEwd, "/f7a.pdf"))
savefig(f7b, string(SAVEwd, "/f7b.pdf"))
savefig(f7c, string(SAVEwd, "/f7c.pdf"))
savefig(f7d, string(SAVEwd, "/f7d.pdf"))



f7 = plot(f7a, f7b, f7c, f7d, layout=(2, 2), size=(pt(84*2), pt(84*2)))
savefig(f7, string(SAVEwd, "/f7.pdf"))

# f8 - density hmc
hmc = runHMC(;epsilon=Lε_tuned/L_tuned,L=L_tuned,IT=100000,
    qs=rand(MvNormal([1,0,2,-2],Diagonal(ones(Dim)).*0.01+zeros( Dim,Dim))))
f8a=density(hmc["SampleQ"][1000:10:100000, 1], xlabel="X₁", color=:black, linewidth=3, alpha=0.7)
f8b=density(hmc["SampleQ"][1000:10:100000, 2], xlabel="X₂", color=:black, linewidth=3, alpha=0.7)
f8c=density(hmc["SampleQ"][1000:10:100000, 3], xlabel="X₃", color=:black, linewidth=3, alpha=0.7)
f8d=density(hmc["SampleQ"][1000:10:100000, 4], xlabel="X₄", color=:black, linewidth=3, alpha=0.7)
f8 = plot(f8a, f8b, f8c, f8d, layout=(2, 2), size=(pt(84*2), pt(84*2)))
savefig(f8, string(SAVEwd, "/f8.pdf"))



# f8 - density hmc
f58a=density(hmc["SampleQ"][1000:10:100000, 1], xlabel="X₁", color=:red, linewidth=3, alpha=0.7)
density!(zzsm[1000:10000, 1], xlabel="θ₁", color=:blue, linewidth=3, alpha=0.7)
f58b=density(hmc["SampleQ"][1000:10:100000, 2], xlabel="X₂", color=:red, linewidth=3, alpha=0.7)
density!(zzsm[1000:10000, 2], xlabel="θ₂", color=:blue, linewidth=3, alpha=0.7)
f58c=density(hmc["SampleQ"][1000:10:100000, 3], xlabel="X₃", color=:red, linewidth=3, alpha=0.7)
density!(zzsm[1000:10000, 3], xlabel="θ₃", color=:blue, linewidth=3, alpha=0.7)
f58d=density(hmc["SampleQ"][1000:10:100000, 4], xlabel="X₄", color=:red, linewidth=3, alpha=0.7)
density!(zzsm[1000:10000, 4], xlabel="θ₄", color=:blue, linewidth=3, alpha=0.7)
f58 = plot(f58a, f58b, f58c, f58d, layout=(2, 2), size=(pt(84*2), pt(84*2)))
savefig(f58, string(SAVEwd, "/f58.pdf"))


zzs = Dict()
hmcs= Dict()
for r in 1:r_ess
    start=rand(MvNormal([1,0,2,-2],Diagonal(ones(Dim)).*0.01+zeros( Dim,Dim)))
    zzs[[r]]=zz(NS=false, B=budget, x0_0=start, tmax=tmax_tuned)
    hmcs[[r]]=runHMC(;epsilon=Lε_tuned/L_tuned,L=L_tuned,IT=Int(round(budget/(L_tuned+1))),qs=start)
    print(r)
end

# f9 - chose nbatches for zz
f9= bs_zz(;try_nb=[5, 10,25, 50, 100, 500])
savefig(f9, string(SAVEwd, "/f9.pdf"))
nbZZ_tuned = 50

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
    label="ZZ", color=:black, xlabel="Dim", ylabel="ESS", legend=false,
    ylims=(0, maximum(hcat(outESS["essZZ"], outESS["essHMC"]))*1.1))
violin!(string.(transpose(1:Dim)), (outESS["essHMC"]), side=:left, linewidth=0,
    label="HMC", color=:gray)
annotate!([(Dim+0.2, maximum(hcat(outESS["essZZ"], outESS["essHMC"]))*1.05, ("ZZ", 14, :black, :center)),
    (Dim+0.2, maximum(hcat(outESS["essZZ"], outESS["essHMC"]))*0.95, ("HMC", 14, :gray, :center))])
savefig(f11, string(SAVEwd, "/f11.pdf"))

# save ESS
JLD.save(string(SAVEwd, "/out.jld"),
    Dict("outESS"=>outESS, "zzs"=>zzs, "hmcs"=>hmcs))
