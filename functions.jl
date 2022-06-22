# Function to compute the gradient of Ux via AD
gradUx(x; Ux=U) = ForwardDiff.gradient(Ux, x)

# rate on one dimension
function rateswitch(t; x0, v0)
    xt= x0+v0*t
    λ = max.(gradUx(xt).*v0, 0.0)
    return λ
end


# global rate
function globalrate(t; x0, v0)
    xt= x0+v0*t
    λ = max.(gradUx(xt).*v0, 0.0)
    return sum(λ)
end

# optimization of the global rate
function myopt(;upper_bound, x0val, v0val, eps=0.000000001, round=true)
    if round==false
        optimΛ  = optimize(t -> - globalrate(t ; x0=x0val, v0=v0val),0,upper_bound)
        Λbar    = - Optim.minimum(optimΛ)
        evals   = Optim.f_calls(optimΛ)
    elseif round==true
        optimΛ  = optimize(t -> - globalrate(t ; x0=x0val, v0=v0val),0,upper_bound, iterations=1, store_trace=true)
        if Optim.converged(optimΛ)==true
            Λbar = - Optim.minimum(optimΛ)
            evals   = Optim.f_calls(optimΛ)
            # println("converged")
        else
            λlower = globalrate(Optim.lower_bound(optimΛ) ; x0=x0val, v0=v0val)
            λcandidate =  - Optim.minimum(optimΛ)
            if (all(y->y==Optim.x_lower_trace(optimΛ)[1], Optim.x_lower_trace(optimΛ)) &&
                λlower >=globalrate(Optim.lower_bound(optimΛ)+eps ; x0=x0val, v0=v0val) &&
                λlower > λcandidate)
                Λbar    = λlower
                evals   = Optim.f_calls(optimΛ)+1
                # println("un converged close to LB")
            else
                λupper = globalrate(Optim.upper_bound(optimΛ) ; x0=x0val, v0=v0val)
                if (all(y->y==Optim.x_upper_trace(optimΛ)[1], Optim.x_upper_trace(optimΛ))&&
                    λupper>=globalrate(Optim.upper_bound(optimΛ)-eps ; x0=x0val, v0=v0val)&&
                    λupper > λcandidate)
                    Λbar = λupper
                    evals   = Optim.f_calls(optimΛ)+1
                else
                    optimΛlong  = optimize(t -> - globalrate(t ; x0=x0val, v0=v0val),0,upper_bound)
                    Λbar = - Optim.minimum(optimΛlong)
                    evals   = Optim.f_calls(optimΛlong)+1
                end
            end
        end
    else
        error("round should be either true or false")
    end
    output=Dict([("Λbar", Λbar), ("evals", evals)])
    return output
end


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
            else        # propose a time for the switch
                tp = rand(Exponential(1/Λ̄))
                if tp >= horizon    # move deterministically if reached the horizon
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


# function to sample points from the skeleton
function zzsample(;N, sk)
    ts  =  sk["SkeletonTime"][1,:]
    vs  =  sk["SkeletonVelocity"][: , :]
    xs  =  sk["SkeletonLocation"][: , :]
    smpl= Array{Float64, 2}(undef, size(xs, 1), N)

    # number of switching times (including 0)
    K   = size(ts,1)
    tm  = (ts[K]/N).*(1:N) # times of the sample
    for i in 1:N
        tm_i=tm[i]
        idx_i = last(findall(i-> i<=tm_i, ts))
        smpl[:, i] = xs[:,idx_i]+vs[:, idx_i].*(tm[i]-ts[idx_i])
    end
    return transpose(smpl)
end



function runHMC(;
    epsilon::Float64,
    L::Int64,
    IT::Int64,
    qs,
    diagnose=false,
    Ux=U
    )
    dim = size(qs)[1]
    μmom = zeros(dim)
    Σmom = Diagonal(ones(dim))+zeros( dim,dim)
    rmom = MvNormal(μmom, Σmom)

    sample  = Array{Float64, 2}(undef, IT+1, dim)
    acccept = zeros( IT+1, 1)
    current_q=qs
    sample[1,:] = current_q
    if diagnose
        qvals      = Array{Float64, 3}(undef, L+1, dim, IT)
        pvals      = Array{Float64, 3}(undef, L+1, dim, IT)
    end

    for i in 1:IT
        # set the initial position
        q = current_q
        p = vec(rand(rmom,1))
        current_p = p

        if diagnose
            # ** plotting **
            qvals      = Array{Float64, 2}(undef, L+1, dim)
            qvals[1,:] = q
            pvals      = Array{Float64, 2}(undef, L+1, dim)
            pvals[1,:] = p
            # ** plotting **
        end

        # half step for the momentum
        p = p - epsilon * gradUx(q)/2
        # alternate full steps positions and momentum
        for l=1:(L-1)
            q=q+epsilon*p
            p=p-epsilon*gradUx(q)
            if diagnose
                qvals[l+1,:] = q
                pvals[l+1,:] = p
            end
        end
        # for the L-th step make a full step for the position
        q=q+epsilon*p
        # and a half step for the momentum
        p = p - epsilon * gradUx(q)/2
        if diagnose
            qvals[L+1,:] = q
            pvals[L+1,:] = p
        end
        # negate the momentum
        p= -p       #not really needed in practice
        # compute the acceptance rate
        current_U  = Ux(current_q)
        current_K  = sum(current_p.^2)/2
        proposed_U = Ux(q)
        proposed_K = sum(p.^2)/2
        if ((rand(Uniform(0,1),1)[1]) < exp(current_U-proposed_U+current_K-proposed_K))
            current_q   = q
            acccept[i+1]= 1
            if diagnose
                plot!(qvals[:,1],qvals[:,dim], color=:green, legend=false)
                plot!([qvals[L+1,1]],[qvals[L+1,dim]], color=:black,seriestype=:scatter,  legend=false)
            end
        else
            acccept[i+1]= 0
            if diagnose
                plot!(qvals[:,1],qvals[:,dim], color=:red, legend=false)
                plot!([qvals[L+1,1]],[qvals[L+1,dim]], color=:black,seriestype=:scatter,  legend=false)
            end
        end
        sample[i+1,:] = current_q

        # global qvals=qvals
        # global pvals=pvals
    end
    output=Dict([("SampleQ", sample), ("accept", acccept)])
    return output
end



function zzsummaries(; dms, # dimension
    sk,  #skeleton of the zigzag as a Dictionary
    B)   #number of batches for ESS
    ts  =  sk["SkeletonTime"][1,:]
    vs  =  sk["SkeletonVelocity"][dms, :]
    xs  =  sk["SkeletonLocation"][dms, :]

    # number of switching times (including 0)
    K   = size(ts,1)

    # first moment
    fms = Array{Float64, 1}(undef, K-1) #factors of the first moment
    for k in 2:K
        fms[k-1] = 0.5*(ts[k]-ts[k-1])*(vs[k-1]*(ts[k]-ts[k-1])+2*xs[k-1])
    end
    k=nothing
    FM = (1/ts[K])*sum(fms)

    # second moment
    sms = Array{Float64, 1}(undef, K-1) #factors of the first moment
    for k in 2:K
        sms[k-1] = ((vs[k-1]*(ts[k]-ts[k-1])+xs[k-1])^3-(xs[k-1]^3))/(3*vs[k-1])
    end
    k=nothing
    SM = (1/ts[K])*sum(sms)

    # Variance
    VAR= SM-FM^2

    # USE OF BATCH MEANS FOR SAMPLE VARIANCE
    FM = Array{Float64, 1}(undef, B)
    for i in 1:B

        # limits of the interval
        tau1_i = (ts[K]/B)*(i-1)
        tau2_i = (ts[K]/B)*(i)

        # plot(sk1["SkeletonTime"][1,:], sk1["SkeletonLocation"][dms,:])
        # vline!([tau1_i, tau2_i])

        # index of the skeleton points contained in the interval
        idx_i = intersect(findall(i-> i>tau1_i, ts),findall(i-> i<=tau2_i, ts))

        ni_i = size(idx_i, 1)+1

        #  lower limit of the integral
        a_i  = Array{Float64, 1}(undef, ni_i)
        a_i[1] =tau1_i
        a_i[2:ni_i] =ts[idx_i]

        #  upper limit of the integral
        b_i  = Array{Float64, 1}(undef, ni_i)
        b_i[1:(ni_i-1)] =ts[idx_i]
        b_i[ni_i] = tau2_i

        # velocities
        v_i  = Array{Float64, 1}(undef, ni_i)
        v_i[1] =vs[last(findall(i-> i<=tau1_i, ts))]
        v_i[2:ni_i] = vs[idx_i]

        # locationsxs
        x_i  = Array{Float64, 1}(undef, ni_i)
        x_i[1] =xs[last(findall(i-> i<=tau1_i, ts))]
        x_i[2:ni_i] = xs[idx_i]

        # initial skeleton point of the current interval
        t0_i = Array{Float64, 1}(undef, ni_i)
        t0_i[1] =ts[last(findall(i-> i<=tau1_i, ts))]
        t0_i[2:ni_i] =ts[idx_i]

        # elements to be summed to obtain the first moment of the batch i
        fms_i = Array{Float64, 1}(undef, ni_i)
        for n in 1:ni_i
            fms_i[n]= -0.5*((a_i[n]-b_i[n])*(v_i[n]*(a_i[n]+b_i[n]-2*t0_i[n])+2*x_i[n]))
        end
        FM[i]= sqrt(B/ts[K])*sum(fms_i)


    end
    # SAMPLE
    SGHAT = 1/(B-1) *sum((FM .- (sum(FM)/B)).^2)

    ESS = ts[K]*VAR/SGHAT
    output = Dict([("FirstMoment", FM), ("SecondMoment", SM), ("SampleVarianceBM", SGHAT), ("EffectiveSampleSize", ESS), ("Dimension", dms)])
    return output
end



function ESSbm(;smpl, nbatches)
    N=(size(smpl)[1])
    m=nbatches
    k=Int64(floor(N/m))
    # batch means
    BM= Array{Float64, 1}(undef, m)
    for j in 1:m
        BM[j]= 1/k *sum(smpl[((j-1)*k+1):(j*k)])
    end
    MU= (1/m) *sum(BM)
    S2batch= (1/(m-1))*sum((BM.-MU).^2)
    Varhatbatch=S2batch/m
    Varhatglobal=(1/(N-1))*sum((smpl.-mean(smpl)).^2)
    ESS = Varhatglobal/Varhatbatch
    return ESS
end




function ESStailprob(;smpl, nbatches)
    N=(size(smpl)[1])
    m=nbatches
    k=Int64(floor(N/m))
    # tailprobs
    smplneg=smpl
    smplneg[findall(i-> i>0, smpl)]=-smpl[findall(i-> i>0, smpl)]
    tailprob=cdf.(TDist(1),smplneg)
    # batch means
    BM= Array{Float64, 1}(undef, m)
    for j in 1:m
        BM[j]= 1/k *sum(tailprob[((j-1)*k+1):(j*k)])
    end
    MU= (1/m) *sum(BM)
    S2batch= (1/(m-1))*sum((BM.-MU).^2)
    Varhatbatch=S2batch/m
    Varhatglobal=(1/(N-1))*sum((tailprob.-mean(tailprob)).^2)
    ESS = Varhatglobal/Varhatbatch
    return ESS
end

#
# # RATES AND FUNCTIONS FOR THE SUBSAPLING
# function rateswitchj(t; j, x0, v0, grad=∇tot, gradj=∇j,  Uxj=Uj)
#     xt= x0+v0*t
#     λ = max.((ForwardDiff.gradient(x -> Uxj(x; j=j), xt)+grad-gradj[j,:]).*v0, 0.0)
#     return λ
# end


# RATES AND FUNCTIONS FOR THE SUBSAPLING
function rateswitchj(t; j, x0, v0, grad=∇tot, gradj=∇j,  Uxj=Uj)
    xt= x0+v0*t
    λ = max.((ForwardDiff.gradient(x -> Uxj(x; j=j), xt)+grad-
            vec(mapslices(mean, gradj[j,:], dims=1))).*v0, 0.0)
    return λ
end


function getMestPOT(;ss, x, v, tmax, J=length(t_vec), D=Dim, SS_s #sub sample size
        )
    # ss=30*20
    # x=[-0.5, 9, -0.5]
    # v=[1,1,1]
    # tmax=0.1
    # J=length(t_vec)
    # D=Dim
    # get the index of the observations used to estimate the max
    idsamp = reshape(sample(1:J,ss*SS_s), ss, SS_s)
    # maximum rate per each dimension (this is unclean might include some nan)
    λ̄uncl     = zeros(ss, D)
    for s in 1:ss
        for d in 1:D
            Optd_s= optimize(t -> -rateswitchj(t; j=idsamp[s,:], x0=x, v0=v)[d], 0, tmax, iterations=1)
            λ̄uncl[s, d]= maximum([-Optim.minimum(Optd_s),
                    rateswitchj(0; j=idsamp[s,:], x0=x, v0=v)[d],
                    rateswitchj(tmax; j=idsamp[s,:], x0=x, v0=v)[d]])
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
        λ̄uncl     = zeros(ss, D)
        for s in 1:ss
            for d in 1:D
                Optd_s= optimize(t -> -rateswitchj(t; j=idsamp[s,:], x0=x, v0=v)[d], 0, tmax, iterations=1)
                λ̄uncl[s, d]= maximum([-Optim.minimum(Optd_s),
                        rateswitchj(0; j=idsamp[s,:], x0=x, v0=v)[d],
                        rateswitchj(tmax; j=idsamp[s,:], x0=x, v0=v)[d]])
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


function zz_w_ss(; NS,      # number of skeleton points
                   x0_0,      # initial location
                   v0_0=missing,    # initial velocity
                   tmax,          # tmax for tuning
                   B=false,        # budget available, if false run until NS
                   # roundopt=true, # whether to use our rounded version of the optimization
                   ε₁ = 1.00e-20,   # boundary error : how close to get before turning
                   ssM =2000,     #number of observations to be used to get the maximum
                   NOBS,           # total number of observation
                   ssS)            #number of subsamples to be taken
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
            M   = getMestPOT(ss=ssM, x=x0i, v=v0i, tmax=horizon, J=NOBS, SS_s=ssS)
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
                        j_0  = sample(1:NOBS, ssS)
                        ar = (rateswitchj(tp; j=j_0, x0=x0i, v0=v0i)[i0])/M[i0]
                        # count evaluations [added to the thinned section]
                        GradEvals[3, k] = GradEvals[3, k] + 1
                        if ar > 1   # if optimization was wrong
                            horizon = tp
                            M   = getMestPOT(ss=ssM, x=x0i, v=v0i, tmax=horizon, J=NOBS, SS_s=ssS)
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
    return output
end
