

function densplot(; Ux=U, xrange=[-1.5,1.5], yrange=[-1.5,1.5], col=false)
    outplot=contour(range(xrange[1], xrange[2], length=100),
            range(yrange[1],yrange[2], length=100),
        [(-exp(Ux([x,y]))) for y in range(xrange[1], xrange[2], length=100), x in range(yrange[1],yrange[2], length=100)],
        color=col)
    return outplot
end

function tmplot(;R = 50, try_tmax, ns=1000, start=missing)
    ltm = length(try_tmax)
    tpp_tmax  = Array{Float64, 2}(undef, R, ltm)
    opt_tmax  = Array{Float64, 2}(undef, R, ltm)
    for t in 1:ltm
        for i in 1:R
            if ismissing(start)
                xstart=rand(Dim)
            else
                xstart=start
            end
            sk1 = zz(; NS=ns, x0_0=xstart,tmax=try_tmax[t])
            tpp_tmax[i,t] =sum(sk1["GradientEvaluations"][2,:])
            opt_tmax[i,t] =sum(sk1["GradientEvaluations"][3,:])
        end
    end
    p1=violin(string.(transpose(try_tmax)), tpp_tmax, color=:blue, alpha=0.7,
        title=("Gradient evaluations for \nthe NHPP via thinning"),
        ylabel="count ∇U(x)", xlabel="tₘₐₓ")
    # plot!(string.(transpose(try_tmax)), transpose(mapslices(mean, tpp_tmax, dims=1)), color=:black, linewidth=2)
    p2=violin(string.(transpose(try_tmax)), opt_tmax, color=:blue, alpha=0.7,
        title=("Gradient evaluations for \nthe optimization routine"),
        ylabel="count ∇U(x)", xlabel="tₘₐₓ")
    # plot!(try_tmax, transpose(mapslices(mean, opt_tmax, dims=1)), color=:black, linewidth=2)
    p3=violin(string.(transpose(try_tmax)), opt_tmax+tpp_tmax, color=:blue, alpha=0.7,
        title=("Total gradient evaluations"),
        ylabel="count ∇U(x)", xlabel="tₘₐₓ")
    # plot!(try_tmax, transpose(mapslices(mean, opt_tmax+tpp_tmax, dims=1)), color=:black, linewidth=2)
    outplot = plot(p1,p2,p3, layout=(1, 3), size = (1200, 400))
    return outplot
end


function tailstart(; xrange, yrange, nside)
    startx=vcat(range(xrange[1], xrange[2], length=nside), xrange[2]*ones(nside),
        reverse(range(xrange[1], xrange[2], length=nside)), ones(nside)*xrange[1])
    starty=vcat(ones(nside)*yrange[1], range(yrange[1], yrange[2], length=nside),
        ones(nside)*yrange[2], reverse(range(yrange[1], yrange[2], length=nside)))
    outplot=plot(0,0, alpha=0, xlabel="X₁", ylabel="X₂", size = (500, 400))
    for r in 1:length(startx)
        zzsk =  zz(; NS=1000, x0_0=[startx[r], starty[r]], tmax=tmax_tuned)
        plot!(zzsk["SK"][:, 2], zzsk["SK"][:, 3])
    end
    return outplot
end


function tailstart_hmc(; xrange, yrange, nside)
    startx=vcat(range(xrange[1], xrange[2], length=nside), xrange[2]*ones(nside),
        reverse(range(xrange[1], xrange[2], length=nside)), ones(nside)*xrange[1])
    starty=vcat(ones(nside)*yrange[1], range(yrange[1], yrange[2], length=nside),
        ones(nside)*yrange[2], reverse(range(yrange[1], yrange[2], length=nside)))
    outplot=plot(0,0, alpha=0, xlabel="X₁", ylabel="X₂", size = (500, 400))
    for r in 1:length(startx)
        hmc = runHMC(;epsilon=Lε_tuned/L_tuned,L=L_tuned,IT=1000,qs=[startx[r], starty[r]])
        plot!(hmc["SampleQ"][:, 1], hmc["SampleQ"][:, 2])
    end
    return outplot
end

function tailend_hmc(; xrange, yrange, nside)
    startx=vcat(range(xrange[1], xrange[2], length=nside), xrange[2]*ones(nside),
        reverse(range(xrange[1], xrange[2], length=nside)), ones(nside)*xrange[1])
    starty=vcat(ones(nside)*yrange[1], range(yrange[1], yrange[2], length=nside),
        ones(nside)*yrange[2], reverse(range(yrange[1], yrange[2], length=nside)))
    outplot=plot(startx,starty,seriestype=:scatter, xlabel="X₁", ylabel="X₂", size = (500, 400))
    return outplot
end

function bs_zz(;try_nb, ZZs = zzs, D=Dim)
    R=length(ZZs)
    zzESS_nb= Array{Float64,3}(undef, R, length(try_nb), D)
    for r in 1:R
        zz_r=ZZs[[r]]
        for i in 1:length(try_nb)
            for d in 1:D
                zzESS_nb[r, i, d] = zzsummaries(dms=d, sk=zz_r, B=try_nb[i])["EffectiveSampleSize"]
            end
        end
    end
    zzmeanESS_nb =  mapslices(mean, zzESS_nb, dims = [1])
    outplot=plot(string.((try_nb)), zzmeanESS_nb[1,:,:], title="ZZ - Average ESS per dimensions",
        xlabel="Number of batches", ylabel="Mean ESS", legend=true,
        labels=string.("dim", transpose(1:Dim)))
    return outplot
end


function bs_hmc(;try_nb, HMCs = hmcs, D=Dim)
    R=length(HMCs)
    hmcESS_nb= Array{Float64,3}(undef, R, length(try_nb), D)
    for r in 1:R
        hmc_r=HMCs[[r]]["SampleQ"]
        for i in 1:length(try_nb)
            for d in 1:D
                hmcESS_nb[r, i, d] = ESSbm(smpl=hmc_r[:, d], nbatches=try_nb[i])
            end
        end
    end
    hmcmeanESS_nb =  mapslices(mean, hmcESS_nb, dims = [1])
    outplot=plot(string.((try_nb)), hmcmeanESS_nb[1,:,:], title="HMC - Average ESS per dimensions",
        xlabel="Number of batches", ylabel="Mean ESS", legend=true,
        labels=string.("dim", transpose(1:Dim)))
    return outplot
end


function ESSsummaries(;ZZs=zzs, D=Dim, nbZZ=nbZZ_tuned,
     HMCs=hmcs, nbHMC=nbHMC_tuned)

     R=length(ZZs)
     zzESS_nb= Array{Float64,2}(undef, R, D)
     for r in 1:R
         zz_r=ZZs[[r]]
         for d in 1:D
             zzESS_nb[r, d] = zzsummaries(dms=d, sk=zz_r, B=nbZZ)["EffectiveSampleSize"]
         end
     end
     sumZZ1 = mapslices(mean, zzESS_nb, dims=1)
     sumZZ2 = mapslices(minimum, zzESS_nb, dims=1)
     sumZZ345 = mapslices( t -> quantile(t,[0.25, 0.50, 0.75]), zzESS_nb, dims=1)
     sumZZ6 = mapslices(maximum, zzESS_nb, dims=1)
     outzz=string(mapslices(summarystats, zzESS_nb, dims=1))

     hmcESS_nb= Array{Float64,2}(undef, R, D)
     for r in 1:R
         hmc_r=HMCs[[r]]["SampleQ"]
         for d in 1:D
             hmcESS_nb[r, d] = ESSbm(smpl=hmc_r[:, d], nbatches=nbHMC)
         end
     end
     sumHMC1 = mapslices(mean, hmcESS_nb, dims=1)
     sumHMC2 = mapslices(minimum, hmcESS_nb, dims=1)
     sumHMC345 = mapslices( t -> quantile(t,[0.25, 0.50, 0.75]), hmcESS_nb, dims=1)
     sumHMC6 = mapslices(maximum, hmcESS_nb, dims=1)
     outhmc=string(mapslices(summarystats, hmcESS_nb, dims=1))

     outmatrixpaper = zeros(6, D*2)
     for d in 1:D
         outmatrixpaper[1,(d-1)*2+1] = sumZZ1[1, d]
         outmatrixpaper[1,(d-1)*2+2] = sumHMC1[1, d]
         outmatrixpaper[2,(d-1)*2+1] = sumZZ2[1, d]
         outmatrixpaper[2,(d-1)*2+2] = sumHMC2[1, d]
         outmatrixpaper[3:5,(d-1)*2+1] = sumZZ345[:, d]
         outmatrixpaper[3:5,(d-1)*2+2] = sumHMC345[:, d]
         outmatrixpaper[6,(d-1)*2+1] = sumZZ6[1, d]
         outmatrixpaper[6,(d-1)*2+2] = sumHMC6[1, d]
     end

     outstring= latexify(round.(outmatrixpaper, digits=2))

     longstring=string("ZZ summaries \n ", outzz, "\n HMC summaries \n ", outhmc)

     output= Dict([("string", outstring), ("essZZ", zzESS_nb), ("essHMC", hmcESS_nb),
             ("longstring", longstring)])
     return output
end




function bs_tp_zz(;try_nb, ZZs = zzs, D=Dim)
    R=length(ZZs)
    zzESS_nb= Array{Float64,3}(undef, R, length(try_nb), D)
    for r in 1:R
        zz_r=ZZs[[r]]
        sm_r= zzsample(;N=1000000, sk=zz_r)
        for i in 1:length(try_nb)
            for d in 1:D
                zzESS_nb[r, i, d] = ESStailprob(smpl=sm_r[:,d], nbatches=try_nb[i])
            end
        end
    end
    zzmeanESS_nb =  mapslices(mean, zzESS_nb, dims = [1])
    outplot=plot(string.((try_nb)), zzmeanESS_nb[1,:,:], title="ZZ - Average ESS per dimensions",
        xlabel="Number of batches", ylabel="Mean ESS", legend=true,
        labels=string.("dim", transpose(1:Dim)))
    return outplot
end


function bs_tp_hmc(;try_nb, HMCs = hmcs, D=Dim)
    R=length(HMCs)
    hmcESS_nb= Array{Float64,3}(undef, R, length(try_nb), D)
    for r in 1:R
        hmc_r=HMCs[[r]]["SampleQ"]
        for i in 1:length(try_nb)
            for d in 1:D
                hmcESS_nb[r, i, d] = ESStailprob(smpl=hmc_r[:, d], nbatches=try_nb[i])
            end
        end
    end
    hmcmeanESS_nb =  mapslices(mean, hmcESS_nb, dims = [1])
    outplot=plot(string.((try_nb)), hmcmeanESS_nb[1,:,:], title="HMC - Average ESS per dimensions",
        xlabel="Number of batches", ylabel="Mean ESS", legend=true,
        labels=string.("dim", transpose(1:Dim)))
    return outplot
end
