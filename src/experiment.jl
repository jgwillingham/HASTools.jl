


mutable struct Experiment
    crystal # intended to be of type 'Slab' from LatticeDynamics.jl
    Tnoz::Real
    date::Date
    k0::Float64
    E0::Float64
    Tsample::Real
    diffractiondata::Array
    tofdata::Array
    function Experiment(crystal, Tnoz, date=today())
        kB = 8.617e-2 # Boltzmann's constant in meV/K
        E0 = (5/2)*kB*Tnoz  # 5/2 comes from forward propulsion at the nozzle
        k0 = 0.642*√Tnoz # just k₀ =√(2mE₀)/ħ (in inverse angstroms)
        new(crystal, Tnoz, date, k0, E0)
    end
end




function getdata(datapath::String; skipto::Int=1)
    data = CSV.File(datapath, skipto=skipto, header=0, normalizenames=true);
    return data
end



function getΔE(exp::Experiment, Δk::Real, angle::Real, specular::Real)
    θin = (180. - specular)/2
    θout = 180 - θin - angle # `angle` is a detector angle
    ΔE = exp.E0 * ( (( sind(θin) + Δk/exp.k0 ) / sind(θout))^2 - 1 )
    return ΔE
end


function scancurve(exp::Experiment, Δklist::Array, angle::Real, specular::Real)
    energy(Δk) = getΔE(exp, Δk, angle, specular)
    Elist = energy.(Δklist)
    return Elist
end


function plotscancurves(datafit, startangle::Real, stopangle::Real, direction::Tuple; n::Int=0)
    specular = datafit.data.specular
    exp = datafit.data.exp
    b1, b2, = exp.crystal.meshReciprocals
    n1, n2 = direction
    G = n1*b1 + n2*b2
    Gnorm = √(G[1]^2 + G[2]^2 + G[3]^2)
    Δklist = [Δk for Δk in 0:0.01:Gnorm/2] # half of the Brillouin zone along this direction
    println("Scan curves for order-$n peak:\n(Calculated with theoretical peak position)")
    plt = plot(legend=false, gridalpha=0.6, minorgrid=true, minorgridalpha=0.2)

    startcurve = scancurve(exp, Δklist .+ n*Gnorm, startangle, specular)
    stopcurve = scancurve(exp, Δklist .+ n*Gnorm, stopangle, specular)
    lowcurve, highcurve = ordercurves(startcurve, stopcurve)
    plot!(Δklist, highcurve, fillrange=[lowcurve highcurve], c=:lightgreen, fillcolor=:lightgreen, alpha=0.5)
    plot!(Δklist, abs.(highcurve), fillrange=[abs.(lowcurve) abs.(highcurve)], c=:lightgreen, fillcolor=:lightgreen, alpha=0.5)

    # capture all kinematically allowed pairs in other half of BZ
    kflipped_startcurve = scancurve(exp, -Δklist .+ n*Gnorm, startangle, specular)
    kflipped_stopcurve = scancurve(exp, -Δklist .+ n*Gnorm, stopangle, specular)
    lowcurve, highcurve = ordercurves(kflipped_startcurve, kflipped_stopcurve)
    plot!(Δklist, highcurve, fillrange=[lowcurve highcurve], c=:lightgreen, fillcolor=:lightgreen, alpha=0.5)
    plot!(Δklist, abs.(highcurve), fillrange=[abs.(lowcurve) abs.(highcurve)], c=:lightgreen, fillcolor=:lightgreen, alpha=0.5)
    plot!([0.], seriestype=:hline, lw=2, c=:black)
    xlabel!("Δk (Å⁻¹)")
    ylabel!("ΔE (meV)")
    xlims!(0, Gnorm/2)
    ylims!(0, Inf)
end


function ordercurves(curve1, curve2)
    lowcurve, highcurve = [], []
    for i in eachindex(curve1)
        if curve1[i] > curve2[i]
            push!(highcurve, curve1[i])
            push!(lowcurve, curve2[i])
        elseif curve1[i] < curve2[i]
            push!(highcurve, curve2[i])
            push!(lowcurve, curve1[i])
        else
            push!(highcurve, curve1[i])
            push!(lowcurve, 0.)
        end
    end
    return lowcurve, highcurve
end
