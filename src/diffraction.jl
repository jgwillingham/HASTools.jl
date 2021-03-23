

mutable struct Diffraction
    exp::Experiment
    data::Array
    angles::Array
    counts::Array
    datafit
    specular::Float64
    function Diffraction(experiment, datapath)
        data = getdata(datapath, skipto=3)
        angles = data.Column1
        counts = data.Column2
        new(experiment, data, angles, counts)
    end
end



function fit!(diff::Diffraction, peakmodel=gaussian_model, p0::Array=[1e3,70.,2.], bgp0=[])
    fullmodel, fullp0 = getfullmodel(peakmodel, p0, bgp0)
    thefit = curve_fit(fullmodel, diff.angles, diff.counts, fullp0)

    num_bgparam = length(bgp0)
    bgparam = thefit.param[end-num_bgparam+1:end]
    datafit = DataFit(diff, thefit, fullmodel, peakmodel, bgparam)
    diff.datafit = datafit

    # Now set tallest peak to specular by default
    if peakmodel == pseudovoigt_model
        peakheights = [peak[1]+peak[4] for peak in datafit.peakparam]
    else
        peakheights = [peak[1] for peak in datafit.peakparam]
    end
    tallest_peak = datafit.peaks[argmax(peakheights)]
    setspecular!(datafit.data, tallest_peak)
    return datafit
end


function setspecular!(diff::Diffraction, peak)
    diff.specular = peak
end


function get_Δk_values(datafit)
    k0 = datafit.data.exp.k0
    specular = datafit.data.specular

    θin  =  (180. - specular)/2
    θout_list = [(180. - θin - peak) for peak in datafit.peaks]

    Δk_list = [k0 * (sind(θout) - sind(θin)) for θout in θout_list]
    return Δk_list
end


function peakanalysis(datafit)
    peaks = datafit.peaks
    Δk_list = get_Δk_values(datafit)
    b1, b2, = datafit.data.exp.crystal.meshReciprocals
    G₁₀ = √(b1[1]^2 + b1[2]^2 + b1[3]^2)  # take norm without using LinearAlgebra.jl
    G₀₁ = √(b2[1]^2 + b2[2]^2 + b2[3]^2)
    G₁₁ = √((b1+b2)[1]^2 + (b1+b2)[2]^2 + (b1+b2)[3]^2)
    ratios₁₀ = round.([Δk/G₁₀ for Δk in Δk_list], digits=2)
    ratios₀₁ = round.([Δk/G₀₁ for Δk in Δk_list], digits=2)
    ratios₁₁ = round.([Δk/G₁₁ for Δk in Δk_list], digits=2)
    # Print results in pretty way
    println("Direction\t       Δk/G")
    println("  (1,0)\t\t", ratios₁₀)
    println("\n  (0,1)\t\t", ratios₀₁)
    println("\n  (1,1)\t\t", ratios₁₁)
end


function get_expected_peaks(datafit, direction::Tuple)
    b1, b2, = datafit.data.exp.crystal.meshReciprocals
    k0 = datafit.data.exp.k0
    n1, n2 = direction
    G = n1*b1 + n2*b2
    Gnorm = √(G[1]^2 + G[2]^2 + G[3]^2)
    θin = (180. - datafit.data.specular)/2
    expected_peaks = []
    for n in -3:3
        try
            θout = asind( sind(θin)+n*Gnorm/k0 ) # from momentum conservation mod recip. vector
            peak = 180. - θin - θout # peak position as detector angle
            push!(expected_peaks, peak)
        catch DomainError
            continue
        end
    end
    exp_peaks = round.(expected_peaks, digits=1)
    return exp_peaks
end
