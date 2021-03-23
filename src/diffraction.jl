

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



function fit!(diff::Diffraction, peakmodel=gaussian_model, p0::Array=[1e3,70.,2.], bgp0=nothing)
    fullmodel, fullp0 = getfullmodel(peakmodel, p0, bgp0)
    thefit = curve_fit(fullmodel, diff.angles, diff.counts, fullp0)
    if bgp0 == nothing
        num_bgparam = 0
    else
        num_bgparam = length(bgp0)
    end
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

    θin  = π/180 * (180. - specular)/2
    θout_list = [π/180 * (180. - specular - peak) for peak in datafit.peaks]

    Δk_list = [k0 * (sin(θin) - sin(θout)) for θout in θout_list]
    return Δk_list
end


function peakanalysis(datafit)
    peaks = datafit.peaks
    Δk_list = get_Δk_values(datafit)
    return Δk_list
end
