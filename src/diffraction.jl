

mutable struct Diffraction
    data::Array
    angles::Array
    counts::Array
    fitparam::Array
    specular::Float64
    function Diffraction(datapath)
        data = getdata(datapath, skipto=3)
        angles = data.Column1
        counts = data.Column2
        new(data, angles, counts)
    end
end



function fit(diff::Diffraction, p0::Array, peakmodel=gaussian_model, bgp0=nothing)
    fullmodel, fullp0 = getfullmodel(peakmodel, p0, bgp0)
    thefit = curve_fit(fullmodel, diff.angles, diff.counts, fullp0)
    if bgp0 == nothing
        num_bgparam = 0
    else
        num_bgparam = length(bgp0)
    end
    bgparam = thefit.param[end-num_bgparam:end]
    datafit = DataFit(diff, thefit, fullmodel, peakmodel, bgparam)
    return datafit
end
