


mutable struct TOF
    data::Array
    times::Array
    counts::Array
    function TOF(datapath)
        data = getdata(datapath, skipto=7)
        times = data.Column1
        counts = data.Column2
        new(data, times, counts)
    end
end




function fit(tof::TOF, p0::Array, peakmodel=gaussian_model, bgp0=nothing)
    fullmodel, fullp0 = getfullmodel(peakmodel, p0, bgp0)
    thefit = curve_fit(fullmodel, tof.times, tof.counts, fullp0)
    if bgp0 == nothing
        num_bgparam = 0
    else
        num_bgparam = length(bgp0)
    end
    bgparam = thefit.param[end-num_bgparam+1:end]
    datafit = DataFit(tof, thefit, fullmodel, peakmodel, bgparam)
    return datafit
end
