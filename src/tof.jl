


mutable struct TOF
    exp::Experiment
    data::Array
    times::Array
    counts::Array
    datafit
    function TOF(experiment, datapath)
        data = getdata(datapath, skipto=7)
        times = data.Column1
        counts = data.Column2
        new(experiment, data, times, counts)
    end
end




function fit!(tof::TOF, peakmodel=gaussian_model, p0::Array=[1e3,60.,3.], bgp0=[])
    fullmodel, fullp0 = getfullmodel(peakmodel, p0, bgp0)
    thefit = curve_fit(fullmodel, tof.times, tof.counts, fullp0)

    num_bgparam = length(bgp0)
    bgparam = thefit.param[end-num_bgparam+1:end]
    datafit = DataFit(tof, thefit, fullmodel, peakmodel, bgparam)
    tof.datafit = datafit
    return datafit
end
