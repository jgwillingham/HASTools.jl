


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



function fit(tof::TOF, p0::Array, model=gaussian_model, bgp0=nothing)
    fullmodel, full_p0 = getfullmodel(model, p0, bgp0)
    thefit = curve_fit(fullmodel, tof.times, tof.counts, full_p0)
    return thefit, fullmodel
end


function plotfit(tof::TOF, toffit)
    thefit, model = toffit
    plot(tof.times, tof.counts, label=nothing, foreground_color_legend=nothing)
    smoothtimes = [t for t in min(tof.times...):0.01:max(tof.times...)]
    plot!(smoothtimes, model(smoothtimes, thefit.param), lw=2)
    xlabel!("t (μs)")
    ylabel!("Counts")
end
