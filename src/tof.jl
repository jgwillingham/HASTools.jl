


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
    fit, model = toffit

    smoothtimes = [t for t in min(tof.times...):0.01:max(tof.times...)]

    stdev = sqrt( sum(fit.resid.^2) / length(fit.resid) )
    plot(smoothtimes, model(smoothtimes, fit.param), lw=2,c=:darkred, legend=false,
        ribbon=stdev, fillalpha=0.5, fillcolor=:orange)
    scatter!(tof.times, tof.counts, label=nothing, c=:black, alpha=0.8, ms=3.5,
        markershape=:star4)

    xlabel!("t (μs)")
    ylabel!("Counts")
    xlims!(min(tof.times...), max(tof.times...))
end
