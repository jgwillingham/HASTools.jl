



struct DataFit{T<:Union{Diffraction, TOF}}
    data::T
    model
    peakmodel
    param::Array
    resid::Array
    peakparam::Array
    bgparam::Array
    peaks::Array
    function DataFit(data, fit, fullmodel, peakmodel, bgparam)
        param = fit.param
        resid = fit.resid
        peaks, peakparam = getPeaks(fit, peakmodel, bgparam)
        new{typeof(data)}(data, fullmodel, peakmodel, param, resid, peakparam, bgparam, peaks)
    end
end



function getPeaks(fit, peakmodel, bgparam)
    if peakmodel == gaussian_model || peakmodel == lorentzian_model
        ppp = 3  # ppp = parameters per peak
    elseif peakmodel == pseudovoigt_model
        ppp = 4
    end
    if bgparam == nothing
        num_bgparam = 0
    else
        num_bgparam = length(bgparam)
    end
    num_peakparam = length(fit.param) - num_bgparam
    peakparam = [fit.param[i:i+ppp-1] for i in 1:ppp:num_peakparam]
    peaks = [peak[2] for peak in peakparam]
    return peaks, peakparam
end




function plotfit(datafit::DataFit{Diffraction})
    smooth_angles = [a for a in min(datafit.data.angles...):0.01:max(datafit.data.angles...)]

    stdev = sqrt( sum(datafit.resid.^2) / length(datafit.resid) )
    plot(smooth_angles, datafit.model(smooth_angles, datafit.param), lw=2,c=:darkred, legend=false,
        ribbon=stdev, fillalpha=0.5, fillcolor=:orange)
    scatter!(datafit.data.angles, datafit.data.counts, label=nothing, c=:black, alpha=0.8, ms=3.5,
        markershape=:star4)

    xlabel!("Angle (ᵒ)")
    ylabel!("Counts")
    xlims!(min(datafit.data.angles...), max(datafit.data.angles...))
end



function plotfit(datafit::DataFit{TOF})
    smoothtimes = [t for t in min(datafit.data.times...):0.01:max(datafit.data.times...)]

    stdev = sqrt( sum(datafit.resid.^2) / length(datafit.resid) )
    plot(smoothtimes, datafit.model(smoothtimes, datafit.param), lw=2,c=:darkred, legend=false,
        ribbon=stdev, fillalpha=0.5, fillcolor=:orange)
    scatter!(datafit.data.times, datafit.data.counts, label=nothing, c=:black, alpha=0.8, ms=3.5,
        markershape=:star4)

    xlabel!("t (μs)")
    ylabel!("Counts")
    xlims!(min(datafit.data.times...), max(datafit.data.times...))
end
