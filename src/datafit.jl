



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

    num_bgparam = length(bgparam)
    num_peakparam = length(fit.param) - num_bgparam
    peakparam = [fit.param[i:i+ppp-1] for i in 1:ppp:num_peakparam]
    peaks = [peak[2] for peak in peakparam]
    return peaks, peakparam
end



function plotdata(diff::Diffraction)
    plot(diff.angles, diff.counts, label=nothing, c=:black, ms=3.5,
        markershape=:star4, ls=:dash, linealpha=0.4, size=(400,250))
    xlabel!("Detector Angle (°)")
    ylabel!("Counts")
    xlims!(min(diff.angles...), max(diff.angles...))
end


function plotdata(tof::TOF)
    plot(tof.times, tof.counts, label=nothing, c=:black, ms=3.5,
        markershape=:star4, ls=:dash, linealpha=0.4, size=(400,250))
    xlabel!("t (μs)")
    ylabel!("Counts")
    xlims!(min(tof.times...), max(tof.times...))
end



function plotfit(datafit::DataFit{Diffraction}; withpeaks=false, withribbon=true)
    smooth_angles = [a for a in min(datafit.data.angles...):0.01:max(datafit.data.angles...)]
    bg = background(smooth_angles, datafit.bgparam)
    fitcurve = datafit.model(smooth_angles, datafit.param)

    if withribbon
        stdev = sqrt( sum(datafit.resid.^2) / length(datafit.resid) )
        ribbon = stdev
    else
        ribbon = nothing
    end
    plt = plot(smooth_angles, fitcurve, lw=2,c=:darkred, label=nothing,
        ribbon=ribbon, fillalpha=0.3, fillcolor=:orange, foreground_color_legend=nothing)
    plot!(datafit.data.angles, datafit.data.counts, label=nothing,
        c=:black, alpha=0.8, ms=3.5, markershape=:star4, ls=:dash, linealpha=0.4)
    ymin, ymax = ylims(plt)

    if withpeaks
        for peakinfo in datafit.peakparam
            peak = datafit.peakmodel(smooth_angles, peakinfo)
            peakangle = round(peakinfo[2], digits=1)
            plot!(smooth_angles, peak.+bg, fillrange=[bg peak.+bg], alpha=0.35,
                label="$peakangle"*"°", palette=:seaborn_dark)
        end
    end

    xlabel!("Detector Angle (°)")
    ylabel!("Counts")
    xlims!(min(datafit.data.angles...), max(datafit.data.angles...))
    ylims!(ymin, ymax)
end



function plotfit(datafit::DataFit{TOF}; withpeaks=false, withribbon=true)
    smoothtimes = [t for t in min(datafit.data.times...):0.01:max(datafit.data.times...)]
    bg = background(smoothtimes, datafit.bgparam)
    fitcurve = datafit.model(smoothtimes, datafit.param)

    if withribbon
        stdev = sqrt( sum(datafit.resid.^2) / length(datafit.resid) )
        ribbon = stdev
    else
        ribbon = nothing
    end
    plt = plot(smoothtimes, fitcurve, lw=2,c=:darkred, label=nothing,
        ribbon=ribbon, fillalpha=0.3, fillcolor=:orange, foreground_color_legend=nothing)
    plot!(datafit.data.times, datafit.data.counts, label=nothing, c=:black, alpha=0.8, ms=3.5,
        markershape=:star4, ls=:dash, linealpha=0.4)
    ymin, ymax = ylims(plt)

    if withpeaks
        for peakinfo in datafit.peakparam
            peak = datafit.peakmodel(smoothtimes, peakinfo)
            peakangle = round(peakinfo[2], digits=1)
            plot!(smoothtimes, peak.+bg, fillrange=[bg peak.+bg], alpha=0.35,
                label="$peakangle μs", palette=:seaborn_dark)
        end
    end
    xlabel!("t (μs)")
    ylabel!("Counts")
    xlims!(min(datafit.data.times...), max(datafit.data.times...))
    ylims!(ymin, ymax)
end
