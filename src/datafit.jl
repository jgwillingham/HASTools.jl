



struct DataFit{T<:Union{Diffraction, TOF}}
    data::T
    model
    param::Array
    resid::Array
    peakparam::Array
    bgparam::Array
    peaks::Array
    function DataFit(data, fit, model, bgparam)
        param = fit.param
        resid = fit.resid
        peaks, peakparam = getPeaks(fit, model, bgparam)
        new{typeof(data)}(data, model, param, resid, peakparam, bgparam, peaks)
    end
end



function getPeaks(fit, model, bgparam)
    if model == gaussian_model || model == lorentzian_model
        ppp = 3  # ppp = parameters per peak
    elseif model == pseudovoigt_model
        ppp = 4
    end
    numbgparam = length(bgparam)
    numpeakparam = length(fit.param) - numbgparam
    peakparam = [fit.param[i:i+ppp-1] for i in 1:ppp:numpeakparam]
    peaks = [peak[2] for peak in peakparam]
    return peaks, peakparam
end
