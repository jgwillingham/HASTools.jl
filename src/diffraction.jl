

mutable struct Diffraction
    data::Array
    angles::Array
    counts::Array
    fitparams::Array
    specular::Float64
    function Diffraction(datapath)
        data = getdata(datapath, skipto=3)
        angles = data.Column1
        counts = data.Column2
        new(data, angles, counts)
    end
end


function fit(diff::Diffraction, p0::Array, model=gaussian_model, bgp0=nothing)
    fullmodel, full_p0 = getfullmodel(model, p0, bgp0)
    thefit = curve_fit(fullmodel, diff.angles, diff.counts, full_p0)
    return thefit, fullmodel
end


function plotfit(diff::Diffraction, diff_fit)
    fit, model = diff_fit

    smooth_angles = [a for a in min(diff.angles...):0.01:max(diff.angles...)]

    stdev = sqrt( sum(fit.resid.^2) / length(fit.resid) )
    plot(smooth_angles, model(smooth_angles, fit.param), lw=2,c=:darkred, legend=false,
        ribbon=stdev, fillalpha=0.5, fillcolor=:orange)
    scatter!(diff.angles, diff.counts, label=nothing, c=:black, alpha=0.8, ms=3.5,
        markershape=:star4)
        
    xlabel!("Angle (ᵒ)")
    ylabel!("Counts")
    xlims!(min(diff.angles...), max(diff.angles...))
end
