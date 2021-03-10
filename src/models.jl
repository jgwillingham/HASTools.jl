


gaussian(x, p...) = p[1] * exp(-(x-p[2])^2/(2*p[3]^2) )
lorentzian(x, p...) = p[1]*(1/π)*p[3] /((x-p[2])^2 + p[3]^2)
pseudovoigt(x, p...) = gaussian(x, p[1:3]...) + lorentzian(x, [p[4], p[2], p[3]]...)

# p0 contains triples A,μ,σ for each Gaussian / Lorentzian in the model 
function gaussian_model(X::Array, p0::Array)
    out = zeros(length(X))
    for i in 1:3:length(p0)
        out += gaussian.(X, p0[i:i+2]...)
    end
    return out
end


function lorentzian_model(X::Array, p0::Array)
    out = zeros(length(X))
    for i in 1:3:length(p0)
        out += lorentzian.(X, p0[i:i+2]...)
    end
    return out
end


function pseudovoigt_model(X::Array, p0::Array)
    out = zeros(length(X))
    for i in 1:4:length(p0)
        out += pseudovoigt.(X,p0[i:i+3]...)
    end
    return out
end


function gaussians(X::Array, p::Array)
    gauss = [gaussian.(X, p[i:i+2]...) for i in 1:3:length(p)]
    return gauss
end


function lorentzians(X::Array, p::Array)
    lorentzes = [lorentzian.(X, p[i:i+2]...) for i in 1:3:length(p)]
    return lorentzes
end


function fit_tof(X::Array, counts::Array, p0::Array, model="lorentzian")
    if model == "lorentzian"
        thefit = curve_fit(lorentzian_model, X, counts, p0)
    elseif model == "gaussian"
        thefit = curve_fit(gaussian_model, X, counts, p0)
    elseif model == "pseudovoigt"
        thefit = curve_fit(pseudovoigt_model, X, counts, p0)
    end
    return thefit
end
