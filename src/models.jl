


#### PEAK MODELS ####

gaussian(x, p...) = p[1] * exp(-(x-p[2])^2/(2*p[3]^2) )
lorentzian(x, p...) = p[1]*p[3]^2 /((x-p[2])^2 + p[3]^2)
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

# p0 contains A,μ,σ,α for each pseudovoigt in the model
function pseudovoigt_model(X::Array, p0::Array)
    out = zeros(length(X))
    for i in 1:4:length(p0)
        out += pseudovoigt.(X,p0[i:i+3]...)
    end
    return out
end



#### BACKGROUND MODELS ####

function constant(X::Array, p::Array)
    out = [p[1] for x in X]
    return out
end


function linear(X::Array, p::Array)
    out = p[2].*X .+ p[1]
    return out
end


function quadratic(X::Array, p::Array)
    out = p[3].*X.^2 .+ p[2].*X .+ p[1]
    return out
end


function background(X::Array, p::Array)
    num_bgparam = length(p)
    if num_bgparam == 0
        bg = zeros(length(X))
    elseif num_bgparam == 1
        bg = constant(X, p)
    elseif num_bgparam == 2
        bg = linear(X, p)
    elseif num_bgparam == 3
        bg = quadratic(X, p)
    end
    return bg
end

#### COMBINE PEAK MODEL AND BACKGROUND MODEL ####

function getfullmodel(peakmodel, p0::Array, bgp0=[])
    num_bgparam = length(bgp0)
    if num_bgparam == 0
        bg(X, p) = constant(X, [0])
    elseif num_bgparam == 1
        bg = constant
    elseif num_bgparam == 2
        bg = linear
    elseif num_bgparam == 3
        bg = quadratic
    end
    fullp0 = append!(p0, bgp0)
    fullmodel(X, p) = peakmodel(X, p[1:end-num_bgparam]) + bg(X, p[end-num_bgparam+1:end])

    return fullmodel, fullp0
end




# Collection of peaks within a model

function gaussians(X::Array, p::Array)
    gauss = [gaussian.(X, p[i:i+2]...) for i in 1:3:length(p)]
    return gauss
end


function lorentzians(X::Array, p::Array)
    lorentzes = [lorentzian.(X, p[i:i+2]...) for i in 1:3:length(p)]
    return lorentzes
end


function pseudovoigts(X::Array, p::Array)
    pvs = [pseudovoigt.(X, p[i:i+3]...) for i in 1:4:length(p)]
    return pvs
end;
