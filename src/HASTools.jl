
module HASTools

using Plots
using LaTeXStrings
using CSV
using LsqFit
using Dates

include("models.jl")
include("experiment.jl")
include("diffraction.jl")
include("tof.jl")
include("datafit.jl")


export gaussian_model, lorentzian_model, pseudovoigt_model,
    gaussians, lorentzians, pseudovoigts, background,
    Experiment, TOF, Diffraction,
    fit, plotfit,
    plotdata

end
