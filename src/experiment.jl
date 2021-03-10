


mutable struct Experiment
    crystal::String
    Tnoz::Real
    date::Tuple
    specular::Float64
    k0::Float64
    E0::Float64
    diffractiondata::Array
    tofdata::Array
    function Experiment(crystal, Tnoz, date=today())
        kB = 8.617e-2 # Boltzmann's constant in meV/K
        E0 = (5/2)*kB*Tnoz  # 5/2 comes from forward propulsion at the nozzle
        k0 = 0.642*√Tnoz # from just solving E0 = ħ^2 k0^2 / 2m (in inverse angstrom)
        specular = nothing
        new(crystal, Tnoz, date, specular, k0, E0)
    end
end
