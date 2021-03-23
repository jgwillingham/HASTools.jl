


mutable struct Experiment
    crystal # intended to be of type 'Slab' from LatticeDynamics.jl
    Tnoz::Real
    date::Date
    k0::Float64
    E0::Float64
    Tsample::Real
    diffractiondata::Array
    tofdata::Array
    function Experiment(crystal, Tnoz, date=today())
        kB = 8.617e-2 # Boltzmann's constant in meV/K
        E0 = (5/2)*kB*Tnoz  # 5/2 comes from forward propulsion at the nozzle
        k0 = 0.642*√Tnoz # just k₀ =√(2mE₀)/ħ (in inverse angstroms)
        new(crystal, Tnoz, date, k0, E0)
    end
end




function getdata(datapath::String; skipto::Int=1)
    data = CSV.File(datapath, skipto=skipto, header=0, normalizenames=true);
    return data
end
