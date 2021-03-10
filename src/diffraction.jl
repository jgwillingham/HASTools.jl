

mutable struct Diffraction
    data::Array
    angles::Array
    counts::Array
    fitparams::Array
    specular::Float64
    function Diffraction(datapath)
        data = getdata(datapath)
        angles = data[!,"angles"]
        counts = data[!,"counts"]
        new(data, angles, counts)
    end
end
