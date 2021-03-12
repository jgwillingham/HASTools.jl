

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
