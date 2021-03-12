


mutable struct TOF
    data::Array
    times::Array
    counts::Array
    function TOF(datapath)
        data = getdata(datapath, skipto=7)
        times = data.Column1
        counts = data.Column2
        new(data, times, counts)
    end
end
