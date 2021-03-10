


mutable struct TOF
    data::Array
    times::Array
    counts::Array
    function TOF(datapath)
        data = getdata(datapath)
        times = data[!,"times"]
        counts = data[!,"counts"]
        new(data, times, counts)
    end
end
