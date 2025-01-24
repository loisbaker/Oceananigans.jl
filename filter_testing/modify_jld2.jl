using JLD2
using JLD2: Group
T = 10
fields_filename = joinpath(@__DIR__, "SW_vort_copy.jld2")
file = jldopen(fields_filename,"r")
iterations = parse.(Int, keys(file["timeseries/t"]))
times = [file["timeseries/t/$iter"] for iter in iterations]
close(file)

file = jldopen(fields_filename,"r+")
Base.delete!(file, "timeseries/t")
g = Group(file, "timeseries/t")
for (i,iter) in enumerate(iterations)
    g["$iter"] = T - times[i]
end
close(file)



