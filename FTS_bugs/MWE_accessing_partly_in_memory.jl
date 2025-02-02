# We want to make a MWE to demonstrate the bugs in FieldTImeSeries
# So far, we think that there's an error trying to read forcing up to the end point, 
# an error accessing the final element the first time, and an error in forcing
# with in partly in InMemory

using Oceananigans
using Oceananigans.Operators
using Oceananigans.OutputReaders
using Oceananigans.OutputReaders: OnDisk
using Oceananigans.Units
using Oceananigans.Utils: Time 
using Printf

function generate_forcing_data!(grid, times, filename)

    forcing_tmp = Field{Center,Center,Center}(grid) 
    forcing_FTS = FieldTimeSeries{Center, Center, Center}(grid, times; backend = OnDisk(), path = filename, name = "forcing")

    for (it, time) in enumerate(forcing_FTS.times)
        @info "writing down data for timestep $it and time $time"
        set!(forcing_tmp, (x, y, z) -> sin(2* π * time) )
        set!(forcing_FTS,forcing_tmp,it) #set!(fts::InMemoryFTS, value, n::Int) = set!(fts[n], value)
    end
end

arch = GPU()
grid = RectilinearGrid(arch, size=(2, 2, 2), extent=(1, 1, 1))
times = 0:0.1:3
filename = joinpath(@__DIR__, "MWE_forcing_file.jld2")

#generate_forcing_data!(grid, times, filename)
N_in_memory = 5
forcing_fts = FieldTimeSeries(filename, "forcing"; backend = InMemory(N_in_memory))


println("Here's our fts: $forcing_fts")
# ERROR NUMBER 1
nt = length(forcing_fts.times)

println("")
println("We can't access time nt - N_in_memory + 2:")
try 
    println(forcing_fts[nt - N_in_memory + 2])
catch e
    println("Accessing element $nt threw an error:")
    println(e)
end
println("")
println("But we can access time nt - N_in_memory+1:")
println(forcing_fts[nt - N_in_memory+1])

println("Now we try nt - N_in_memory + 2  again - no error")
try 
    println(forcing_fts[nt - N_in_memory + 2])
catch e
    println("Accessing element $nt threw an error:")
    println(e)
end


println("Here's our fts: $forcing_fts")
# u_forcing = FieldTimeSeries(filename, "u_forcing"; backend = InMemory(4))
# print("About to define model")
# model = NonhydrostaticModel(; grid, forcing=(; u=u_forcing))
# time_step!(model, 1)

# Current fts status given by backend: InMemory(a,b), where b is number in memory, and a is maybe the start?
# For InMemory(5) and length 31, if we request [28] we get an error, then request 27 we get no error and backend (27,5)
# then request 28 we get backend (27,5) again, but it knows 28. 

# This might not be a problem if going forwards, as the final few will be in memory. Problem comes from requesting [28],
# so it tries to load in 28,29,30,31,32 and 32 doesn't exist. Can load 27,28,29,30,31 though. 


# With fix:
