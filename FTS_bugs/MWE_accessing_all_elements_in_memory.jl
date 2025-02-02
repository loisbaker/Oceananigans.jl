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
using Oceananigans.Fields: index_binary_search
using Oceananigans.OutputReaders: interpolating_time_indices
using Printf

function generate_forcing_data!(grid, times, filename)

    forcing_tmp = Field{Center,Center,Center}(grid) 
    forcing_FTS = FieldTimeSeries{Center, Center, Center}(grid, times; backend = OnDisk(), path = filename, name = "forcing")

    for (it, time) in enumerate(forcing_FTS.times)
        @info "writing down data for timestep $it and time $time"
        set!(forcing_tmp, (x, y, z) -> sin(2* π * time) )
        set!(forcing_FTS,forcing_tmp,it) 
    end
end

arch = GPU()
grid = RectilinearGrid(arch, size=(2, 2, 2), extent=(1, 1, 1))
times = 0:0.1:3.1
filename = joinpath(@__DIR__, "MWE_forcing_file.jld2")

#generate_forcing_data!(grid, times, filename)
forcing_fts = FieldTimeSeries(filename, "forcing"; backend = InMemory(4))
println(forcing_fts)

for i in 0:0.01:3.2
    println("")
    println("Next")
    println("Time = $i")
    v = maximum(forcing_fts[Time(i)])
    println(v)
    
    println(forcing_fts.backend)
end
# This generates weird results for backend InMemory(n)

# Let's just check that we don't error at the end of the interval because we're trying to load in too many indices 
# We should actually try to fix this error (only load to the end of the dataset), as it's the cause of the other error too, but might not occur naturally. 
# But we've fixed the weirdness hopefully by moving the order in field_tim_series indexing and going for the n2-1 index rather than the n2 so
# that we can still interpolate the gap. 

# Then check what's actually erroring in the update_velocities! in the perform_filtering 