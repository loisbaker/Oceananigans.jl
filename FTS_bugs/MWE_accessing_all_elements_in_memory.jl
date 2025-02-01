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
# grid = RectilinearGrid(arch, size=(2, 2, 2), extent=(1, 1, 1))
# times = 0:0.1:3
filename = joinpath(@__DIR__, "MWE_forcing_file.jld2")

#generate_forcing_data!(grid, times, filename)
N_in_memory = 5
forcing_fts = FieldTimeSeries(filename, "forcing"; backend = InMemory())
model = NonhydrostaticModel(; grid, tracers=(:c), forcing=(; c=forcing_fts))

