using Oceananigans
using Oceananigans.Operators
using Oceananigans.OutputReaders
using Oceananigans.OutputReaders: OnDisk
using Oceananigans.Units
using Oceananigans.Utils: Time 
using Printf
using CairoMakie


print("Hello")
arch = GPU()
grid = RectilinearGrid(arch, size=(2, 2, 2), extent=(1, 1, 1))
times = 0:0.1:3
filename = "MWE_file.jld2"

u_forcing_tmp = Field{Center,Center,Center}(grid) 
u_forcing = FieldTimeSeries{Face, Center, Center}(grid, times; backend = OnDisk(), path = filename, name = "u_forcing")

for (t, time) in enumerate(u_forcing.times)
    @info "writing down data for timestep $t and time $time"
    set!(u_forcing_tmp, (x, y,z) -> sin(π * x) * time)
    set!(u_forcing, u_forcing_tmp, t)

    #set!(u_forcing[t], (x, y, z) -> sin(π * x) * time)
end

u_forcing = FieldTimeSeries(filename, "u_forcing"; backend = InMemory(4))
print("About to define model")
model = NonhydrostaticModel(; grid, forcing=(; u=u_forcing))
time_step!(model, 1)

# # Make sure the field time series updates correctly
# u_forcing = FieldTimeSeries{Face, Center, Center}(grid, 0:0.1:4; backend = InMemory(2))

# model = NonhydrostaticModel(; grid, forcing=(; u=u_forcing))
# time_step!(model, 2)
# time_step!(model, 2)