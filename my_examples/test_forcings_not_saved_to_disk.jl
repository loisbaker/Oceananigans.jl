# include("../test/dependencies_for_runtests.jl")
using Oceananigans
using Oceananigans.Operators
using Oceananigans.OutputReaders
using Oceananigans.OutputReaders: OnDisk
using Oceananigans.Units
using Oceananigans.Utils: Time 
using Printf
using CairoMakie
using Oceananigans.BoundaryConditions: ImpenetrableBoundaryCondition
using Oceananigans.Fields: Field
using Oceananigans.Forcings: MultipleForcings

print("Hello")
arch = GPU()
grid = RectilinearGrid(arch, size=(2, 2, 2), extent=(1, 1, 1))

u_forcing = FieldTimeSeries{Face, Center, Center}(grid, 0:0.1:3)

for (t, time) in enumerate(u_forcing.times)
    set!(u_forcing[t], (x, y, z) -> sin(π * x) * time)
end

model = NonhydrostaticModel(; grid, forcing=(; u=u_forcing))
time_step!(model, 1)

# Make sure the field time series updates correctly
u_forcing = FieldTimeSeries{Face, Center, Center}(grid, 0:0.1:4; backend = InMemory(2))

model = NonhydrostaticModel(; grid, forcing=(; u=u_forcing))
time_step!(model, 2)
time_step!(model, 2)

    # @test u_forcing.backend.start == 4

    # return true
