
# Remember we've hacked the source code in solution_and_tracer_tendencies.jl to give us modified shallow water

using Oceananigans
using Oceananigans.Models: ShallowWaterModel
using MAT
using Printf

# ## Two-dimensional domain 

# The shallow water model is two-dimensional and uses grids that are `Flat`
# in the vertical direction. We use length scales non-dimensionalized by the width
# of the Bickley jet.

grid = RectilinearGrid(GPU(),size = (256, 256),
                       x = (0, 2π),
                       y = (0, 2π),
                       topology = (Periodic, Periodic, Flat))

# ## Building a `ShallowWaterModel`
#
# We build a `ShallowWaterModel` with the `WENO` advection scheme,
# 3rd-order Runge-Kutta time-stepping, non-dimensional Coriolis, and
# gravitational acceleration

Fr = 0.3
Ro = 0.4 # These match geostrophic balance of the ICs, so don't change them without generating new ICs

gravitational_acceleration = 1/Fr^2
coriolis = FPlane(f=1/Ro)

model = ShallowWaterModel(; grid, coriolis, gravitational_acceleration,
                          timestepper = :RungeKutta3,
                           momentum_advection = WENO())


f = coriolis.f
g = gravitational_acceleration

# Get initial condition from mat file
file = matopen("uvh_256_Fr_0_3_Ro_0_4_wave_0_5.mat")
u_i = read(file, "u") 
v_i = read(file, "v")
h_i = read(file, "h")
close(file)


uh_i = u_i.*h_i # dot is important for elementwise array multiplication!
vh_i = v_i.*h_i

set!(model, uh = uh_i, vh = vh_i, h= h_i )

## Build velocities and vorticity
uh, vh, h = model.solution
u = Field(uh / h) # This is fine when these things aren't arrays, otherwise need elementwise division
v = Field(vh / h)
ω = Field(∂x(v) - ∂y(u))
#dudt = Field(∂t(u))

## Running a simulation
simulation = Simulation(model, Δt = 1e-3, stop_time = 20)

function progress(sim)
    model = sim.model
    uh, vh, h = model.solution
    @info @sprintf("Simulation time: %s, max(|uh|, |vh|, |h|): %.2e, %.2e, %.2e \n", 
                   prettytime(sim.model.clock.time), 
                   maximum(abs, uh), maximum(abs, vh), 
                   maximum(abs, h))

    return nothing
end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))

# Build the `output_writer` for the two-dimensional fields to be output.
# Output every `t = 1.0`.

fields_filename= joinpath(@__DIR__, "SW_vort")

# We'll output both datasets, just because we have the movie working for the netcdf but want the jld2 for reading forcing
simulation.output_writers[:fields_jld2] = JLD2OutputWriter(model, (; ω,u,v),
                                                        filename = fields_filename,
                                                        schedule = TimeInterval(0.1),
                                                        overwrite_existing = true)

simulation.output_writers[:fields_nc] = NetCDFOutputWriter(model, (; ω,u,v ),
                                                        filename = fields_filename,
                                                        schedule = TimeInterval(0.1),
                                                        overwrite_existing = true)
# And finally run the simulation.
run!(simulation)
