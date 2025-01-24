# Remember we've hacked the source code in solution_and_tracer_tendencies.jl to give us modified shallow water

using Oceananigans
using Oceananigans.Models: HydrostaticFreeSurfaceModel
using Printf
using Oceananigans.Units: Time
using Oceananigans.Units

arch = GPU()
Nx, Ny, Nz = 10, 10, 10

grid = RectilinearGrid(arch, size = (Nx, Ny),
                        x = (0, 2π),
                        y = (0, 2π),
                        halo = (4, 4),
                        topology = (Periodic, Periodic, Flat))

tracers= (:c)

# Test
forcing_func(x,y,t,c) = sin(x) + c

c_forcing = Forcing(forcing_func,field_dependencies = :c)
forcing = (;c=c_forcing)
velocities = PrescribedVelocityFields(u=0,v=0)
model = HydrostaticFreeSurfaceModel(;grid, tracers= tracers, velocities =velocities,forcing=forcing)

## Running a `Simulation`
simulation = Simulation(model, Δt = 1e-3, stop_time = 0.5) # Should be able to go to t = 10, look into this

function progress(sim)
    model = sim.model
    @info @sprintf("Simulation time: %s",prettytime(sim.model.clock.time))
     return nothing
 end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))

# And finally run the simulation.
run!(simulation)

# So the combination is GPU, prescribed velocities, and field dependent forcing
#We can do any two of these but not all three. For a 2D grid and HydrostaticFreeSurface model, we do need prescribed velocities otherwise the grid isn't valid. 

