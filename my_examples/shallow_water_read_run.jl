

# Remember we've hacked the source code in solution_and_tracer_tendencies.jl to give us modified shallow water

using Oceananigans
using Oceananigans.Models.HydrostaticFreeSurfaceModels: HydrostaticFreeSurfaceModel, PrescribedVelocityFields
using MAT
using Printf

fields_filename = "SW_vort.jld2"
arch = GPU()
print(fields_filename)
v_t = FieldTimeSeries(fields_filename, "v"; architecture=arch, backend=InMemory(2))
u_t = FieldTimeSeries(fields_filename, "u"; architecture=arch, backend=InMemory(2))
grid = u_t.grid

model = NonhydrostaticModel(; grid, 
                          timestepper = :RungeKutta3,
                          tracers = nothing,
                          closure = nothing,
                          buoyancy= nothing)
# u_i = interior(u_t[1])
# set!(model, (u = u_i)) works


#set!(model, u = u_t[1])



# simulation = Simulation(model, Δt = 1e-3, stop_time = 10)

# function progress(sim)
#     model = sim.model
#     uh, vh, h = model.solution
#     @info @sprintf("Simulation time: %s, max(|uh|, |vh|, |h|): %.2e, %.2e, %.2e \n", 
#                    prettytime(sim.model.clock.time), 
#                    maximum(abs, uh), maximum(abs, vh), 
#                    maximum(abs, h))

#     return nothing
# end
# simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))

# # Build the `output_writer` for the two-dimensional fields to be output.
# # Output every `t = 1.0`.

# fields_filename = joinpath(@__DIR__, "SW_vort_forced.nc")
# simulation.output_writers[:fields] = NetCDFOutputWriter(model, (; ω,u),
#                                                         filename = fields_filename,
#                                                         schedule = TimeInterval(1),
#                                                         overwrite_existing = true)


# # And finally run the simulation.

# run!(simulation)





