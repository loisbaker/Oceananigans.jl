
using Oceananigans
using Oceananigans.Models: LagrangianFilter
using Oceananigans.Models.LagrangianFiltering: set_times_on_disk!, create_tracers, set_forcing_params, create_forcing, sum_tracers, update_velocities!
using Printf
using Oceananigans.Units: Time


fields_filename = joinpath(@__DIR__, "SW_vort.jld2")
arch = GPU()
N = 2

const T = set_times_on_disk!(fields_filename, direction="backward")

# Load in the velocities. Need to work out the best backend. 
u_t = FieldTimeSeries(fields_filename, "u"; architecture=arch, backend=InMemory()) # Works in memory(), in memory(2) errors at the end
v_t = FieldTimeSeries(fields_filename, "v"; architecture=arch, backend=InMemory())
ω_t = FieldTimeSeries(fields_filename, "ω"; architecture=arch, backend=InMemory())

# Some weird bug in FieldTimeSeries. Errors the first time we try to access the final element, but not after
# function try_fts_error(fts,index)

#     try
#         fts[index]
#     catch 
#         println("errored")
#     else
#         println("no error")
#     end
#     return nothing
# end

# try_fts_error(u_t,201)
# try_fts_error(v_t,201)
# try_fts_error(ω_t,201)

# try_fts_error(u_t,201)
# try_fts_error(v_t,201)
# try_fts_error(ω_t,201)


grid = u_t.grid

forcing_params = set_forcing_params(N=N,freq_c=2)

tracers= create_tracers(forcing_params)

forcing = create_forcing(ω_t, forcing_params)

# Define model 
model = LagrangianFilter(;grid, tracers= tracers, forcing=forcing)

set!(model, u=-u_t[1], v= -v_t[1])

u = model.velocities.u
v = model.velocities.v
ω = Field(∂x(v) - ∂y(u))
g_total = sum_tracers(model,forcing_params)

## Running a `Simulation`
simulation = Simulation(model, Δt = 1e-3, stop_time = T) 

simulation.callbacks[:update_velocities] = Callback(update_velocities!, parameters = (u_t, v_t))


function progress(sim)
    @info @sprintf("Simulation time: %s, max(|u|, |v|), max(|g_total|), min(|g_total|): %.2e, %.2e, %.2e, %.2e \n", 
                   prettytime(sim.model.clock.time), 
                   maximum(abs, u), maximum(abs, v),maximum(abs, g_total),minimum(abs, g_total))
   
     return nothing
 end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))


output_filename = joinpath(@__DIR__, "backward_LF.nc")
simulation.output_writers[:fields] = NetCDFOutputWriter(model, (; g_total,ω),
                                                        filename = output_filename,
                                                        schedule = TimeInterval(0.1),
                                                        overwrite_existing = true)


# And finally run the simulation.
start = time()
run!(simulation)
computation_time = time() - start
println("Computation time: ", computation_time, " seconds")

# Work out why the backwards is working in the separate script and not here


# # Now, run it backwards
# set_times_on_disk!(fields_filename, direction="backward");

# # Load in the velocities. Need to work out the best backend. Possibly no need to reload here?

# u_t = FieldTimeSeries(fields_filename, "u"; architecture=arch, backend=InMemory())
# v_t = FieldTimeSeries(fields_filename, "v"; architecture=arch, backend=InMemory())
# ω_t = FieldTimeSeries(fields_filename, "ω"; architecture=arch, backend=InMemory())

# # Recreate forcing with backward scalar timeseries. No need to reset tracers. Might not be needed either
# forcing = create_forcing(ω_t, forcing_params)

# # Define model 
# model = LagrangianFilter(;grid, tracers= tracers, forcing=forcing)

# # Set initial velocities
# set!(model, u=-u_t[Time(0)], v= -v_t[Time(0)])

# u = model.velocities.u
# v = model.velocities.v
# ω = Field(∂x(v) - ∂y(u))
# g_total = sum_tracers(model,forcing_params)

# ## Running a `Simulation`
# simulation = Simulation(model, Δt = 1e-3, stop_time = T) 

# simulation.callbacks[:update_velocities] = Callback(update_velocities!, parameters = (-u_t, -v_t))

# simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))

# output_filename = joinpath(@__DIR__, "backward_LF.nc")

# simulation.output_writers[:fields] = NetCDFOutputWriter(model, (; g_total,ω),
#                                                         filename = output_filename,
#                                                         schedule = TimeInterval(0.1),
#                                                         overwrite_existing = true)

# # And finally run the simulation.
# start = time()
# run!(simulation)
# computation_time = time() - start
# println("Computation time: ", computation_time, " seconds")                                                        