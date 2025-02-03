
using Oceananigans
using Oceananigans.Models: LagrangianFilter
using Oceananigans.Models.LagrangianFiltering: set_data_on_disk!, create_tracers, set_forcing_params, create_forcing, sum_tracers, update_velocities!
using Printf
using Oceananigans.Units: Time

fields_filename = joinpath(@__DIR__, "SW_vort.jld2")
T_start = 0
T_end = 10

#TODO 
# Output sensible filtered fields
# Get FieldTimeSeries working correctly
# Add xi maps
# Add multiple tracer functionality

arch = GPU()

N = 1 # Filter order

T = set_data_on_disk!(fields_filename, direction="forward", T_start = T_start, T_end = T_end)

# Load in the velocities. Need to work out the best backend. Ideally InMemory(2), but this has errors. AND it really messes up the field. We don't trust the partly in memory functionality, needs fixing
u_t = FieldTimeSeries(fields_filename, "u"; architecture=arch, backend=InMemory(4))
v_t = FieldTimeSeries(fields_filename, "v"; architecture=arch, backend=InMemory(4))
ω_t = FieldTimeSeries(fields_filename, "ω"; architecture=arch, backend=InMemory(4))

grid = u_t.grid

forcing_params = set_forcing_params(N=N,freq_c=2)

tracers= create_tracers(forcing_params)

forcing = create_forcing(ω_t, forcing_params)

# Define model 
model = LagrangianFilter(;grid, tracers= tracers, forcing=forcing)

u = model.velocities.u
v = model.velocities.v
ω = Field(@at (Center,Center,Center) ∂x(v) - ∂y(u))
g_total = sum_tracers(model,forcing_params)

# Running a `Simulation`
simulation = Simulation(model, Δt = 1e-3, stop_time = T) 

simulation.callbacks[:update_velocities] = Callback(update_velocities!, parameters = (u_t, v_t))


function progress(sim)
    @info @sprintf("Simulation time: %s, max(|u|, |v|), max(|g_total|), min(|g_total|): %.2e, %.2e, %.2e, %.2e \n", 
                   prettytime(sim.model.clock.time), 
                   maximum(abs, u), maximum(abs, v),minimum(abs, u),minimum(abs, v))             
    println(u_t.backend)
     return nothing
 end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(50))


output_filename = joinpath(@__DIR__, "forward_LF.nc")
simulation.output_writers[:fields] = NetCDFOutputWriter(model, (; u,v,ω,g_total),
                                                        filename = output_filename,
                                                        schedule = TimeInterval(0.1),
                                                        overwrite_existing = true)


# And finally run the simulation.
start = time()
run!(simulation)
computation_time = time() - start
println("Computation time: ", computation_time, " seconds")


# Now, run it backwards. Switch the data direction on disk
set_data_on_disk!(fields_filename, direction="backward", T_start=T_start, T_end = T_end);

# Reset time
simulation.model.clock.time = 0

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

