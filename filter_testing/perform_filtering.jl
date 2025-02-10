using Oceananigans
using Oceananigans.Models: LagrangianFilter
using Oceananigans.Models.LagrangianFiltering
using Oceananigans.Units: Time
using Printf

fields_filename = joinpath(@__DIR__, "SW_vort.jld2")
T_start = 0
T_end = 10

arch = GPU()

# Set filter order
N = 2 

# Manipulate data on disk to have correct order 
T = set_data_on_disk!(fields_filename, direction="forward", T_start = T_start, T_end = T_end)

# Define tracers to filter
filter_tracer_names = ("ω","u")

# Define velocities to use for filtering
velocity_names = ("u","v")

# Load in saved data from simulation
saved_velocities, saved_tracers, grid = load_data(fields_filename, filter_tracer_names, velocity_names = velocity_names, architecture=arch, backend=InMemory(4))

# Set filtering parameters
filter_params = set_filter_params(N=N,freq_c=2)

# Create all the tracers we'll need to solve for
tracers = create_tracers(filter_tracer_names, velocity_names, filter_params)

# Create forcing for these tracers
forcing = create_forcing(tracers, saved_tracers, filter_tracer_names, velocity_names, filter_params)

# Define model 
model = LagrangianFilter(;grid, tracers= tracers, forcing=forcing)

# Define some outputs
u = model.velocities.u
v = model.velocities.v
ω = Field(@at (Center,Center,Center) ∂x(v) - ∂y(u))
filtered_outputs = create_output_fields(model, filter_tracer_names, velocity_names, filter_params)

# Define the filtering simulation 
simulation = Simulation(model, Δt = 1e-3, stop_time = T) 

simulation.callbacks[:update_velocities] = Callback(update_velocities!, parameters = saved_velocities)


function progress(sim)
    @info @sprintf("Simulation time: %s, max(|u|, |v|), max(|g_total|), min(|g_total|): %.2e, %.2e, %.2e, %.2e \n", 
                   prettytime(sim.model.clock.time), 
                   maximum(abs, u), maximum(abs, v),minimum(abs, u),minimum(abs, v))             
     return nothing
 end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(50))


output_filename = joinpath(@__DIR__, "forward_LF.nc")
simulation.output_writers[:fields] = NetCDFOutputWriter(model, filtered_outputs,
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

simulation.output_writers[:fields] = NetCDFOutputWriter(model, filtered_outputs,
                                                        filename = output_filename,
                                                        schedule = TimeInterval(0.1),
                                                        overwrite_existing = true)

# And finally run the simulation.
start = time()
run!(simulation)
computation_time = time() - start
println("Computation time: ", computation_time, " seconds")                                                        

