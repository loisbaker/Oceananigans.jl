using Oceananigans
using Oceananigans.Units
using Oceananigans.Models: LagrangianFilter
using Oceananigans.Models.LagrangianFiltering
using Oceananigans.Units: Time
using Printf

fields_filename = joinpath(@__DIR__, "lee_wave_sim.jld2")
T_start = 0
T_end = 5days
dt = 50seconds
arch = GPU()

# Set filter order and cut-off frequency
# Amplitude of frequency response of filter will be squared Butterworth order 2^N
N = 1
freq_c = 1e-4

# Manipulate data on disk to have correct order 
T = set_data_on_disk!(fields_filename, direction="forward", T_start = T_start, T_end = T_end)

# Define tracers to filter
filter_tracer_names = ("b",)

# Define velocities to use for filtering
velocity_names = ("u","w")

# Load in saved data from simulation
saved_velocities, saved_tracers, grid = load_data(fields_filename, filter_tracer_names, velocity_names = velocity_names, architecture=arch, backend=InMemory(4))

# Set filtering parameters
filter_params = set_filter_params(N=N,freq_c=freq_c)

# Create all the tracers we'll need to solve for
tracers = create_tracers(filter_tracer_names, velocity_names, filter_params)

# Create forcing for these tracers
forcing = create_forcing(tracers, saved_tracers, filter_tracer_names, velocity_names, filter_params)

# Define model 
model = LagrangianFilter(;grid, tracers = tracers, forcing = forcing)

# Define our outputs
u = model.velocities.u
w = model.velocities.w
filtered_outputs = create_output_fields(model, filter_tracer_names, velocity_names, filter_params)

# Define the filtering simulation 
simulation = Simulation(model, Δt = dt, stop_time = T) 

simulation.callbacks[:update_velocities] = Callback(update_velocities!, parameters = saved_velocities)


function progress(sim)
    @info @sprintf("Simulation time: %s, max(|u|, |v|), min(|u|, |v|): %.2e, %.2e, %.2e, %.2e \n", 
                   prettytime(sim.model.clock.time), 
                   maximum(abs, u), maximum(abs, w), minimum(abs, u),minimum(abs, w))             
     return nothing
 end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(50))


output_filename = joinpath(@__DIR__, "forward_LF_LW.nc")
simulation.output_writers[:fields] = NetCDFOutputWriter(model, filtered_outputs,
                                                        filename = output_filename,
                                                        schedule = TimeInterval(5hours),
                                                        overwrite_existing = true)


# And run the simulation.
run!(simulation)

# Now, run it backwards. Switch the data direction on disk
set_data_on_disk!(fields_filename, direction="backward", T_start=T_start, T_end = T_end);

# Reset time
simulation.model.clock.time = 0

output_filename = joinpath(@__DIR__, "backward_LF_LW.nc")

simulation.output_writers[:fields] = NetCDFOutputWriter(model, filtered_outputs,
                                                        filename = output_filename,
                                                        schedule = TimeInterval(5hours),
                                                        overwrite_existing = true)

#And run the simulation.
run!(simulation)
