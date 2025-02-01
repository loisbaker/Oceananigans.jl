
using Oceananigans
using Oceananigans.Models: LagrangianFilter
using Oceananigans.Models.LagrangianFiltering: set_data_on_disk!, create_tracers, set_forcing_params, create_forcing, sum_tracers#, update_velocities!
using Printf
using Oceananigans.Units: Time

fields_filename = joinpath(@__DIR__, "SW_vort.jld2")
# T_start = 1
# T_end = 19

#TODO 
# Output sensible filtered fields
# Get FieldTimeSeries working correctly
# Add xi maps
# Add multiple tracer functionality

# Forcing with the in memory fts seems fine, but the update velocities throws the bounds error. Maybe its not remembering? 
# Need to also check if they are loading in correctly - potentially not. Check the interpolation Time(). Something very weird going on with reversing access
arch = GPU()

N = 1 # Filter order
T = 10
#T = set_data_on_disk!(fields_filename, direction="forward", T_start = T_start, T_end = T_end)

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
gC1 = model.tracers.gC1
gS1 = model.tracers.gS1

# Running a `Simulation`
simulation = Simulation(model, Δt = 1e-3, stop_time = T) 

function update_velocities!(sim)
    model = sim.model
    time = sim.model.clock.time
    # u_fts, v_fts = fts_velocities
    set!(model, u=u_t[Time(time)], v= v_t[Time(time)])
    return nothing
end

simulation.callbacks[:update_velocities] = Callback(update_velocities!)#, parameters = (u_t, v_t))


function progress(sim)
    @info @sprintf("Simulation time: %s, max(|u|, |v|), max(|g_total|), min(|g_total|): %.2e, %.2e, %.2e, %.2e \n", 
                   prettytime(sim.model.clock.time), 
                   maximum(abs, u), maximum(abs, v),minimum(abs, u),minimum(abs, v))
    println(Time(sim.model.clock.time))               
    println(u_t.backend)
     return nothing
 end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(50))


# output_filename = joinpath(@__DIR__, "forward_LF.nc")
# simulation.output_writers[:fields] = NetCDFOutputWriter(model, (; u,v,ω,g_total,gC1,gS1),
#                                                         filename = output_filename,
#                                                         schedule = TimeInterval(0.1),
#                                                         overwrite_existing = true)


# And finally run the simulation.
start = time()
run!(simulation)
computation_time = time() - start
println("Computation time: ", computation_time, " seconds")


# Now, run it backwards. Switch the data direction on disk
# set_data_on_disk!(fields_filename, direction="backward", T_start=T_start, T_end = T_end);

# #Load in the velocities. Need to work out the best backend. 
# u_t = FieldTimeSeries(fields_filename, "u"; architecture=arch, backend=InMemory())
# v_t = FieldTimeSeries(fields_filename, "v"; architecture=arch, backend=InMemory())
# ω_t = FieldTimeSeries(fields_filename, "ω"; architecture=arch, backend=InMemory())

# tracers= create_tracers(forcing_params)

# forcing = create_forcing(ω_t, forcing_params)

# # Define model 
# model = LagrangianFilter(;grid, tracers= tracers, forcing=forcing)

# u = model.velocities.u
# v = model.velocities.v
# ω = Field(@at (Center,Center,Center) ∂x(v) - ∂y(u))
# g_total = sum_tracers(model,forcing_params)
# gC1 = model.tracers.gC1
# gS1 = model.tracers.gS1

# # Running a `Simulation`
# simulation = Simulation(model, Δt = 1e-3, stop_time = T) 

# simulation.callbacks[:update_velocities] = Callback(update_velocities!, parameters = (u_t, v_t))

# simulation.callbacks[:progress] = Callback(progress, IterationInterval(500))


# # Reset time
# #simulation.model.clock.time = 0

# output_filename = joinpath(@__DIR__, "backward_LF.nc")

# simulation.output_writers[:fields] = NetCDFOutputWriter(model, (; g_total,ω,gC1,gS1),
#                                                         filename = output_filename,
#                                                         schedule = TimeInterval(0.1),
#                                                         overwrite_existing = true)

# # And finally run the simulation.
# start = time()
# run!(simulation)
# computation_time = time() - start
# println("Computation time: ", computation_time, " seconds")                                                        

# Weirdness with going back to access 
# [ Info: Simulation time: 50.000 ms, max(|u|, |v|), max(|g_total|), min(|g_total|): 9.55e-01, 5.52e-01, 2.50e-05, 3.76e-07 
# Time{Float64}(0.0499999999999999)
# InMemory{Int64}(1, 4)
# [ Info: Simulation time: 100.000 ms, max(|u|, |v|), max(|g_total|), min(|g_total|): 9.92e-01, 5.75e-01, 1.12e-05, 6.73e-06 
# Time{Float64}(0.09999999999999987)
# InMemory{Int64}(1, 4)
# [ Info: Simulation time: 150.000 ms, max(|u|, |v|), max(|g_total|), min(|g_total|): 9.71e-01, 5.62e-01, 6.83e-05, 1.48e-06 
# Time{Float64}(0.1499999999999999)
# InMemory{Int64}(1, 4)
# [ Info: Simulation time: 200.000 ms, max(|u|, |v|), max(|g_total|), min(|g_total|): 1.01e+00, 5.81e-01, 3.57e-05, 5.67e-06 
# Time{Float64}(0.19999999999999996)
# InMemory{Int64}(1, 4)
# [ Info: Simulation time: 250 ms, max(|u|, |v|), max(|g_total|), min(|g_total|): 9.72e-01, 6.00e-01, 2.06e-05, 7.60e-06 
# Time{Float64}(0.25)
# InMemory{Int64}(1, 4)
# [ Info: Simulation time: 300.000 ms, max(|u|, |v|), max(|g_total|), min(|g_total|): 1.14e+00, 6.71e-01, 4.17e-06, 8.03e-05 
# Time{Float64}(0.30000000000000004)
# InMemory{Int64}(5, 4)
# [ Info: Simulation time: 350.000 ms, max(|u|, |v|), max(|g_total|), min(|g_total|): 9.82e-01, 6.34e-01, 1.30e-06, 1.05e-05 
# Time{Float64}(0.3500000000000001)
# InMemory{Int64}(4, 4)
# [ Info: Simulation time: 400.000 ms, max(|u|, |v|), max(|g_total|), min(|g_total|): 1.01e+00, 6.51e-01, 3.61e-05, 6.70e-07 
# Time{Float64}(0.40000000000000013)
# InMemory{Int64}(4, 4)
# [ Info: Simulation time: 450.000 ms, max(|u|, |v|), max(|g_total|), min(|g_total|): 1.02e+00, 6.29e-01, 2.23e-05, 3.78e-07 
# Time{Float64}(0.4500000000000002)
# InMemory{Int64}(4, 4)
# [ Info: Simulation time: 500.000 ms, max(|u|, |v|), max(|g_total|), min(|g_total|): 1.07e+00, 6.40e-01, 6.34e-05, 3.74e-05 
# Time{Float64}(0.5000000000000002)
# InMemory{Int64}(4, 4)
# [ Info: Simulation time: 550.000 ms, max(|u|, |v|), max(|g_total|), min(|g_total|): 1.07e+00, 6.22e-01, 1.99e-05, 1.11e-05 
# Time{Float64}(0.5500000000000003)
# InMemory{Int64}(4, 4)
# [ Info: Simulation time: 600.000 ms, max(|u|, |v|), max(|g_total|), min(|g_total|): 1.18e+00, 6.49e-01, 3.54e-05, 8.96e-08 
# Time{Float64}(0.6000000000000003)
# InMemory{Int64}(8, 4)

# Look into what status: time=0.0 means

# Logic isn't right
# Start at u_t[Time(0)] -- InMemory{Int64}(1, 4)
#u_t[Time(0.2)] -- InMemory{Int64}(1, 4)
#u_t[Time(0.3)] -- InMemory{Int64}(1, 4)
# Note that u_t.times[4] = 0.3, so we're on the edge of what's in memory
#u_t[Time(0.31)] -- InMemory{Int64}(5, 4), so we've moved to knowing 0.4 onwards - not good for knowing 0.31 but doesn't error, does get it wrong!!
#u_t[Time(0.32)] -- InMemory{Int64}(4, 4), like it corrects

# This demonstrates the error well. Need to correct logic of partly in memory fts, as this is also loading in the wrong values.
for i in 0:0.01:1
    v = maximum(u_t[Time(i)])
    println("")
    println("Next")
    println(v)
    println("Time = $i")
    println(u_t.backend)
end
 # So need to fix this and also how it goes for values where in memory would be out of range
 # This might well fix things 