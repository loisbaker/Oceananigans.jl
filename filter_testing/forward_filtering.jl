
using Oceananigans
using Oceananigans.Models: LagrangianFilter
using Oceananigans.Models.LagrangianFiltering: set_times_on_disk!, create_tracers, set_forcing_params, create_forcing
using Printf
using Oceananigans.Units: Time


fields_filename = joinpath(@__DIR__, "SW_vort.jld2")
arch = GPU()
const T = set_times_on_disk!(fields_filename, direction="forward")
N = 2

# Load in the velocities. Need to work out the best backend. 
u_t = FieldTimeSeries(fields_filename, "u"; architecture=arch, backend=InMemory(2))
v_t = FieldTimeSeries(fields_filename, "v"; architecture=arch, backend=InMemory(2))
ω_t = FieldTimeSeries(fields_filename, "ω"; architecture=arch, backend=InMemory(2))
grid = u_t.grid

forcing_params = set_forcing_params(N=N,freq_c=2)

tracers= create_tracers(forcing_params)

forcing = create_forcing(ω_t, forcing_params)

# Define model 
model = LagrangianFilter(;grid, tracers= tracers, forcing=forcing)

set!(model, u=u_t[1], v= v_t[1])

u = model.velocities.u
v = model.velocities.v
ω = Field(∂x(v) - ∂y(u))
gC1 = model.tracers.gC1
# gS = model.tracers.gS
# g_total = forcing_params.a*gC + forcing_params.b*gS

#g_total = CenterField(grid)
# for i in 1:N
#     g_total .+ (forcing_params.a[i]*model.tracers[Symbol("gC$i")] + forcing_params.b[i]*model.tracers[Symbol("gS$i")])
# end
#g_total = forcing_params.a[1]*model.tracers[:gC1] + forcing_params.b[1]*model.tracers[:gS1]
## Running a `Simulation`
simulation = Simulation(model, Δt = 1e-3, stop_time = T) 

function update_velocities!(sim)
    model = sim.model
    time = sim.model.clock.time
    set!(model, u=u_t[Time(time)], v= v_t[Time(time)])
    return nothing
end

simulation.callbacks[:update_velocities] = Callback(update_velocities!)


function progress(sim)
    model = sim.model
    u = model.velocities.u
    v = model.velocities.v
    g_total = model.tracers.gC1
    @info @sprintf("Simulation time: %s, max(|u|, |v|), max(|gC|), min(|gC|): %.2e, %.2e, %.2e, %.2e \n", 
                   prettytime(sim.model.clock.time), 
                   maximum(abs, u), maximum(abs, v),maximum(abs, g_total),minimum(abs, g_total))
   
     return nothing
 end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))


output_filename = joinpath(@__DIR__, "forward_LF.nc")
simulation.output_writers[:fields] = NetCDFOutputWriter(model, (; gC1,ω),
                                                        filename = output_filename,
                                                        schedule = TimeInterval(0.1),
                                                        overwrite_existing = true)


# # And finally run the simulation.
start = time()
run!(simulation)
computation_time = time() - start
println("Computation time: ", computation_time, " seconds")

