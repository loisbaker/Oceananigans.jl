
using Oceananigans
using Oceananigans.Models: LagrangianFilter
using Oceananigans.Models.LagrangianFiltering: set_times_on_disk!#, create_tracers, set_forcing_params
using Printf
using Oceananigans.Units: Time

fields_filename = joinpath(@__DIR__, "SW_vort.jld2")
arch = GPU()
const T = set_times_on_disk!(fields_filename, direction="forward")
N = 1

# Load in the velocities. Need to work out the best backend. 
u_t = FieldTimeSeries(fields_filename, "u"; architecture=arch, backend=InMemory(2))
v_t = FieldTimeSeries(fields_filename, "v"; architecture=arch, backend=InMemory(2))
ω_t = FieldTimeSeries(fields_filename, "ω"; architecture=arch, backend=InMemory(2))
grid = u_t.grid

tracers= (:gC,:gS)

# Let's do this for N=1 for now
function set_forcing_params(;freq_c=1)
    N = 1
    a = (freq_c/2^N)*sin(pi/(2^(N+1)))
    b = (freq_c/2^N)*cos(pi/(2^(N+1)))
    c = freq_c*sin(pi/(2^(N+1)))
    d = freq_c*cos(pi/(2^(N+1)))
    forcing_params = (a=a,b=b,c=c,d=d)
    return forcing_params
end

forcing_params = set_forcing_params(freq_c=2)

# Set forcing
# gC_forcing_func(x,y,t,gC,gS,p) = -p.c*gC - p.d*gS
# gS_forcing_func(x,y,t,gC,gS,p) = -p.c*gS + p.d*gC
# forcing_params = set_forcing_params(N=N,freq_c=2)
# tracers= create_tracers(N=N)

# Make a function factory to create forcings 
# Function factory example
function make_gC_forcing_func(i)
    # Return a new function 
    return (x,y,t,gC,gS,p) -> -p.c[i]*gC - p.d[i]*gS
    
end


# Set forcing - this works
gC_forcing_func(x,y,t,a,b,p) = -p.c*a - p.d*b
gS_forcing_func(x,y,t,a,b,p) =  -p.c*b + p.d*a

gC_forcing1 = ω_t
gC_forcing2 = Forcing(gC_forcing_func, parameters =forcing_params, field_dependencies = (:gC,:gS))
gS_forcing = Forcing(gS_forcing_func, parameters =forcing_params, field_dependencies = (:gC,:gS))

forcing=(; gC=(gC_forcing1,gC_forcing2), gS=gS_forcing)
# closure = ScalarDiffusivity(κ=0.1)


# Define model 
model = LagrangianFilter(;grid, tracers= tracers, forcing=forcing)

set!(model, u=u_t[1], v= v_t[1])

u = model.velocities.u
v = model.velocities.v
ω = Field(∂x(v) - ∂y(u))
gC = model.tracers.gC
gS = model.tracers.gS
g_total = forcing_params.a*gC + forcing_params.b*gS

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
    gC= model.tracers.gC
    ω = Field(∂x(v) - ∂y(u))
    @info @sprintf("Simulation time: %s, max(|u|, |v|), max(|gC|), min(|gC|): %.2e, %.2e, %.2e, %.2e \n", 
                   prettytime(sim.model.clock.time), 
                   maximum(abs, u), maximum(abs, v),maximum(abs, gC),minimum(abs, gC))
   
     return nothing
 end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))


output_filename = joinpath(@__DIR__, "forward_LF.nc")
simulation.output_writers[:fields] = NetCDFOutputWriter(model, (; g_total,ω),
                                                        filename = output_filename,
                                                        schedule = TimeInterval(0.1),
                                                        overwrite_existing = true)


# # And finally run the simulation.
start = time()
run!(simulation)
computation_time = time() - start
println("Computation time: ", computation_time, " seconds")

