
using Oceananigans
using Oceananigans.Models: ShallowWaterModel, HydrostaticFreeSurfaceModel, NonhydrostaticModel
using MAT
using Printf
using Oceananigans.Units: Time
using JLD2
using JLD2: Group

function set_times_on_disk!(filename; direction = "backward")
    jldopen(fields_filename,"r+") do file

        if !haskey(file, "direction")
            # It should currently be forward, because direction hasn't been assigned yet
            file["direction"] = "forward"
        end    
        current_direction = file["direction"]
        println("Current direction is $current_direction")
        if current_direction != direction
            Base.delete!(file, "direction")
            file["direction"] = direction
            iterations = parse.(Int, keys(file["timeseries/t"]))
            times = [file["timeseries/t/$iter"] for iter in iterations]
            T = maximum(times)
            Base.delete!(file, "timeseries/t")
            g = Group(file, "timeseries/t")
            for (i,iter) in enumerate(reverse(iterations))
                g["$iter"] = times[i]
            end

            # We should also reorder all elements so they are indexed correctly.
            #Assuming the actual iteration number doesn't matter, just the order in the jld2 group
            for var in keys(file["timeseries"])[1:end-1] # not t
                for iter in reverse(keys(file["timeseries/$var"])[2:end]) # not serialized
                    data = file["timeseries/$var/$iter"]
                    Base.delete!(file, "timeseries/$var/$iter")
                    file["timeseries/$var/$iter"] = data
                end
            end

            println("New direction is $direction")
        else
            iterations = parse.(Int, keys(file["timeseries/t"]))
            times = [file["timeseries/t/$iter"] for iter in iterations]
            T = maximum(times)
            println("No need to update direction")
        end
        return T    
    end
    
end

# Something is happening in this bit to make the file weird in terms of maximum velocity

fields_filename = joinpath(@__DIR__, "SW_vort.jld2")
const T = set_times_on_disk!(fields_filename, direction="backward")

arch = GPU()

# Load in the velocities. Need to work out the best backend. 
u_t = FieldTimeSeries(fields_filename, "u"; architecture=arch, backend=InMemory())
v_t = FieldTimeSeries(fields_filename, "v"; architecture=arch, backend=InMemory())
ω_t = FieldTimeSeries(fields_filename, "ω"; architecture=arch, backend=InMemory())
grid = u_t.grid
print("Sanity check u  = $(maximum(u_t[Time(0)]))")

# Define the forcing for the tracers

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
gC_forcing_func(x,y,t,gC,gS,p) = -p.c*gC - p.d*gS
gS_forcing_func(x,y,t,gC,gS,p) = -p.c*gS + p.d*gC

gC_forcing1 = ω_t
gC_forcing2 = Forcing(gC_forcing_func, parameters =forcing_params, field_dependencies = (:gC,:gS))
gS_forcing = Forcing(gS_forcing_func, parameters =forcing_params, field_dependencies = (:gC,:gS))

forcing=(; gC=(gC_forcing1,gC_forcing2), gS=gS_forcing)

# Define model (still solving for momentum, need to disable this)
model = NonhydrostaticModel(;grid, tracers= tracers, forcing=forcing)

set!(model, u=-u_t[Time(0)], v= -v_t[Time(0)])

u = model.velocities.u
v = model.velocities.v
ω = Field(∂x(v) - ∂y(u))
gC = model.tracers.gC
gS = model.tracers.gS
g_total = forcing_params.a*gC + forcing_params.b*gS

# ## Running a `Simulation`
simulation = Simulation(model, Δt = 1e-3, stop_time = T) 

function update_velocities!(sim)
    model = sim.model
    time = sim.model.clock.time
    set!(model, u=-u_t[Time(time)], v= -v_t[Time(time)])
    return nothing
end

simulation.callbacks[:update_velocities] = Callback(update_velocities!)


function progress(sim)
    model = sim.model
    u = model.velocities.u
    v = model.velocities.v
    gC = model.tracers.gC
    ω = Field(∂x(v) - ∂y(u))
    @info @sprintf("Simulation time: %s, max(|u|, |v|), max(|gC|), min(|gC|): %.2e, %.2e, %.2e, %.2e \n", 
                   prettytime(sim.model.clock.time), 
                   maximum(abs, u), maximum(abs, v),maximum(abs, gC),minimum(abs, gC))
   
     return nothing
 end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))


output_filename = joinpath(@__DIR__, "backward_LF.nc")
simulation.output_writers[:fields] = NetCDFOutputWriter(model, (; ω,g_total),
                                                        filename = output_filename,
                                                        schedule = TimeInterval(0.1),
                                                        overwrite_existing = true)
# simulation.output_writers[:fields] = NetCDFOutputWriter(model, (; u,v,ω),
#                                                         filename = output_filename,
#                                                         schedule = IterationInterval(4),
#                                                         overwrite_existing = true)

# # # And finally run the simulation.
# start = time()
run!(simulation)
# computation_time = time() - start
# println("Computation time: ", computation_time, " seconds")

# # Simulation time: 0 seconds, max(|u|, |v|), max(|gC|), min(|gC|): 3.18e+01, 2.02e+01,