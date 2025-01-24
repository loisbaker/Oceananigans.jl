using Oceananigans
using Oceananigans.Operators
using Oceananigans.OutputReaders
using Oceananigans.OutputReaders: OnDisk
using Oceananigans.OutputReaders: GPUAdaptedFieldTimeSeries
using Oceananigans.Units
using Oceananigans.Utils: Time 
using Printf
using Oceananigans.Architectures: on_architecture
using CairoMakie

include("generate_forcing_data.jl")

#####
##### Periodic channel with time-dependent surface boundary fluxes 
#####

stop_time = 365days
forcing_frequency = 1day

arch = GPU()

# Defining the grid 
grid = LatitudeLongitudeGrid(arch;
                             size = (60, 60, 10), 
                         latitude = (15, 75), 
                        longitude = (0, 60), 
                         topology = (Periodic, Bounded, Bounded),
                             halo = (4, 4, 4),
                                z = (-1000, 0))

#####
##### Generate Forcing Data
#####

# Create a time series of atmospheric data with timestep 1 day
# running for 1 year. (This step is not necessary if the input file is already available)
times = range(0, stop_time, step = forcing_frequency)
forcing_file = "forcing_data.jld2"

# T_forcing_tmp = Field((Center, Center, Center), grid) 
# T_forcing = FieldTimeSeries{Center, Center, Center}(grid, times; backend = OnDisk(), path = forcing_file, name = "T_forcing")

# for (t, time) in enumerate(T_forcing.times)
#     @info "writing down data for timestep $t and time $time"
#     set!(T_forcing_tmp, (x, y, z) -> sin(π * x) * time )
#     set!(T_forcing,T_forcing_tmp,t)
# end

generate_forcing_data!(grid, times, forcing_file)


T_forcing = FieldTimeSeries(forcing_file, "T_forcing"; backend = InMemory(2))
T_top_bc = FluxBoundaryCondition(0) 
u_top_bc = FluxBoundaryCondition(0)


T_bcs = FieldBoundaryConditions(top = T_top_bc)
u_bcs = FieldBoundaryConditions(top = u_top_bc)

# #####
# ##### Physical and Numerical Setup
# #####

momentum_advection = VectorInvariant(vorticity_scheme = WENO(), 
                                      vertical_scheme = WENO())

tracer_advection = WENO()

buoyancy = SeawaterBuoyancy(equation_of_state = LinearEquationOfState(), constant_salinity = 35)

free_surface = SplitExplicitFreeSurface(grid, cfl = 0.7)

coriolis = HydrostaticSphericalCoriolis()

convective_adjustment  = ConvectiveAdjustmentVerticalDiffusivity(convective_κz = 0.5)
vertical_diffusivity   = VerticalScalarDiffusivity(κ = 1e-5, ν = 1e-4)
horizontal_diffusivity = HorizontalScalarBiharmonicDiffusivity(κ = 1e9, ν = 1e9)

closure = (convective_adjustment, vertical_diffusivity, horizontal_diffusivity)

#####
##### Create the model and initial conditions
#####

model = HydrostaticFreeSurfaceModel(; 
                                     grid, 
                                     momentum_advection,
                                     tracer_advection,
                                     tracers = :T,
                                     forcing = (;T=T_forcing),
                                     free_surface, buoyancy, coriolis, closure,
                                     boundary_conditions = (T = T_bcs, u = u_bcs))

Tᵢ(x, y, z) = 2 * (1 + z / 1000)
                                    
set!(model, T = Tᵢ)

#####
##### Simulation setup
#####

simulation = Simulation(model; Δt = 15minutes, stop_time)

function progress(sim)
    model = sim.model
    u, v, w = model.velocities
    T = model.tracers.T 
    @info @sprintf("Simulation time: %s, max(|u|, |v|, |w|, |T|): %.2e, %.2e, %.2e, %.2e \n", 
                   prettytime(sim.model.clock.time), 
                   maximum(abs, u), maximum(abs, v), 
                   maximum(abs, w), maximum(abs, T))

    return nothing
end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

output_file = "seasonal_baroclinic_channel.jld2"

simulation.output_writers[:fields] = JLD2OutputWriter(model, merge(model.velocities, model.tracers);
                                                      schedule = TimeInterval(1day), 
                                                      filename = output_file,
                                                      overwrite_existing = true)

run!(simulation)

# #####
# ##### Visualize the simulation!!
# #####

# # T_series = FieldTimeSeries(output_file, "T")
# # u_series = FieldTimeSeries(output_file, "u")
# # v_series = FieldTimeSeries(output_file, "v")

# # iter = Observable(1)

# # Tt = @lift(interior(T_series[$iter], :, :, 10))
# # ut = @lift(interior(u_series[$iter], :, :, 10))
# # vt = @lift(interior(v_series[$iter], :, :, 10))

# # fig = Figure()
# # ax = Axis(fig[1, 1], title = "top T")
# # heatmap!(ax, Tt, colormap = :thermal, colorrange = (0, 1))
# # ax = Axis(fig[2, 1], title = "top u")
# # heatmap!(ax, ut, colormap = :viridis, colorrange = (-0.1, 0.1))
# # ax = Axis(fig[2, 2], title = "top v")
# # heatmap!(ax, vt, colormap = :viridis, colorrange = (-0.1, 0.1))

# # CairoMakie.record(fig, "results.mp4", 1:length(T_series), framerate = 10) do i
# # #CairoMakie.record(fig, "results.mp4", 1:365, framerate = 10) do i
# #     @info "frame $i of $(length(T_series))"
# #     iter[] = i
# # end
