# # An unstable Bickley jet in Shallow Water model 
#
# This example uses Oceananigans.jl's `ShallowWaterModel` to simulate
# the evolution of an unstable, geostrophically balanced, Bickley jet
# The example is periodic in ``x`` with flat bathymetry and
# uses the conservative formulation of the shallow water equations.
# The initial conditions superpose the Bickley jet with small-amplitude perturbations.
# See ["The nonlinear evolution of barotropically unstable jets," J. Phys. Oceanogr. (2003)](https://doi.org/10.1175/1520-0485(2003)033<2173:TNEOBU>2.0.CO;2)
# for more details on this problem.
#
# The mass transport ``(uh, vh)`` is the prognostic momentum variable
# in the conservative formulation of the shallow water equations,
# where ``(u, v)`` are the horizontal velocity components and ``h``
# is the layer height. 
#
# ## Install dependencies
#
# First we make sure that we have all of the packages that are required to
# run the simulation.
#
# ```julia
# using Pkg
# pkg"add Oceananigans, NCDatasets, Polynomials, CairoMakie"
# ```

# Remember we've hacked the source code in solution_and_tracer_tendencies.jl to give us modified shallow water

using Oceananigans
using Oceananigans.Models: ShallowWaterModel
using MAT
using Printf

# ## Two-dimensional domain 
#
# The shallow water model is two-dimensional and uses grids that are `Flat`
# in the vertical direction. We use length scales non-dimensionalized by the width
# of the Bickley jet.

grid = RectilinearGrid(GPU(),size = (256, 256),
                       x = (0, 2π),
                       y = (0, 2π),
                       topology = (Periodic, Periodic, Flat))

# ## Building a `ShallowWaterModel`
#
# We build a `ShallowWaterModel` with the `WENO` advection scheme,
# 3rd-order Runge-Kutta time-stepping, non-dimensional Coriolis, and
# gravitational acceleration

Fr = 0.3
Ro = 0.4 # These match geostrophic balance of the ICs
gravitational_acceleration = 1/Fr^2
coriolis = FPlane(f=1/Ro)

model = ShallowWaterModel(; grid, coriolis, gravitational_acceleration,
                          timestepper = :RungeKutta3,
                           momentum_advection = WENO())

# ## Background state and perturbation
#
# The background velocity ``ū`` and free-surface ``η̄`` correspond to a
# geostrophically balanced Bickely jet with maximum speed of ``U`` and maximum 
# free-surface deformation of ``Δη``,

f = coriolis.f
g = gravitational_acceleration

# Get initial condition from mat file
file = matopen("uvh_256_Fr_0_3_Ro_0_4_wave_0_5.mat")
u_i = read(file, "u") 
v_i = read(file, "v")
h_i = read(file, "h")

# u_i = copy(transpose(read(file, "ur"))) # Just transpose gives a datatype the GPOU doesn't like
# v_i = copy(transpose(read(file,"vr")))
# h_i = copy(transpose(read(file,"h")))
close(file)

# # We first set a "clean" initial condition without noise for the purpose of discretely
# # calculating the initial 'mean' vorticity,

uh_i = u_i.*h_i
vh_i = v_i.*h_i

set!(model, uh = uh_i, vh = vh_i, h= h_i )

## Build velocities and vorticity
uh, vh, h = model.solution
u = Field(uh / h) # This is fine when these things aren't arrays, otherwise need elementwise division
v = Field(vh / h)
ω = Field(∂x(v) - ∂y(u))

## Running a `Simulation`

simulation = Simulation(model, Δt = 1e-3, stop_time = 10)

function progress(sim)
    model = sim.model
    uh, vh, h = model.solution
    @info @sprintf("Simulation time: %s, max(|uh|, |vh|, |h|): %.2e, %.2e, %.2e \n", 
                   prettytime(sim.model.clock.time), 
                   maximum(abs, uh), maximum(abs, vh), 
                   maximum(abs, h))

    return nothing
end
simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))

# Build the `output_writer` for the two-dimensional fields to be output.
# Output every `t = 1.0`.

fields_filename = joinpath(@__DIR__, "SW_vort.jld2")
simulation.output_writers[:fields] = JLD2OutputWriter(model, (; ω,u),
                                                        filename = fields_filename,
                                                        schedule = TimeInterval(1),
                                                        overwrite_existing = true)


# And finally run the simulation.

run!(simulation)

# Now try to load back in the data


u_t = FieldTimeSeries(fields_filename, "u")
w_t = FieldTimeSeries(fields_filename, "ω")
times = u_t.times




# # Visualize the results
# # Load required packages to read output and plot.

# using NCDatasets, Printf, CairoMakie
# nothing #hide

# # # Define the coordinates for plotting.
# xω, yω = xnodes(ω), ynodes(ω)
# xu, yu = xnodes(u), ynodes(u)
# nothing #hide

# # Read in the `output_writer` for the two-dimensional fields and then create an animation 
# # showing both the perturbation vorticities.

# fig = Figure(size = (660, 750))

# axis_kwargs = (xlabel = "x", ylabel = "y")
# ax_ω  = Axis(fig[2, 1]; title = "Total vorticity, ω", axis_kwargs...)
# # ax_ω′ = Axis(fig[2, 3]; title = "Horizontal velocity, u", axis_kwargs...)

# n = Observable(1)

# ds = NCDataset(simulation.output_writers[:fields].filepath, "r")

# times = ds["time"][:]

# ω = @lift ds["ω"][:, :, 1, $n]
# # hm_ω = heatmap!(ax_ω, x, y, ω, colorrange = (-1, 1), colormap = :balance)
# hm_ω = heatmap!(ax_ω, xω, yω, ω, colorrange = (-1, 1), colormap = :balance)
# Colorbar(fig[2, 2], hm_ω)

# # u = @lift ds["u"][:, :, 1, $n]
# # # hm_ω′ = heatmap!(ax_ω′, x, y, ω′, colorrange = (-1, 1), colormap = :balance)
# # hm_ω′ = heatmap!(ax_ω′, xu, yu, u, colormap = :balance)
# # Colorbar(fig[2, 4], hm_ω′)

# title = @lift @sprintf("t = %.1f", times[$n])
# fig[1, 1:2] = Label(fig, title, fontsize=24, tellwidth=false)

# current_figure() #hide
# fig

# # Finally, we record a movie.

# frames = 1:length(times)

# record(fig, "SW_vort.mp4", frames, framerate=12) do i
#     n[] = i
# end
# nothing #hide

# # ![](shallow_water_Bickley_jet.mp4)

# # It's always good practice to close the NetCDF files when we are done.

# close(ds)
