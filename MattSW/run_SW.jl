"""
Simulations of the shallow water equations for the two-basin problem.

The basin shapes are set using the contours of the function

f(x, y) = |z² - 1|    where    z = x + iy,

This ensures that the problem can be transformed simply under the
conformal mapping transformation

σ + iτ = log(z² - 1).

The resulting (linearised) PV equation can now be solved for the
topographic Rossby wave modes in (σ, τ) space.

This code include the initial value problem with a Gaussian height
profile initial condition, and the forced problem with a Gaussian
forcing in the x velocity. These lines can be commented out as
required. Note that we solve the equations in conservative form so
our fields are h, uh and vh, rather than h, u and v.

Matthew N. Crowe

Project with Michael Nguyen and Ted Johnson.

"""

using Oceananigans
using Oceananigans.Models: ShallowWaterModel
using Oceananigans.Models.ShallowWaterModels: ShallowWaterScalarDiffusivity
using NCDatasets

# Set some physical and numerical parameters:

f, g = 10, 1                          # Coriolis parameter and gravity
Nx, Ny = 1024, 1024                   # number of gridpoints
σ₀, σ₁, σ₂, b = 0.6, -0.5, -0.5, 1    # definitions for bottom bathymetry
H₀, L = 0.1, 0.6                      # height and width of initial perburbation
ν = 1e-5                              # viscosity (not always needed for WENO)
δ = 0.01                              # smoothing lengthscale
parameters = ()                       # additional parameters for forcing
T = 100 / f                           # simulation stop time
saves = 400                           # number of saves

# Define the bottom topography:

σ(x, y) = log((x^2 - y^2 - 1)^2 + 4x^2*y^2) / 2

H(x, y) = -exp(-b*(σ₀ - σ₂)) + (σ(x, y) > σ₀ ? exp(-b*(σ₀ - σ₂)) : (x > 0 ? (σ(x, y) > σ₂ ?
    exp(-b*(σ(x, y) - σ₂)) : 1) : (σ(x, y) > σ₁ ? exp(-b*(σ(x, y) - σ₂)) : exp(-b*(σ₁ - σ₂)))))

Smooth(F, x, y) = (F(x + δ, y + δ) + F(x + δ, y) + F(x + δ, y - δ) +
                   F(x, y + δ)     + F(x, y)     + F(x, y - δ)     +
                   F(x - δ, y + δ) + F(x - δ, y) + F(x - δ, y - δ)) / 9

Hs(x, y) = Smooth(H, x, y)

# Set the initial height and volume fluxes:

h₀(x, y)  = exp(-b*(σ₀ - σ₂)) + H₀ * exp(-(x^2 + y^2)/L^2)
uh₀(x, y) = 0                 + g / f * H₀ * 2y/L^2 * exp(-(x^2 + y^2)/L^2) * h₀(x, y)
vh₀(x, y) = 0                 - g / f * H₀ * 2x/L^2 * exp(-(x^2 + y^2)/L^2) * h₀(x, y)

# Create grid:

grid = RectilinearGrid(GPU(),
                       size = (Nx, Ny),
                       x = (-2, 2),
                       y = (-2, 2),
                       topology = (Bounded, Bounded, Flat)) # Periodic

# Define forcing functions in uh and vh evolution equations:

uh_forcing(x, y, t, h, p) = 0
vh_forcing(x, y, t, h, p) = 0 # 0.1 * exp(-(x^2 + (y-sqrt(2))^2) / 0.2^2) * h

uh_Forcing = Forcing(uh_forcing, parameters=parameters, field_dependencies = :h)
vh_Forcing = Forcing(vh_forcing, parameters=parameters, field_dependencies = :h)

# Create model:

model = ShallowWaterModel(;
                          grid,
                          coriolis = FPlane(f=f),
                          gravitational_acceleration = g,
                          timestepper = :RungeKutta3,
                          momentum_advection = WENO(order=5), # slow, doesn't need viscosity
                          #momentum_advection = Centered(),   # fast, but less stable
                          closure = ShallowWaterScalarDiffusivity(ν=ν),
                          bathymetry = Hs,
                          forcing=(uh=uh_Forcing, vh=vh_Forcing),
                          )

# Set initial condition:

set!(model, h = h₀, uh = uh₀, vh = vh₀)

# Define vertical vorticity:

uh, vh, h = model.solution
B = model.bathymetry

u, v = uh / h, vh / h
ω = Field(∂x(v) - ∂y(u))

# Set timestep using surface wave speed:

Δt = 0.5 * (grid.Lx / grid.Nx) / sqrt(g * h₀(grid.Lx/2, grid.Ly/2))

# Define simulation, callbacks and saves:

simulation = Simulation(model, Δt = Δt, stop_time = T)

callback(sim) = @info "Iteration: $(sim.model.clock.iteration), " *
    "$(round(sim.model.clock.time / T * 100, sigdigits = 3))% complete, " *
    "t = $(round(sim.model.clock.time, sigdigits = 3))"

simulation.callbacks[:disp_time] = Callback(callback, IterationInterval(100))

fields_filename_nc = joinpath(@__DIR__, "shallow_water.nc")
fields_filename_jld2 = joinpath(@__DIR__, "shallow_water.jld2")
simulation.output_writers[:fields] = NetCDFWriter(model,
                                                        (; ω, u, v, h, B),
                                                        filename = fields_filename_nc,
                                                        schedule = TimeInterval(T/saves),
                                                        overwrite_existing = true)


simulation.output_writers[:fields_jld2] = JLD2Writer(model, (; ω,u,v,h),
                                                        filename = fields_filename_jld2,
                                                        schedule = TimeInterval(T/saves),
                                                        overwrite_existing = true)                                                        
# Run simulation:

run!(simulation)
