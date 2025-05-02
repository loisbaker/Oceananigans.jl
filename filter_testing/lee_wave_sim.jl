using Oceananigans
using Oceananigans.Units
using Oceananigans.ImmersedBoundaries: PartialCellBottom
using Printf
#using SpecialFunctions

# Set grid params
Nx=250
Nz=250
H=3kilometers
Fr=0.25
N = 1e-3
k = 0.005
h₀=25meters
f = 1e-4
stop_time=5days

# set other constants
Nsqr = N^2
x_range = 3 * pi / k
U_const = N * h₀ / Fr

# create grid
underlying_grid = RectilinearGrid(GPU(),
    size = (Nx, Nz),
    x = (-x_range, x_range),
    z = (-H, 0),
    halo = (4, 4),
    topology = (Periodic, Flat, Bounded)
)

hill(x) = h₀ * (1 + sin(k * x))
bottom(x) = -H + hill(x)

grid = ImmersedBoundaryGrid(underlying_grid, PartialCellBottom(bottom))

# v forcing
coriolis = FPlane(f = f)
v_forcing_eq(x, z, t, p) = p.U_const * p.coriolis.f
v_forcing = Forcing(v_forcing_eq, parameters = (; U_const, coriolis))

ν_diff = κ_diff = 1
horizontal_diffusivity = HorizontalScalarDiffusivity(ν = ν_diff, κ = κ_diff)

model = NonhydrostaticModel(;
    grid,
    coriolis = coriolis,
    buoyancy = BuoyancyTracer(),
    tracers = (:b),
    closure = horizontal_diffusivity,
    forcing = (; v=v_forcing)
)

# field constants and timestepping
uᵢ(x, z) = U_const

bᵢ(x, z) = Nsqr * z
dx = max(H / Nz, 2 * x_range / Nx)
Δt = dx / (2 * U_const)
println(Δt)
set!(model, u = uᵢ, b=bᵢ )

simulation = Simulation(model; Δt, stop_time)
wizard = TimeStepWizard(max_change=1.1, cfl=0.7)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(4))

wall_clock = Ref(time_ns())

function progress(sim)
    elapsed = 1e-9 * (time_ns() - wall_clock[])

    msg = @sprintf(
        "iteration: %d, time: %s, wall time: %s, max|w|: %6.3e, m s⁻¹\n",
        iteration(sim),
        prettytime(sim),
        prettytime(elapsed),
        maximum(abs, sim.model.velocities.w)
    )
 
    wall_clock[] = time_ns()
    @info msg
    return nothing
end

add_callback!(simulation, progress, name = :progress, IterationInterval(1))

b = model.tracers.b
u, v, w = model.velocities

fields_filename= joinpath(@__DIR__, "lee_wave_sim")
save_fields_interval = 1hours

simulation.output_writers[:fields_jld2] = JLD2OutputWriter(model, (; b,u,w),
                                                        filename = fields_filename,
                                                        schedule = TimeInterval(save_fields_interval),
                                                        overwrite_existing = true)

simulation.output_writers[:fields_nc] = NetCDFOutputWriter(model, (; b,u,v ),
                                                        filename = fields_filename,
                                                        schedule = TimeInterval(save_fields_interval),
                                                        overwrite_existing = true)


run!(simulation)