
using Oceananigans
using Oceananigans.Models: LagrangianFilter
using Oceananigans.Models.LagrangianFiltering: set_data_on_disk!, create_tracers, set_forcing_params, create_forcing, sum_tracers, update_velocities!
using Printf
using Oceananigans.Units: Time
arch = GPU()

fields_filename = joinpath(@__DIR__, "SW_vort_copy.jld2")
T_start = 0.1
T_end = 4.9


T = set_data_on_disk!(fields_filename, direction="forward", T_start=T_start, T_end = T_end)

# Load in the velocities. Need to work out the best backend. 
u_t = FieldTimeSeries(fields_filename, "u"; architecture=arch, backend=InMemory())
v_t = FieldTimeSeries(fields_filename, "v"; architecture=arch, backend=InMemory())
ω_t = FieldTimeSeries(fields_filename, "ω"; architecture=arch, backend=InMemory())

# Now, run it backwards. Switch the data direction on disk
set_data_on_disk!(fields_filename, direction="backward", T_start=T_start, T_end = T_end);

# Load in the velocities. Need to work out the best backend. 
u_t = FieldTimeSeries(fields_filename, "u"; architecture=arch, backend=InMemory())
v_t = FieldTimeSeries(fields_filename, "v"; architecture=arch, backend=InMemory())
ω_t = FieldTimeSeries(fields_filename, "ω"; architecture=arch, backend=InMemory())

