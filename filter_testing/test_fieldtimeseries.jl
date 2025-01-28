using Oceananigans
using Oceananigans.Models: LagrangianFilter
using Oceananigans.Models.LagrangianFiltering: set_data_on_disk!, create_tracers, set_forcing_params, create_forcing, sum_tracers, update_velocities!
using Printf
using Oceananigans.Units: Time
arch = GPU()

