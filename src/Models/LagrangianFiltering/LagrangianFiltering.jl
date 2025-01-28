module LagrangianFiltering

export LagrangianFilter

using DocStringExtensions

using KernelAbstractions: @index, @kernel

using Oceananigans.Utils
using Oceananigans.Grids
using Oceananigans.Solvers

using Oceananigans.DistributedComputations
using Oceananigans.DistributedComputations: reconstruct_global_grid, Distributed
using Oceananigans.Grids: XYRegularRG, XZRegularRG, YZRegularRG, XYZRegularRG
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid
using Oceananigans.Utils: sum_of_velocities

import Oceananigans: fields, prognostic_fields 
import Oceananigans.Advection: cell_advection_timescale
import Oceananigans.TimeSteppers: step_lagrangian_particles!

export set_data_on_disk!, create_tracers, set_forcing_params
export c_div_U

#####
##### LagrangianFilter definition
#####
include("filter_utils.jl")
include("lagrangian_filtering_advection_operators.jl")
include("lagrangian_filter.jl")
include("show_lagrangian_filter.jl")
include("set_lagrangian_filter.jl")

#####
##### AbstractModel interface
#####

function cell_advection_timescale(model::LagrangianFilter)
    grid = model.grid
    velocities = model.velocities
    return cell_advection_timescale(grid, velocities)
end

"""
    fields(model::LagrangianFilter)

Return a flattened `NamedTuple` of the fields in `model.velocities`, `model.tracers`, and any
auxiliary fields for a `LagrangianFilter` model.
"""
fields(model::LagrangianFilter) = merge(model.velocities,
                                           model.tracers,
                                           model.auxiliary_fields)

"""
    prognostic_fields(model::LagrangianFilter)

Return a flattened `NamedTuple` of the prognostic fields associated with `LagrangianFilter`.
"""
prognostic_fields(model::LagrangianFilter) = merge(model.velocities, model.tracers)

# Unpack model.particles to update particle properties. See Models/LagrangianParticleTracking/LagrangianParticleTracking.jl
step_lagrangian_particles!(model::LagrangianFilter, Δt) = step_lagrangian_particles!(model.particles, model, Δt)

include("update_lagrangian_filter_state.jl")
include("lagrangian_filter_tendency_kernel_functions.jl")
include("compute_lagrangian_filter_tendencies.jl")
include("compute_lagrangian_filter_buffer_tendencies.jl")

end # module

