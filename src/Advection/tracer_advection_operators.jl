
@inline _advective_tracer_flux_x(args...) = advective_tracer_flux_x(args...)
@inline _advective_tracer_flux_y(args...) = advective_tracer_flux_y(args...)
@inline _advective_tracer_flux_z(args...) = advective_tracer_flux_z(args...)

#####
##### Fallback tracer fluxes!
#####

# Fallback for `nothing` advection
@inline _advective_tracer_flux_x(i, j, k, grid, ::Nothing, args...) = zero(grid)
@inline _advective_tracer_flux_y(i, j, k, grid, ::Nothing, args...) = zero(grid)
@inline _advective_tracer_flux_z(i, j, k, grid, ::Nothing, args...) = zero(grid)

#####
##### Tracer advection operator
#####

"""
    div_uc(i, j, k, grid, advection, U, c)

Calculate the divergence of the flux of a tracer quantity ``c`` being advected by
a velocity field, ``𝛁⋅(𝐯 c)``,

```
1/V * [δxᶜᵃᵃ(Ax * u * ℑxᶠᵃᵃ(c)) + δyᵃᶜᵃ(Ay * v * ℑyᵃᶠᵃ(c)) + δzᵃᵃᶜ(Az * w * ℑzᵃᵃᶠ(c))]
```
which ends up at the location `ccc`.
"""
@inline function div_Uc(i, j, k, grid, advection, U, c)
    return 1/Vᶜᶜᶜ(i, j, k, grid) * (δxᶜᵃᵃ(i, j, k, grid, _advective_tracer_flux_x, advection, U.u, c) +
                                    δyᵃᶜᵃ(i, j, k, grid, _advective_tracer_flux_y, advection, U.v, c) +
                                    δzᵃᵃᶜ(i, j, k, grid, _advective_tracer_flux_z, advection, U.w, c))
end

# Fallbacks for zero velocities, zero tracer and `nothing` advection
@inline div_Uc(i, j, k, grid, advection, ::ZeroU, c) = zero(grid)
@inline div_Uc(i, j, k, grid, advection, U, ::ZeroField) = zero(grid)
@inline div_Uc(i, j, k, grid, advection, ::ZeroU, ::ZeroField) = zero(grid)

@inline div_Uc(i, j, k, grid, ::Nothing, U, c) = zero(grid)
@inline div_Uc(i, j, k, grid, ::Nothing, ::ZeroU, c) = zero(grid)
@inline div_Uc(i, j, k, grid, ::Nothing, U, ::ZeroField) = zero(grid)
@inline div_Uc(i, j, k, grid, ::Nothing, ::ZeroU, ::ZeroField) = zero(grid)

# """
#     c_div_U(i, j, k, grid, advection, U, c)

# Calculate the tracer multiplied by divergence of the velocity field, incase velocity is non-divergent (for example,
# when using PrescribedVelocities) ``c𝛁⋅𝐯``,

# ```
# c[i, j, k] * divᶜᶜᶜ(U)
# ```
# which ends up at the location `ccc`.
# """

# # LB added: Try to implement an extra term incase velocity is non-divergent
# @inline function c_div_U(i, j, k, grid, advection, U, c)
#     return c[i, j, k] * divᶜᶜᶜ(i, j, k, grid, U.u, U.v, U.w)
# end

# # Fallbacks for zero velocities, zero tracer and `nothing` advection
# @inline c_div_U(i, j, k, grid, advection, ::ZeroU, c) = zero(grid)
# @inline c_div_U(i, j, k, grid, advection, U, ::ZeroField) = zero(grid)
# @inline c_div_U(i, j, k, grid, advection, ::ZeroU, ::ZeroField) = zero(grid)

# @inline c_div_U(i, j, k, grid, ::Nothing, U, c) = zero(grid)
# @inline c_div_U(i, j, k, grid, ::Nothing, ::ZeroU, c) = zero(grid)
# @inline c_div_U(i, j, k, grid, ::Nothing, U, ::ZeroField) = zero(grid)
# @inline c_div_U(i, j, k, grid, ::Nothing, ::ZeroU, ::ZeroField) = zero(grid)
