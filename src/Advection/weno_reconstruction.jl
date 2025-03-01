#####
##### Weighted Essentially Non-Oscillatory (WENO) advection scheme
#####

struct WENO{N, FT, XT, YT, ZT, PP, CA, SI} <: AbstractUpwindBiasedAdvectionScheme{N, FT}
    
    "Coefficient for ENO reconstruction on x-faces" 
    coeff_xᶠᵃᵃ::XT
    "Coefficient for ENO reconstruction on x-centers"
    coeff_xᶜᵃᵃ::XT
    "Coefficient for ENO reconstruction on y-faces"
    coeff_yᵃᶠᵃ::YT
    "Coefficient for ENO reconstruction on y-centers"
    coeff_yᵃᶜᵃ::YT
    "Coefficient for ENO reconstruction on z-faces"
    coeff_zᵃᵃᶠ::ZT
    "Coefficient for ENO reconstruction on z-centers"
    coeff_zᵃᵃᶜ::ZT

    "Bounds for maximum-principle-satisfying WENO scheme"
    bounds :: PP

    "Advection scheme used near boundaries"
    buffer_scheme :: CA
    "Reconstruction scheme used for symmetric interpolation"
    advecting_velocity_scheme :: SI

    function WENO{N, FT}(coeff_xᶠᵃᵃ::XT, coeff_xᶜᵃᵃ::XT,
                         coeff_yᵃᶠᵃ::YT, coeff_yᵃᶜᵃ::YT, 
                         coeff_zᵃᵃᶠ::ZT, coeff_zᵃᵃᶜ::ZT,
                         bounds::PP, buffer_scheme::CA,
                         advecting_velocity_scheme :: SI) where {N, FT, XT, YT, ZT, PP, CA, SI}

            return new{N, FT, XT, YT, ZT, PP, CA, SI}(coeff_xᶠᵃᵃ, coeff_xᶜᵃᵃ, 
                                                      coeff_yᵃᶠᵃ, coeff_yᵃᶜᵃ, 
                                                      coeff_zᵃᵃᶠ, coeff_zᵃᵃᶜ,
                                                      bounds, buffer_scheme, advecting_velocity_scheme)
    end
end

"""
    WENO([FT=Float64;] 
         order = 5,
         grid = nothing, 
         bounds = nothing)
               
Construct a weighted essentially non-oscillatory advection scheme of order `order`.

Keyword arguments
=================

- `order`: The order of the WENO advection scheme. Default: 5
- `grid`: (defaults to `nothing`)

Examples
========
```jldoctest
julia> using Oceananigans

julia> WENO()
WENO(order=5)
 Boundary scheme:
    └── WENO(order=3)
 Symmetric scheme:
    └── Centered(order=4)
 Directions:
    ├── X regular
    ├── Y regular
    └── Z regular
```

```jldoctest
julia> using Oceananigans

julia> Nx, Nz = 16, 10;

julia> Lx, Lz = 1e4, 1e3;

julia> chebychev_spaced_z_faces(k) = - Lz/2 - Lz/2 * cos(π * (k - 1) / Nz);

julia> grid = RectilinearGrid(size = (Nx, Nz), halo = (4, 4), topology=(Periodic, Flat, Bounded),
                              x = (0, Lx), z = chebychev_spaced_z_faces);

julia> WENO(grid; order=7)
WENO(order=7)
 Boundary scheme:
    └── WENO(order=5)
 Symmetric scheme:
    └── Centered(order=6)
 Directions:
    ├── X regular
    ├── Y regular
    └── Z stretched
```
"""
function WENO(FT::DataType=Float64; 
              order = 5,
              grid = nothing, 
              bounds = nothing)
    
    if !(grid isa Nothing) 
        FT = eltype(grid)
    end

    mod(order, 2) == 0 && throw(ArgumentError("WENO reconstruction scheme is defined only for odd orders"))

    if !isnothing(bounds)
        @warn "Bounds preserving WENO is experimental."
    end

    if order < 3
        # WENO(order=1) is equivalent to UpwindBiased(order=1)
        return UpwindBiased(FT; order=1)
    else
        N = Int((order + 1) ÷ 2)
        weno_coefficients = compute_reconstruction_coefficients(grid, FT, :WENO; order = N)
        advecting_velocity_scheme = Centered(FT; grid, order = order - 1)
        buffer_scheme = WENO(FT; grid, order=order-2, bounds) 
    end

    return WENO{N, FT}(weno_coefficients..., bounds, buffer_scheme, advecting_velocity_scheme)
end

WENO(grid, FT::DataType=Float64; kwargs...) = WENO(FT; grid, kwargs...)

# Flavours of WENO
const PositiveWENO = WENO{<:Any, <:Any, <:Any, <:Any, <:Any, <:Tuple}

Base.summary(a::WENO{N}) where N = string("WENO(order=", N*2-1, ")")

Base.show(io::IO, a::WENO{N, FT, RX, RY, RZ, PP}) where {N, FT, RX, RY, RZ, PP} =
    print(io, summary(a), " \n",
              a.bounds isa Nothing ? "" : " Bounds : \n    └── $(a.bounds) \n",
              " Boundary scheme: ", "\n",
              "    └── ", summary(a.buffer_scheme) , "\n",
              " Symmetric scheme: ", "\n",
              "    └── ", summary(a.advecting_velocity_scheme) , "\n",
              " Directions:", "\n",
              "    ├── X $(RX == Nothing ? "regular" : "stretched") \n",
              "    ├── Y $(RY == Nothing ? "regular" : "stretched") \n",
              "    └── Z $(RZ == Nothing ? "regular" : "stretched")" )

Adapt.adapt_structure(to, scheme::WENO{N, FT, XT, YT, ZT, PP}) where {N, FT, XT, YT, ZT, PP} =
     WENO{N, FT}(Adapt.adapt(to, scheme.coeff_xᶠᵃᵃ), Adapt.adapt(to, scheme.coeff_xᶜᵃᵃ),
                 Adapt.adapt(to, scheme.coeff_yᵃᶠᵃ), Adapt.adapt(to, scheme.coeff_yᵃᶜᵃ),
                 Adapt.adapt(to, scheme.coeff_zᵃᵃᶠ), Adapt.adapt(to, scheme.coeff_zᵃᵃᶜ),
                 Adapt.adapt(to, scheme.bounds),
                 Adapt.adapt(to, scheme.buffer_scheme),
                 Adapt.adapt(to, scheme.advecting_velocity_scheme))

on_architecture(to, scheme::WENO{N, FT, XT, YT, ZT, PP}) where {N, FT, XT, YT, ZT, PP} =
    WENO{N, FT}(on_architecture(to, scheme.coeff_xᶠᵃᵃ), on_architecture(to, scheme.coeff_xᶜᵃᵃ),
                on_architecture(to, scheme.coeff_yᵃᶠᵃ), on_architecture(to, scheme.coeff_yᵃᶜᵃ),
                on_architecture(to, scheme.coeff_zᵃᵃᶠ), on_architecture(to, scheme.coeff_zᵃᵃᶜ),
                on_architecture(to, scheme.bounds),
                on_architecture(to, scheme.buffer_scheme),
                on_architecture(to, scheme.advecting_velocity_scheme))

# Retrieve precomputed coefficients (+2 for julia's 1 based indices)
@inline retrieve_coeff(scheme::WENO, r, ::Val{1}, i, ::Type{Face})   = @inbounds scheme.coeff_xᶠᵃᵃ[r+2][i] 
@inline retrieve_coeff(scheme::WENO, r, ::Val{1}, i, ::Type{Center}) = @inbounds scheme.coeff_xᶜᵃᵃ[r+2][i] 
@inline retrieve_coeff(scheme::WENO, r, ::Val{2}, i, ::Type{Face})   = @inbounds scheme.coeff_yᵃᶠᵃ[r+2][i] 
@inline retrieve_coeff(scheme::WENO, r, ::Val{2}, i, ::Type{Center}) = @inbounds scheme.coeff_yᵃᶜᵃ[r+2][i] 
@inline retrieve_coeff(scheme::WENO, r, ::Val{3}, i, ::Type{Face})   = @inbounds scheme.coeff_zᵃᵃᶠ[r+2][i] 
@inline retrieve_coeff(scheme::WENO, r, ::Val{3}, i, ::Type{Center}) = @inbounds scheme.coeff_zᵃᵃᶜ[r+2][i] 
