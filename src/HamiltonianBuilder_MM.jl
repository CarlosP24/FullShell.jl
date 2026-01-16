
# Modified Hollow Core Hamiltonian for a multi-mJ model. Valid for a SIS with radii mismatch.
# IMPORTANT NOTICE: THIS CODE USES B INSTEAD OF Φ so LP areas can be different.

"""
    Params_mm

Parameter structure for multi-mJ (multi-angular-momentum) nanowire Hamiltonians.

This structure is designed for building Hamiltonians that explicitly include multiple
angular momentum channels (mJ sites) in a single tight-binding model. Unlike `Params`,
this uses magnetic field B directly rather than normalized flux Φ to allow for different
effective areas in the semiconductor and superconductor regions.

# Fields
## Material and lattice parameters
- `ħ2ome::Float64 = 76.1996`: ℏ²/(2m_e) in meV⋅nm²
- `πoΦ0::Float64 = 1.519257e-3`: π/Φ₀ in units of 1/(T⋅nm²)
- `μBΦ0::Float64 = 119.6941183`: Bohr magneton times magnetic flux quantum
- `m0::Float64 = 0.023`: Effective mass in units of electron mass
- `a0::Float64 = 5`: Lattice constant
- `t::Float64`: Hopping parameter (computed from ħ2ome, m0, a0)
- `echarge::Float64 = 1`: Elementary charge

## Geometry
- `R::Float64 = 70`: Inner radius of the cylinder
- `w::Float64 = 10`: Width of the semiconductor shell 
- `d::Float64 = 10`: Width of the superconductor shell
- `L::Float64 = 0`: Length of the cylinder 

## Multi-mJ specific
- `num_mJ::Int = 5`: Number of angular momentum channels to include
- `hops0::Bool = false`: Whether to include hopping between mJ channels
- `range_hop_m::Int = 3 * num_mJ`: Range for mJ-channel hoppings (in units of a0)

## Physical parameters
- `α::Float64 = 0`: Rashba spin-orbit coupling strength
- `μ::Float64 = 0`: Chemical potential
- `g::Float64 = 0`: g-factor
- `τΓ::Float64 = 1`: Coupling strength to superconductor
- `B::Float64 = 0`: Magnetic field (Tesla)
- `θ::Float64 = 0`: Angle of magnetic field with respect to wire axis

## Superconducting parameters
- `Δ0::ComplexF64 = 0.23`: Bulk superconducting gap 
- `ξd::Float64 = 70`: Superconducting coherence length 
- `shell::String = "Usadel"`: Type of superconductor ("Usadel" or "Ballistic")
- `iω::Float64 = 1e-4`: Small imaginary frequency for regularization

## Legacy/unused parameters
- `σ::Float64 = 0`: Noise parameter (unused, kept for compatibility)
- `preα::Float64 = 0`: Legacy parameter
- `Φ::Float64 = 0`: Legacy flux parameter
- `ishollow::Bool = true`: Legacy parameter
- `Vmax::Float64 = 0`: Legacy dome profile parameter
- `Vmin::Float64 = 0`: Legacy dome profile parameter
- `Vexponent::Float64 = 0`: Legacy dome profile parameter

# Examples
```julia
# Create multi-mJ parameters
p = Params_mm(R=70, num_mJ=5, B=0.1, α=0.2)

# Build multi-mJ Hamiltonian
hSM, hSC, params = build_cyl_mm(p)
```
"""
@with_kw struct Params_mm @deftype Float64
    ħ2ome = 76.1996
    πoΦ0 = 1.519257e-3                    # πoΦ0 = π/\Phi_0 in 1/(T * nm^2). \Phi / \Phi_0 = π / \Phi_0 * B * RLP^2 = πoΦ0 * B * RLP^2
    μBΦ0 = 119.6941183                  #Bohr magneton times magnetic flux quantum
    m0 = 0.023
    a0 = 5
    t = ħ2ome/(2m0*a0^2)
    echarge = 1
    R = 70
    w = 10
    d = 10
    num_mJ = 5
    α = 0
    μ = 0
    g = 0
    τΓ = 1
    B = 0
    θ = 0
    Δ0::ComplexF64 = 0.23
    ξd = 70
    L = 0
    σ = 0                               #Noise parameter, unused here but useful for further constructions
    shell::String = "Usadel"
    iω = 1e-4
    hops0::Bool = false
    range_hop_m = 3 * num_mJ

    # unneccesary here, but needed for legacy code 
    preα = 0
    Φ = 0                               #flux normalized to the flux quantum always
    ishollow::Bool = true
    Vmax = 0                            #dome profile parameters
    Vmin = 0
    Vexponent = 0
end

"""
    build_cyl_mm(; nforced=nothing, phaseshifted=false, kw...)
    build_cyl_mm(p::Params_mm; nforced=nothing, phaseshifted=false)

Build a multi-mJ (multi-angular-momentum) cylindrical nanowire Hamiltonian.

This function constructs a Hamiltonian that explicitly includes multiple angular momentum
channels (mJ) as discrete lattice sites. This approach is useful when inter-channel coupling is needed, such as in a radii-mismatch Josephson junction or to introduce profile disorder.

Unlike `build_cyl`, which uses normalized flux Φ/Φ₀, this function uses the magnetic field
B directly to allow different effective loop areas in different regions.

# Arguments
- `p::Params_mm`: Parameter structure for multi-mJ model
- `nforced::Union{Nothing,Int}`: Force a specific winding number
- `phaseshifted::Bool`: Apply spatially-varying superconducting phase
- `kw...`: Keyword arguments to override `Params_mm` defaults

# Returns
- `hSM`: Parametric Hamiltonian for the semiconductor
- `hSC`: Parametric Hamiltonian with superconducting proximity effect
- `p`: The `Params_mm` structure used

# Physical Model
The lattice includes:
- One spatial dimension (along the wire)
- `2*num_mJ + 1` transverse sites representing different mJ channels
- Empty hopping between mJ channels (controlled by `hops0` and `range_hop_m`)

# Examples
```julia
# Build with default parameters
hSM, hSC, params = build_cyl_mm()

# Build with custom parameters
hSM, hSC, params = build_cyl_mm(R=70, num_mJ=5, B=0.1, α=0.2, μ=0.5)

# Use a Params_mm structure
p = Params_mm(R=70, w=10, d=10, num_mJ=7, B=0.05)
hSM, hSC, params = build_cyl_mm(p)
```
"""
build_cyl_mm(; nforced = nothing, phaseshifted = false, kw...) = build_cyl_mm(Params_mm(; kw...); nforced, phaseshifted)

function build_cyl_mm(p::Params_mm; nforced = nothing, phaseshifted = false)
    @unpack πoΦ0, μBΦ0, a0, t, echarge, R, w, d, num_mJ, α, μ, g, τΓ, B, θ, Δ0, ξd, shell, iω, hops0, range_hop_m  = p

    # Lattice
    # Includes sites along the length of the wire + mJ sites in the transverse direction.

    RLP2 = (R + d/2)^2
    R = floor(R/a0)*a0
    Rav(w) = R - w/2 
    lat = LP.square(; a0) |> supercell((1,0), region = r -> abs(r[2]/a0) <= num_mJ)

    # Model
    # Kinetic term
    # Allow t-hopping only through the length dimension
    ishopz(dr) = iszero(dr[2])
    ishopm(dr) = iszero(dr[1])
    p2 = @onsite((r; μ = μ) -> σ0τz * (2.0 * t - μ)) + hopping((r, dr) -> -t * σ0τz; range = a0, region = (r, dr) -> ishopz(dr))

    # Lienear SOC 
    rashba = @hopping((r, dr; α = α) -> α * (im * dr[1] / (2a0^2)) * σyτz; range = a0, region = (r, dr) -> ishopz(dr))

    # g - Zeeman 
    zeeman = @onsite((; B = B, θ = θ, g = g) -> σzτ0 * 0.5 * g * μBΦ0 * πoΦ0 * B * cos(θ))

    # Magnetic field
    eAφ(B; θ = θ, w = w) = echarge * 0.5 * B * cos(θ) * Rav(w) * πoΦ0
    Φ(B; θ = θ) = B * cos(θ) * RLP2 * πoΦ0
    n(B; θ = θ) = round(Int, Φ(B; θ))
    mJ(r, B; θ = θ) = r[2]/a0 + ifelse(iseven(n(B; θ)), 0.5, 0)

    J(r, B; θ = θ) = mJ(r, B; θ) * σ0τ0 - 0.5 * σzτ0 - 0.5 * n(B; θ) * σ0τz

    gauge = @onsite((r; B = B, α = α, θ = θ, w = w) -> 
        σ0τz * (σ0τz * eAφ(B; θ, w) + J(r, B; θ)/Rav(w))^2 * t * a0^2 -
        σzτz * (σ0τz * eAφ(B; θ, w) + J(r, B; θ)/Rav(w)) * α
    )

    mJ_hop = hopping((r, dr) -> 0*σ0τ0, range = range_hop_m * a0, region = (r, dr) -> ishopm(dr))

    # SM hamiltonian 

    model = p2 + rashba + zeeman + gauge
    if hops0
        model = model + mJ_hop
    end
    hSM = lat |> hamiltonian(model; orbitals = Val(4))

    # Superconductor
    Λ(B, θ) = pairbreaking(Φ(B), n(B; θ), Δ0, ξd, R, d; θ = θ)

    if shell == "Usadel"
        ΣS = ΣS3DUsadel
    elseif shell == "Ballistic"
        ΣS = ΣS3DBallistic
    else
        ΣS = ΣΔ
    end

    ΣS! = @onsite!((o, r; ω = 0, B = B, τΓ = τΓ, θ = θ) ->
            o +  τΓ * ΣS(Δ0, Λ(B, θ), ω);
    )

    # Superconductor phase 
    PhaseShift! = @onsite!((o, r; phase = 0) -> 
        conj(Uphase(phase)) * o * Uphase(phase);
    )

    hSC = hSM |> ΣS!


    if phaseshifted
        hSC = hSC |> PhaseShift!
    end

    return hSM, hSC, p
end

function bandwidth(p::Params_mm)
    @unpack ħ2ome, μ, m0, a0 = p
    return max(abs(4*ħ2ome/(2m0*a0^2) - μ), abs(-4*ħ2ome/(2m0*a0^2) - μ))
end

function get_itip(wire::Params_mm)
    @unpack R, d, πoΦ0, Δ0, ξd, R, d, θ = wire
    RLP2 = (R + d/2)^2
    Φ(B) = B * RLP2 * πoΦ0
    n(B, θ) = round(Int, Φ(B * cos(θ)))
    Λ(B, θ) = pairbreaking(Φ(B), n(B, θ), Δ0, ξd, R, d; θ = θ)
    return (B; θ = θ) -> real(itip(Δ0, Λ(B, θ))) * 0.99
end

function get_Φ(wire::Params_mm)
    @unpack R, d, πoΦ0  = wire
    RLP2 = (R + d/2)^2
    return B -> B * RLP2 * πoΦ0
end

function get_B(wire::Params_mm)
    @unpack R, d, πoΦ0  = wire
    RLP2 = (R + d/2)^2
    return Φ -> Φ / (RLP2 * πoΦ0)
end

function get_Ω(wire::Params_mm)
    @unpack Δ0, ξd, R, d  = wire
    Φ = get_Φ(wire)
    return (B; θ = 0) -> Ω(pairbreaking(Φ(B), round(Int, Φ(B * cos(θ))), Δ0, ξd, R, d; θ = θ), Δ0)
end

"""
    build_harmonic_deformations(wire::Params_mm, harmonics::Dict{Int, Complex})

Constructs a hopping function that incorporates harmonic deformations into the Hamiltonian for a multimode wire.

# Arguments
- `wire::Params_mm`: Structure containing wire parameters-
- `harmonics::Dict{Int, Complex}`: Dictionary mapping harmonic indices to complex deformation amplitudes.

# Returns
- A hopping function with harmonic deformation effects, suitable for use with the `@hopping!` macro.

# Details
The function defines several utilities for Hamiltonian construction, including effective radius calculations, flux, and angular momentum operators. Harmonic deformations are encoded via the `harmonics` dictionary and incorporated into the hopping terms through `k_mixer` and `α_mixer`, which modify the kinetic and spin-orbit contributions, respectively.

The returned hopping function includes these deformation effects and restricts hopping to the appropriate region using the `ishopm` utility.

"""
function build_harmonic_deformations(wire::Params_mm, harmonics::Dict{Int, Complex})
    @unpack R, w, d, πoΦ0, echarge, a0, t, θ = wire

    # Hamiltonian utilities
    RLP2 = (R + d/2)^2
    Rav(w) = R - w/2
    eAφ(B; θ = θ, w = w) = echarge * 0.5 * B * cos(θ) * Rav(w) * π * πoΦ0
    Φ(B; θ = θ) = B * cos(θ) * RLP2 * πoΦ0
    n(B; θ = θ) = round(Int, Φ(B * cos(θ)))
    mJ(r, B; θ = θ) = r[2]/a0 + ifelse(iseven(n(B; θ)), 0.5, 0)
    J(r, B; θ = θ) = mJ(r, B; θ = θ) * σ0τ0 - 0.5 * σzτ0 - 0.5 * n(B; θ) * σ0τz

    # Mode mixing utilities
    ℓ(dr) = round(Int, dr[2]/a0 |> abs)
    ishopm(dr) = iszero(dr[1])
    δRp(dr) = get(harmonics, ℓ(dr), 0)
    δR(dr) = ifelse(dr[2] > 0, δRp(dr), conj(δRp(dr)))

    # Mode mixers
    k_mixer(r, dr, B; θ = θ, w = w) = - t * a0^2 * δR(dr) * (J(r, B; θ)^2 / Rav(w)^2 + 0.25 * ℓ(dr)^2 * σ0τ0 / Rav(w)^2 - eAφ(B; θ, w)^2 * σ0τ0) * σ0τz
    α_mixer(r, dr, B; θ = θ, w = w) = - (α / 2) * δR(dr) * (J(r, B; θ) / Rav(w) - eAφ(B; θ, w) * σ0τz) * σzτz

    return @hopping!((t, r, dr; B = 0, θ = θ, w = w) ->
        t + k_mixer(r, dr, B; θ, w) + α_mixer(r, dr, B; θ, w),
        range = 2 * a0 * length(harmonics), region = (r, dr) -> ishopm(dr)
    )
end