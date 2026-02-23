
"""
    Params

Parameter structure for full-shell nanowire Hamiltonians.

This structure contains all physical parameters needed to construct a cylindrical
semiconductor-superconductor nanowire Hamiltonian. All quantities are in units of
meV and nm unless otherwise specified.

# Fields
## Material parameters
- `ħ2ome::Float64 = 76.1996`: ℏ²/(2m_e) in meV⋅nm²
- `μBΦ0::Float64 = 119.6941183`: Bohr magneton times magnetic flux quantum
- `m0::Float64 = 0.023`: Effective mass in units of electron mass
- `g::Float64 = 12`: g-factor
- `P::Float64 = 919.7`: Kane parameter
- `Δg::Float64 = 417`: Semiconductor bandgap
- `Δs::Float64 = 390`: Split-off gap
- `preα::Float64`: Rashba prefactor computed from bulk parameters

## Lattice parameters
- `a0::Float64 = 5`: Lattice constant
- `t::Float64`: Hopping parameter, computed from ħ2ome, m0, and a0
- `R::Float64 = 70`: Inner radius of the cylinder 
- `w::Float64 = 10`: Width of the semiconductor shell 
- `d::Float64 = 10`: Width of the superconductor shell 
- `L::Float64 = 100`: Length of the cylinder
- `ishollow::Bool = true`: Whether the cylinder is hollow (true) or has filled interior

## Electrostatic profile
- `Vmax::Float64 = 0`: Maximum electrostatic potential 
- `Vmin::Float64`: Minimum electrostatic potential, defaults to Vmax
- `Vexponent::Float64 = 2`: Exponent for dome-shaped potential profile

## Superconducting parameters
- `Δ0::ComplexF64 = 0.2`: Bulk superconducting gap 
- `ξd::Float64 = 70`: Superconducting coherence length 
- `shell::String = "Usadel"`: Type of superconductor ("Usadel", "Ballistic", or other)

## Tunable parameters
- `α::Float64 = 0`: Rashba spin-orbit coupling strength (field independent)
- `μ::Float64 = 0`: Chemical potential
- `τΓ::Float64 = 1`: semiconductor-superconductor couplqing strength
- `Φ::Float64 = 1`: Magnetic flux through cylinder in units of flux quantum Φ₀
- `θ::Float64 = 0`: Angle of magnetic field with respect to wire axis
- `Z::Int = 0`: Quantum number related to angular momentum

## Advanced parameters
- `echarge::Float64 = 1`: Elementary charge (in natural units)
- `num_mJ::Int = 5`: Number of angular momentum channels (for multi-mJ models)
- `iω::Float64 = 1e-5`: Small imaginary frequency for regularization
- `Zs::Union{UnitRange, Vector{Int}, Int} = -5:5`: Range of Z quantum numbers
- `πoΦ0::Float64 = 1.519267e-3`: π/Φ₀ in units of 1/(T⋅nm²)
- `hops0::Bool = false`: Whether to include zero-range hoppings
- `range_hop_m::Int = 0`: Range for angular momentum hoppings

# Examples
```julia
# Create parameters with defaults
p = Params()

# Create parameters with custom values
p = Params(R = 100, μ = 0.5, α = 0.2, Δ0 = 0.25)

# Use with build_cyl
hSM, hSC, params = build_cyl(p)
```
"""
@with_kw struct Params @deftype Float64 #Units: nm, meV
    ħ2ome = 76.1996
    μBΦ0 = 119.6941183                  #Bohr magneton times magnetic flux quantum
    m0 = 0.023
    g = 12
    P = 919.7                           #Bulk SOC parameters
    Δg = 417  
    Δs = 390  
    preα = P^2/3 * (1/Δg^2 - 1/(Δg + Δs)^2)  #α = preα \times Rashba things
    a0 = 5
    az = a0
    t = ħ2ome/(2m0*a0^2)
    echarge = 1
    R = 70                              #radius of the cylinder
    w = 10                              #width of the semiconductor
    d = 10                              #width of the superconductor
    Vmax = 0                            #dome profile parameters
    Vmin = Vmax
    Vexponent = 2
    Δ0::ComplexF64 = 0.2 
    ξd = 70 
    α = 0                               #parameter defaults
    μ = 0
    τΓ = 1
    Φ = 1                               #flux normalized to the flux quantum always
    θ = 0                               #angle of the flux
    ishollow::Bool = true
    L = 100                             #length of the cylinder
    shell::String = "Usadel"
    Z = 0
    bandbottom::Bool = false            # whether to shift the energy zero to the band bottom 

    # unneccesary here, but needed for legacy code
    num_mJ = 5
    iω = 1e-5
    Zs::Union{UnitRange, Vector{Int}, Int} = -5:5
    πoΦ0 = 1.519267e-3                  # πoΦ0 = π/\Phi_0 in 1/(T * nm^2). \Phi / \Phi_0 = π / \Phi_0 * B * RLP^2 = πoΦ0 * B * RLP^2
    hops0::Bool = false
    range_hop_m = 0
end

# Hamiltonian constructor 
ΣS3DUsadel(Δ0, Λ, ω;) = - Δ0 *(uUsadel(Δ0, Λ, ω) * σ0τ0 - σ0τx) / sqrt(complex(1-uUsadel(Δ0, Λ, ω)^2))  # This self-energy respects causality for Re[ω] < 0.
ΣS3DUsadel_old(Δ0, Λ, ω;) = - Δ0 *(uUsadel_old(Δ0, Λ, ω) * σ0τ0 - σ0τx) / sqrt(complex(1-uUsadel_old(Δ0, Λ, ω)^2))
ΣS3DBallistic(Δ0, Λ, ω;) = - Δ0 *((ω/Ω(Λ, Δ0)) * σ0τ0 - σ0τx) / sqrt(complex(1-(ω/Ω(Λ, Δ0))^2))   # This self-energy respects causality for all ω.
#ΣS3DUsadel(Δ0, Λ, ω;) = - Δ0 *(usimple(Δ0, Λ, ω) * σ0τ0 - σ0τx) / sqrt(complex(1-usimple(Δ0, Λ, ω)^2))
ΣΔ(Δ0, Λ, ω;) = (ΔD(Λ, Δ0)^(2/3) - Λ^(2/3))^(3/2) * σ0τx
Uphase(phase) = exp(im * phase * σ0τz /2)

"""
    build_cyl(; nforced=nothing, phaseshifted=false, kw...)
    build_cyl(p::Params; nforced=nothing, phaseshifted=false)

Build a full-shell cylindrical nanowire Hamiltonian with proximity-induced superconductivity.

This function constructs two Hamiltonians for a semiconductor nanowire with a superconducting
shell: one for the normal state (`hSM`) and one with the superconducting proximity effect (`hSC`).
The model includes kinetic energy, Rashba spin-orbit coupling, Zeeman coupling, orbital magnetic
field effects, and a proximity-induced pairing potential from the superconductor.

# Arguments
- `p::Params`: Parameter structure containing all physical parameters
- `nforced::Union{Nothing,Int}`: If specified, forces a particular winding number instead of computing from flux
- `phaseshifted::Bool`: If true, applies a spatially-varying superconducting phase
- `kw...`: Keyword arguments to override default `Params` values

# Returns
- `hSM`: Parametric Hamiltonian for the semiconductor in normal state
- `hSC`: Parametric Hamiltonian with superconducting proximity effect
- `p`: The `Params` structure used (useful when using keyword arguments)

# Physical Model
The Hamiltonian includes:
1. **Kinetic term**: Tight-binding with effective mass `m0`
2. **Electrostatic potential**: Dome-shaped profile V(ρ) = Vmax + (Vmin - Vmax)(ρ/R)^Vexponent
3. **Rashba SOC**: Linear in momentum with strength α + preα⋅∂V/∂ρ
4. **Zeeman coupling**: g-factor coupling to magnetic field
5. **Orbital magnetic field**: Vector potential in cylindrical gauge
6. **Proximity effect**: Self-energy from superconductor via Usadel or ballistic formalism

# Examples
```julia
# Build with default parameters
hSM, hSC, params = build_cyl()

# Build with custom parameters
hSM, hSC, params = build_cyl(R=100, μ=0.5, α=0.2, Δ0=0.25, Φ=0.5)

# Use a Params structure
p = Params(R=70, w=10, d=10, Δ0=0.23, μ=0.75)
hSM, hSC, params = build_cyl(p)

# Force a specific winding number
hSM, hSC, params = build_cyl(R=70, Φ=1.5, nforced=1)
```

# Notes
- Default units are meV for energy and nm for length
- The returned Hamiltonians are parametric and can be evaluated at different parameter values
- The superconducting self-energy type is controlled by the `shell` parameter:
  - "Usadel": Diffusive superconductor (default)
  - "Ballistic": Clean superconductor
  - Other: Simplified gap suppression model
"""
build_cyl(; nforced = nothing, phaseshifted = false, kw...) = build_cyl(Params(; kw...); nforced, phaseshifted)

function build_cyl(p::Params; nforced = nothing, phaseshifted = false)
    @unpack μBΦ0, ħ2ome, m0, g, preα, a0, az, echarge, R, w, d, Vmax, Vmin, Vexponent, Δ0, ξd, α, μ, τΓ, Φ, θ, Z, ishollow, shell, bandbottom = p
    
    t = ħ2ome/(2m0*a0^2)
    tz = ħ2ome/(2m0*az^2)

    # Lattice
    RLP2 = (R + d/2)^2
    R = floor(R/a0)*a0
    area_LP = π * RLP2 

    A0 = max(a0, az)

    lat = if ishollow
            Rav = R - w/2
            LP.square(; a0) |> supercell((1, 0)) |> Quantica.transform!(r -> r + SA[0, Rav])
          else
            lattice(sublat((0., 0.)); bravais = SA[az 0.; 0. a0]') |> supercell((1, 0), region = r -> max(a0, R - w) <= r[2] <= R)
          end

    # Model

    # Kinetic term (DLL-FDM discretization)
    p2 = if ishollow
      @onsite((r; μ = μ) -> σ0τz * (tz * 2.0 - μ)) + hopping((r, dr) -> -tz * σ0τz; range = A0)
    else
      @onsite((r; μ = μ) -> σ0τz * (ifelse(r[2] ≈ a0, 2.0 * tz + 1.5 * t, 2.0 * tz + 2.0 * t) - μ)) + hopping((r, dr) -> -tz * σ0τz * iszero(dr[2]) - t * σ0τz * iszero(dr[1]) * r[2]/sqrt(r[2]^2 - 0.25*dr[2]^2); range = A0)
    end
    # p2 = if ishollow
    #   @onsite((r; μ = μ) -> σ0τz * (tz * 2.0 - μ)) + hopping((r, dr) -> -tz * σ0τz; range = A0)
    # else
    #   @onsite((r; μ = μ) -> σ0τz * (ifelse(r[2] ≈ a0, 0.5 * t, 0.0) - μ)) + hopping((r, dr) -> -tz * σ0τz * iszero(dr[2]) - t * σ0τz * iszero(dr[1]) * r[2]/sqrt(r[2]^2 - 0.25*dr[2]^2); range = A0)
    # end
    # Dome profile
    V(ρ, v0, v1) = v0 + (v1 - v0) * (ρ/R)^Vexponent
    dϕ(ρ, v0, v1) = - (Vexponent/R) * (v1 - v0) * (ρ/R)^(Vexponent - 1) # ϕ = -V
    potential = @onsite((r; Vmax = Vmax, Vmin = Vmin, ) -> σ0τz * V(r[2], Vmax, Vmin))

    # Linear SOC
    rashba = @hopping((r, dr; α = α, preα = preα, Vmin = Vmin, Vmax = Vmax) -> (α + preα * dϕ(r[2], Vmax, Vmin)) * (im * dr[1] / (2 * az^2)) * σyτz; range = A0)

    # g - Zeeman
    zeeman = @onsite((; Φ = Φ, θ = θ) -> σzτ0 * 0.5 * g * μBΦ0 * Φ * cos(θ) / area_LP)

    # Magnetic field
    eAφ(r, Φ; θ = θ) = echarge * 0.5 * π * Φ * r[2] * cos(θ) / area_LP
    n(Φ; θ = θ) = ifelse(nforced === nothing, round(Int, Φ * cos(θ)), nforced)
    mJ(Z, Φ; θ = θ) = Z + ifelse(iseven(n(Φ; θ)), 0.5, 0.0)
    J(Z, Φ; θ = θ) = mJ(Z, Φ; θ) * σ0τ0 - 0.5 * σzτ0 - 0.5 * n(Φ; θ) * σ0τz
    gauge = @onsite((r; Φ = Φ, θ = θ, Z = Z, α = α, preα = preα, Vmax = Vmax, Vmin = Vmin) ->
          σ0τz * (σ0τz * eAφ(r, Φ; θ) + J(Z, Φ; θ) / r[2])^2 * t * a0^2 -
          σzτz * (σ0τz * eAφ(r, Φ; θ) + J(Z, Φ; θ) / r[2]) * (α + preα * dϕ(r[2], Vmax, Vmin))
    )

    # SM hamiltonian 

    hSM = lat |> hamiltonian(p2 + potential + rashba + zeeman + gauge; orbitals = Val(4))

    # Superconductor
    Λ(Φ, θ) = pairbreaking(Φ, n(Φ; θ), Δ0, ξd, R, d; θ = θ)

    if shell == "Usadel"
      ΣS = ΣS3DUsadel
    elseif shell == "Ballistic"
      ΣS = ΣS3DBallistic
    elseif shell == "Usadel_old"
      ΣS = ΣS3DUsadel_old
    else
      ΣS = ΣΔ
    end

    ΣS! = @onsite!((o, r; ω = 0, Φ = Φ, τΓ = τΓ, θ = θ) ->
          o +  τΓ * a0 / (2π * R) *  ΣS(Δ0, Λ(Φ, θ), ω);
          region = ishollow ? Returns(true) : r -> r[2] > R - a0/2
    )

    # Superconductor phase 
    PhaseShift! = @onsite!((o, r; phase = 0) -> 
                  conj(Uphase(phase)) * o * Uphase(phase);
                  region = ishollow ? Returns(true) : r -> r[2] > R - a0/2
    )

    hSC = hSM |> ΣS!

    if phaseshifted
      hSC = hSC |> PhaseShift!
    end

    if bandbottom
      h0_SM = hSM(; μ = 0, Φ = 0, preα = 0, α = 0)[] |> Array
      E_bottom = minimum(real(eigvals(h0_SM)))  # lowest eigenvalue of SM, not abs!
      E_bottom! = @onsite!((o, r;) -> o - E_bottom * σ0τz; region = Returns(true))
      hSM = hSM |> E_bottom!
      hSC = hSC |> E_bottom!
    end

    return hSM, hSC, p
end

"""
    bandwidth(p::Params)

Calculate the bandwidth of the semiconductor nanowire.

Returns the maximum energy difference between the chemical potential and the band edges,
considering the tight-binding dispersion.

# Arguments
- `p::Params`: Parameter structure

# Returns
- Bandwidth in meV
"""
function bandwidth(p::Params)
    @unpack ħ2ome, μ, m0, a0 = p
    return max(abs(4*ħ2ome/(2m0*a0^2) - μ), abs(-4*ħ2ome/(2m0*a0^2) - μ))
end

"""
    get_itip(wire::Params)

Get a function that returns the imaginary frequency threshold as a function of flux.

This function returns a closure that computes 0.99×itip(Φ), where itip is the
imaginary frequency above which the Usadel solution becomes unphysical. This is useful
for determining valid integration paths in the complex frequency plane.

# Arguments
- `wire::Params`: Parameter structure

# Returns
- A function `Φ -> iω_threshold` that maps flux to the imaginary frequency threshold
"""
function get_itip(wire::Params)
  @unpack R, d, Δ0, ξd, R, d, θ  = wire
  n(Φ; θ = θ) = round(Int, Φ * cos(θ))
  Λ(Φ; θ = θ) = pairbreaking(Φ, n(Φ; θ), Δ0, ξd, R, d; θ = θ)
  return Φ -> real(itip(Δ0, Λ(Φ; θ))) * 0.99
end

"""
    get_Φ(wire::Params)

Get a function to convert magnetic field to flux.

Returns a closure that converts a magnetic field B (in Tesla) to the corresponding
flux through the nanowire in units of the flux quantum Φ₀.

# Arguments
- `wire::Params`: Parameter structure

# Returns
- A function `B -> Φ/Φ₀` that maps magnetic field to normalized flux
"""
function get_Φ(wire::Params)
  @unpack R, d, πoΦ0  = wire
  RLP2 = (R + d/2)^2
  return B -> B * RLP2 * πoΦ0
end

"""
    get_B(wire::Params)

Get a function to convert flux to magnetic field.

Returns a closure that converts flux (in units of flux quantum Φ₀) to the corresponding
magnetic field B (in Tesla).

# Arguments
- `wire::Params`: Parameter structure

# Returns
- A function `Φ/Φ₀ -> B` that maps normalized flux to magnetic field
"""
function get_B(wire::Params)
  @unpack R, d, πoΦ0 = wire
  RLP2 = (R + d/2)^2
  return Φ -> Φ / (RLP2 * πoΦ0)
end

"""
    get_Ω(wire::Params)

Get a function that returns the effective gap Ω as a function of flux.

The effective gap Ω(Λ, Δ₀) accounts for the suppression of the superconducting gap
due to magnetic field and pair-breaking effects.

# Arguments
- `wire::Params`: Parameter structure

# Returns
- A function `Φ -> Ω` that maps flux to the effective superconducting gap in meV
"""
function get_Ω(wire::Params)
  @unpack Δ0, ξd, R, d, θ = wire
  return (Φ; θ = θ) -> Ω(pairbreaking(Φ, round(Int, Φ * cos(θ)), Δ0, ξd, R, d; θ), Δ0)
end