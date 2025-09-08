
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

    # unneccesary here, but needed for legacy code
    num_mJ = 5
    iω = 1e-5
    Zs::Union{UnitRange, Vector{Int}, Int} = -5:5
    πoΦ0 = 1.519267e-3                  # πoΦ0 = π/\Phi_0 in 1/(T * nm^2). \Phi / \Phi_0 = π / \Phi_0 * B * RLP^2 = πoΦ0 * B * RLP^2
    hops0::Bool = false
    range_hop_m = 0
end

# Hamiltonian constructor 
ΣS3DUsadel(Δ0, Λ, ω;) = - Δ0 *(uUsadel(Δ0, Λ, ω) * σ0τ0 - σ0τx) / sqrt(complex(1-uUsadel(Δ0, Λ, ω)^2))
ΣS3DUsadel_old(Δ0, Λ, ω;) = - Δ0 *(uUsadel_old(Δ0, Λ, ω) * σ0τ0 - σ0τx) / sqrt(complex(1-uUsadel_old(Δ0, Λ, ω)^2))
ΣS3DBallistic(Δ0, Λ, ω;) = - Δ0 *((ω/Ω(Λ, Δ0)) * σ0τ0 - σ0τx) / sqrt(complex(1-(ω/Ω(Λ, Δ0))^2))
#ΣS3DUsadel(Δ0, Λ, ω;) = - Δ0 *(usimple(Δ0, Λ, ω) * σ0τ0 - σ0τx) / sqrt(complex(1-usimple(Δ0, Λ, ω)^2))
ΣΔ(Δ0, Λ, ω;) = (ΔD(Λ, Δ0)^(2/3) - Λ^(2/3))^(3/2) * σ0τx
Uphase(phase) = exp(im * phase * σ0τz /2)

build_cyl(; nforced = nothing, phaseshifted = false, kw...) = build_cyl(Params(; kw...); nforced, phaseshifted)

function build_cyl(p::Params; nforced = nothing, phaseshifted = false)
    @unpack μBΦ0, m0, g, preα, a0, t, echarge, R, w, d, Vmax, Vmin, Vexponent, Δ0, ξd, α, μ, τΓ, Φ, θ, Z, ishollow, shell = p

    # Lattice
    RLP2 = (R + d/2)^2
    R = floor(R/a0)*a0
    area_LP = π * RLP2 

    lat = if ishollow
      
            Rav = R - w/2
            LP.square(; a0) |> supercell((1, 0)) |> Quantica.transform!(r -> r + SA[0, Rav])
          else
            LP.square(; a0) |> supercell((1, 0), region = r -> max(a0, R - w) <= r[2] <= R)
          end

    # Model

    # Kinetic term
    p2 = @onsite((r; μ = μ) -> σ0τz *(t * ifelse(r[2] ≈ a0, 2.0 + 1.5, 2.0 + 2.0*!ishollow)  - μ)) + hopping((r, dr) -> -t * σ0τz * ifelse(iszero(dr[1]), r[2]/sqrt(r[2]^2 - 0.25*dr[2]^2), 1); range = a0)

    # Dome profile
    V(ρ, v0, v1) = v0 + (v1 - v0) * (ρ/R)^Vexponent
    dϕ(ρ, v0, v1) = - (Vexponent/R) * (v1 - v0) * (ρ/R)^(Vexponent - 1) # ϕ = -V
    potential = @onsite((r; Vmax = Vmax, Vmin = Vmin, ) -> σ0τz * V(r[2], Vmax, Vmin))

    # Linear SOC
    rashba = @hopping((r, dr; α = α, preα = preα, Vmin = Vmin, Vmax = Vmax) -> (α + preα * dϕ(r[2], Vmax, Vmin)) * (im * dr[1] / (2 * a0^2)) * σyτz; range = a0)

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
          o +  τΓ * ΣS(Δ0, Λ(Φ, θ), ω);
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

    return hSM, hSC, p
end

function bandwidth(p::Params)
    @unpack ħ2ome, μ, m0, a0 = p
    return max(abs(4*ħ2ome/(2m0*a0^2) - μ), abs(-4*ħ2ome/(2m0*a0^2) - μ))
end

function get_itip(wire::Params)
  @unpack R, d, Δ0, ξd, R, d, θ  = wire
  n(Φ; θ = θ) = round(Int, Φ * cos(θ))
  Λ(Φ; θ = θ) = pairbreaking(Φ, n(Φ; θ), Δ0, ξd, R, d; θ = θ)
  return Φ -> real(itip(Δ0, Λ(Φ; θ))) * 0.99
end

function get_Φ(wire::Params)
  @unpack R, d, πoΦ0  = wire
  RLP2 = (R + d/2)^2
  return B -> B * RLP2 * πoΦ0
end

function get_B(wire::Params)
  @unpack R, d, πoΦ0 = wire
  RLP2 = (R + d/2)^2
  return Φ -> Φ / (RLP2 * πoΦ0)
end

function get_Ω(wire::Params)
  @unpack Δ0, ξd, R, d, θ = wire
  return (Φ; θ = θ) -> Ω(pairbreaking(Φ, round(Int, Φ * cos(θ)), Δ0, ξd, R, d; θ), Δ0)
end