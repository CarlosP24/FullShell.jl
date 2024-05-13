
@with_kw struct Params @deftype Float64 #Units: nm, meV
    ħ2ome = 76.1996
    μBΦ0 = 119.6941183                  #Bohr magneton times magnetic flux quantum
    m0 = 0.023
    g = 12
    P = 919.7                           #Bulck SOC parameters
    Δg = 417  
    Δs = 390  
    preα = P^2/3 * (1/Δg^2 - 1/(Δg + Δs)^2)  #α = preα \times Rashba things
    a0 = 5
    t = ħ2ome/(2m0*a0^2)
    echarge = 1
    R = 70                              #radius of the cylinder
    w = 10                              #width of the semiconductor
    d = 10                              #wdith of the superconductor
    Vmax = 0                            #dome profile parameters
    Vmin = 0
    Vexponent = 2
    Δ0::ComplexF64 = 0.2 
    ξd = 70 
    α = 0                               #parameter defaults
    μ = 0
    τΓ = 1
    Φ = 1                               #flux normalized to the flux quantum always
    ishollow::Bool = true
    L = 100                             #length of the cylinder
end

# Hamiltonian constructor 

ΣS3DUsadel(Δ0, Λ, ω;) = -(uUsadel(Δ0, Λ, ω) * σ0τ0 - σ0τx) / sqrt(complex(1-uUsadel(Δ0, Λ, ω)^2))

build_cyl(; nforced = nothing, kw...) = build_cyl(Params(; kw...); nforced,)

function build_cyl(p::Params; nforced = nothing)
    @unpack μBΦ0, m0, g, preα, a0, t, echarge, R, w, d, Vmax, Vmin, Vexponent, Δ0, ξd, α, μ, τΓ, Φ, ishollow = p 

    # Lattice

    R = floor(R/a0)*a0
    lat = if ishollow
            Rav = R - w/2
            LP.square(; a0) |> supercell((1, 0)) |> Quantica.transform!(r -> r + SA[0, Rav])
          else
            LP.square(; a0) |> supercell((1, 0), region = r -> max(a0, R - w) <= r[2] <= R)
          end

    # Model

    # Kinetic term
    p2 = @onsite((r; μ = μ) -> σ0τz *(t * ifelse(r[2] ≈ a0, 2.0 + 1.5, 2.0 + 2.0*!ishollow) - μ)) + hopping((r, dr) -> -t * σ0τz * ifelse(iszero(dr[1]), r[2]/sqrt(r[2]^2 - 0.25*dr[2]^2), 1); range = a0)

    # Dome profile
    V(ρ, v0, v1) = v0 + (v1 - v0) * (ρ/R)^Vexponent
    dϕ(ρ, v0, v1) = - (Vexponent/R) * (v1 - v0) * (ρ/R)^(Vexponent - 1) # ϕ = -V
    potential = @onsite((r; Vmax = Vmax, Vmin = Vmin) -> σ0τz * V(r[2], Vmax, Vmin))

    # Linear SOC
    rashba = @hopping((r, dr; α = α, preα = preα, Vmin = Vmin, Vmax = Vmax) -> (α + preα * dϕ(r[2], Vmax, Vmin)) * (im * dr[1] / (2 * a0^2)) * σyτz; range = a0)

    # g - Zeeman
    zeeman = @onsite((; Φ = Φ) -> σzτ0 * 0.5 * g * μBΦ0 * Φ / area_LP)

    # Magnetic field
    area_LP = π * (R + d/2)^2 
    eAφ(r, Φ) = echarge * 0.5 * π * Φ * r[2] / area_LP
    n(Φ) = ifelse(nforced === nothing, round(Int, Φ), nforced)
    mJ(Z, Φ) = Z + ifelse(iseven(n(Φ)), 0.5, 0.0)
    J(Z, Φ) = mJ(Z, Φ) * σ0τ0 - 0.5 * σzτ0 - 0.5 * n(Φ) * σ0τz 
    gauge = @onsite((r; Φ = Φ, Z = 0, α = α, preα = preα, Vmax = Vmax, Vmin = Vmin) ->
          σ0τz * (σ0τz * eAφ(r, Φ) + J(Z, Φ) / r[2])^2 * t * a0^2 -
          σzτz * (σ0τz * eAφ(r, Φ) + J(Z, Φ) / r[2]) * (α + preα * dϕ(r[2], Vmax, Vmin)) 
    ) 

    # SM hamiltonian 

    hSM = lat |> hamiltonian(p2 + potential + rashba + zeeman + gauge; orbitals = Val(4))

    # Superconductor
    Λ(Φ) = pairbreaking(Φ, n(Φ), Δ0, ξd, R, d)
    ΦLP(Φ) = LP_lobe(n(Φ), ξd, R, d)
    ΣS! = @onsite!((o, r; ω = 0, Φ = Φ, τΓ = τΓ) ->
          o +  τΓ * Δ0 * ifelse(is_in_lobe(Φ, ΦLP(Φ)...), ΣS3DUsadel(Δ0, Λ(Φ), ω), 1im * σ0τ0);
          region = ishollow ? Returns(true) : r -> r[2] > R - a0/2
    )

    hSC = hSM |> ΣS!

    return hSM, hSC, p
end