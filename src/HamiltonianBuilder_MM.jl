
# Modified Hollow Core Hamiltonian for a multi-mJ model. Valid for a SIS with radii mismatch.
# IMPORTANT NOTICE: THIS CODE USES B INSTEAD OF Φ so LP areas can be different.

@with_kw struct Params_mm @deftype Float64
    ħ2ome = 76.1996
    conv = 1.5193e-3 # Magnetic field in T to flux prefactor
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
    Δ0::ComplexF64 = 0.23
    ξd = 70
    L = 0
    σ = 0                               #Noise parameter, unused here but useful for further constructions
    Usadel = true

    # unneccesary here, but needed for legacy code 
    preα = 0
    Φ = 0                               #flux normalized to the flux quantum always
    ishollow::Bool = true
    Vmax = 0                            #dome profile parameters
    Vmin = 0
    Vexponent = 0
end

build_cyl_mm(; nforced = nothing, phaseshifted = false, kw...) = build_cyl_mm(Params_mm(; kw...); nforced, phaseshifted)

function build_cyl_mm(p::Params_mm; nforced = nothing, phaseshifted = false)
    @unpack conv,μBΦ0, a0, t, echarge, R, w, d, num_mJ, α, μ, g, τΓ, B, Δ0, ξd, Usadel  = p

    # Lattice
    # Includes sites along the length of the wire + mJ sites in the transverse direction.

    R = floor(R/a0)*a0
    Rav = R - w/2 
    lat = LP.square(; a0) |> supercell((1,0), region = r -> abs(r[2]/a0) <= num_mJ)

    # Model
    # Kinetic term
    # Allow t-hopping only through the length dimension
    ishopz(dr) = iszero(dr[2])
    p2 = @onsite((; μ = μ) -> σ0τz * (2.0 * t - μ)) + hopping((r, dr) -> -t * σ0τz; range = a0, region = (r, dr) -> ishopz(dr))

    # Lienear SOC 
    rashba = @hopping((r, dr; α = α) -> α * (im * dr[1] / (2a0^2)) * σyτz; range = a0, region = (r, dr) -> ishopz(dr))

    # g - Zeeman 
    zeeman = @onsite((; B = B) -> σzτ0 * 0.5 * g * μBΦ0 * conv * B)

    # Magnetic field
    area_LP = π * (R + d/2)^2
    eAφ(B) = echarge * 0.5 * B * Rav * π * conv
    Φ(B) = B * area_LP * conv
    n(B) = round(Int, Φ(B))
    mJ(r, B) = r[2]/a0 + ifelse(iseven(n(B)), 0.5, 0)
    J(r, B) = mJ(r, B)*σ0τ0 - 0.5*σzτ0 - 0.5*n(B)*σ0τz

    gauge = @onsite((r; B = B, α = α) -> 
        σ0τz * (σ0τz * eAφ(B) + J(r, B)/Rav)^2 * t * a0^2 -
        σzτz * (σ0τz * eAφ(B) + J(r, B)/Rav) * α
    )

    # SM hamiltonian 

    hSM = lat |> hamiltonian(p2 + rashba + zeeman + gauge; orbitals = Val(4))

    # Superconductor
    Λ(B) = pairbreaking(Φ(B), n(B), Δ0, ξd, R, d)

    if Usadel 
        ΣS = ΣS3DUsadel
      else
        ΣS = ΣΔ
      end

    ΣS! = @onsite!((o, r; ω = 0, B = B, τΓ = τΓ) ->
            o +  τΓ * ΣS(Δ0, Λ(B), ω);
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