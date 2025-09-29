
# Modified Hollow Core Hamiltonian for a multi-mJ model. Valid for a SIS with radii mismatch.
# IMPORTANT NOTICE: THIS CODE USES B INSTEAD OF Φ so LP areas can be different.

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
    zeeman = @onsite((; B = B, θ = θ) -> σzτ0 * 0.5 * g * μBΦ0 * πoΦ0 * B * cos(θ))

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