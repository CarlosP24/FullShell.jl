using Quantica, FullShell
using LinearAlgebra, Parameters

p = Params(; ishollow = false, R = 70, w = 70, a0 = 2)
@unpack μBΦ0, ħ2ome, m0, g, preα, a0, az, echarge, R, w, d, Vmax, Vmin, Vexponent, Δ0, ξd, α, μ, τΓ, Φ, θ, Z, ishollow, shell, bandbottom = p

nforced = nothing
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

    sp = spectrum(hSM(; μ = 0, Φ = 0, preα = 0), 0)

    ϵs, ψs = sp 
    ϵs .|> real .|> abs |> minimum