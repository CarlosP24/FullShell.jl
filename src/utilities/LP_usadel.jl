include("AbrikosovSolver.jl")

"""
    pairbreaking(Φ, n, Δ0, ξd, R, d; θ=0)

Calculate the pair-breaking parameter Λ due to magnetic field in a cylindrical nanowire.

The pair-breaking parameter quantifies the suppression of superconductivity due to the
magnetic field. It accounts for both parallel and perpendicular components of the field
and includes geometric effects from the cylindrical shell geometry.

# Arguments
- `Φ::Real`: Magnetic flux through the cylinder in units of flux quantum Φ₀
- `n::Int`: Winding number (integer closest to Φ⋅cos(θ))
- `Δ0::Number`: Bulk superconducting gap (meV)
- `ξd::Real`: Superconducting coherence length (nm)
- `R::Real`: Inner radius of the cylinder (nm)
- `d::Real`: Width of the superconductor shell (nm)
- `θ::Real = 0`: Angle between magnetic field and wire axis (radians)

# Returns
- `Λ::Real`: Pair-breaking parameter in meV, guaranteed to be at least Δ₀×10⁻³

# Formula
The pair-breaking parameter is computed as:
```
Λ_∥ = 4(Φ_∥ - n)² + (d²/R_LP²)(Φ_∥² + n²/3)
Λ_⊥ = 4Φ_⊥²
Λ = (ξ_d² Δ₀)/(1.76π R_LP²) (Λ_∥ + Λ_⊥)
```
where R_LP = R + d/2, Φ_∥ = Φcos(θ), and Φ_⊥ = Φsin(θ).

# Examples
```julia
# Calculate pair-breaking for parallel field
Λ = pairbreaking(1.5, 2, 0.25, 70, 70, 10)

# With angled field
Λ = pairbreaking(1.0, 1, 0.25, 70, 70, 10; θ=π/4)
```
"""
function pairbreaking(Φ, n, Δ0, ξd, R, d; θ = 0)
    Φpar = Φ * cos(θ)
    Φper = Φ * sin(θ)
    RLP = R + d/2
    Λpar =  (4 * (Φpar - n)^2 + d^2 / RLP^2 * (Φpar^2 + (n^2)/3))
    Λper = 4 * Φper^2
    Λ = ξd^2 * Δ0 / (1.76 * π * RLP^2) * (Λpar +  Λper)
    return maximum([real(Λ), real(Δ0) * 1e-3])
end

"""
    ΔD(Λ, Δ0)

Calculate the suppressed superconducting pairing in the presence of pair-breaking.

This function computes the reduced pairing Δ_d in the diffusive superconductor due to
pair-breaking effects characterized by Λ, using the Abrikosov solution.

# Arguments
- `Λ::Real`: Pair-breaking parameter (meV)
- `Δ0::Real`: Bulk superconducting pairing (meV)

# Returns
- `Δ_d::Real`: Suppressed pairing (meV)

# Notes
- Uses the `ΔΛ` function from AbrikosovSolver.jl
- Returns a safe value if Δ_d ≤ 0 to avoid numerical issues
"""
function ΔD(Λ, Δ0)
    Δd = ΔΛ(real(Λ), real(Δ0))
    if !(Δd > 0)
        # Δd = ΔΛ(real(Δ0/2 - imag(ω)), real(Δ0))
        Δd = ΔΛ(real(Δ0/2 - 1e-5), real(Δ0))
    end
    return Δd
end

"""
    Ω(Λ, Δ0)

Calculate the effective superconducting gap Ω.

This function computes the superconducting gap including pair-breaking.
# Arguments
- `Λ::Real`: Pair-breaking parameter (meV)
- `Δ0::Real`: Bulk superconducting gap (meV)

# Returns
- `Ω::Real`: Superconducting gap (meV)

# Formula
```
Ω = [(Δ_d^{2/3} - Λ^{2/3})^{3/2}]
```
where Δ_d is the suppressed gap from `ΔD(Λ, Δ0)`.
"""
function Ω(Λ, Δ0)
    return real(complex((complex(ΔD(Λ, Δ0))^(2/3) -  complex(Λ)^(2/3)))^(3/2))
end

# function usimple(Δ0, Λ, ω)
#     return (ω/Ω(Λ, Δ0, ω)) * sqrt(complex(1 - (ω/Ω)^2))
# end

# Modified ps_sign so uUsadel is analitical in the full complex plane. 
# WARNING: it is defined only for NEGATIVE frequencies.
function ps_sign(x)
    return real(x) < 0 ? -1.0 + imag(x)*im*1.0 : 1.0 + imag(x)*im*1.0
    #return sign(real(x)) + iω*1im
end

function uUsadel_old(Δ0, Λ, ω)
    #ω = ifelse(real(ω) == 0, imag(ω) + imag(ω)*im, ω)
    #Δd = Δ0 * (1 - π/4 * Λ/Δ0 - π^2/32 * (Λ/Δ0)^2 - π^3/96 * (Λ/Δ0)^3)
    Δd = ΔD(Λ, Δ0)
    pep = complex(-Δd^6 + 3 * Δd^4 * (Λ^2 + ω^2) + (Λ^2 + ω^2)^3 - 3 * Δd^2 * (Λ^4 - 16 * Λ^2 * ω^2 + ω^4) + 6 * Δd * Λ * sqrt(complex(-3*(Δd^2 - Λ^2)^3 * ω^2 + 9 * (Δd^4 + 7 * Δd^2 * Λ^2 + Λ^4) * ω^4 + 9 * (-Δd^2 + Λ^2) * ω^6 + 3 * ω^8)))^(1/3)
    nun = complex(ω^2 - Δd^2 + Λ^2)
    rai = sqrt(complex(ω^2 - 2 * nun / 3 + nun^2 / (3 * pep) + pep / 3))
    usaIn = (ω + ps_sign(ω) * rai - ps_sign(ω) * sqrt(complex(2 * ω^2 - 4 * nun/3 - nun^2 /(3 * pep) - pep/ 3 - ps_sign(ω) * 2 * (Δd^2 + Λ^2) * ω / rai))) /(2 * Δd) 
    return  usaIn
end

function uUsadel(Δ0, Λ, ω)
    Δd = ΔD(Λ, Δ0)
    pep = complex(-Δd^6 + 3 * Δd^4 * (Λ^2 + ω^2) + (Λ^2 + ω^2)^3 - 3 * Δd^2 * (Λ^4 - 16 * Λ^2 * ω^2 + ω^4) + 6 * Δd * Λ * sqrt(complex(-3*(Δd^2 - Λ^2)^3 * ω^2 + 9 * (Δd^4 + 7 * Δd^2 * Λ^2 + Λ^4) * ω^4 + 9 * (-Δd^2 + Λ^2) * ω^6 + 3 * ω^8)))^(1/3)
    nun = complex(ω^2 - Δd^2 + Λ^2)
    rai = sqrt(complex(ω^2 - 2 * nun / 3 + nun^2 / (3 * pep) + pep / 3))
    usa = complex(ω / (2Δd) - rai / (2Δd) + sqrt(complex(2ω^2 - 4 * nun /3 - nun^2 / (3 * pep) - pep /3 + 2 * (Δd^2 + Λ^2) * ω / rai)) / (2Δd))
    # Guarantee analyticality in the full complex plane
    usa = real(usa) + abs(imag(usa)) * 1im
    return usa
end

"""
    itip(Δ0, Λ)

Calculate the threshold imaginary frequency for the Usadel solution.

This function computes the imaginary frequency iω_tip above which the Usadel solution
becomes unphysical (complex or multivalued). It's crucial for determining valid
integration paths in the complex frequency plane.

# Arguments
- `Δ0::Number`: Bulk superconducting pairing (meV)
- `Λ::Real`: Pair-breaking parameter (meV)

# Returns
- `iω_tip::Real`: Threshold imaginary frequency (meV)

# Notes
- Returns a small positive value (1e-5) if the calculation yields a complex result
- This indicates the integration path is not valid for the given parameters
- Logs an error message if the result is complex

# Examples
```julia
Λ = pairbreaking(1.0, 1, 0.25, 70, 70, 10)
ω_threshold = itip(0.25, Λ)
```
"""
function itip(Δ0, Λ)
    Δd = ΔD(Λ, Δ0)
    x = Λ/Δd
    p1 = -4 + 4*x^2
    n2 = 2 * 2^(1/3) * (1 - 47 * x^2 + x^4)
    d2 = try
        (2 + 345 * x^2 - 345 * x^4 - 2 * x^6 + 9 * sqrt(3) * x * (1 + x^2) * sqrt(8 + 359 * x^2 + 8 * x^4))^(1/3) 
    catch
        @error "Itip is complex. Integration path not valid"
    end
    p2 = n2 / d2
    p3 = 2^(2/3) * d2
    iωm =  sqrt(complex((p1 + p2 + p3) / 6)) * Δd
    return ifelse(isreal(iωm), iωm, 1e-5)
end

"""
    LP_lobe(n, ξd, R, d)

Calculate the Little-Parks lobe boundaries for a given winding number.

The Little-Parks effect causes oscillations in the superconducting pairing
as a function of magnetic flux. This function computes the flux range (lobe) where
superconductivity survives for a given winding number n.

# Arguments
- `n::Int`: Winding number
- `ξd::Real`: Superconducting coherence length (nm)
- `R::Real`: Inner radius of the cylinder (nm)
- `d::Real`: Width of the superconductor shell (nm)

# Returns
- `ΦLPa::Real`: Lower boundary of the Little-Parks lobe (in units of Φ₀)
- `ΦLPb::Real`: Upper boundary of the Little-Parks lobe (in units of Φ₀)

# Notes
The lobe is centered around the winding number n and has a width that depends on
the coherence length and geometry.

# Examples
```julia
# Get lobe boundaries for n=1
Φ_min, Φ_max = LP_lobe(1, 70, 70, 10)

# Check if flux is in the lobe
Φ = 1.2
in_lobe = is_in_lobe(Φ, Φ_min, Φ_max)
```
"""
function LP_lobe(n, ξd, R, d)
    RLP = R + d/2
    pre = 1.76 * π * exp(-π/4) / 4
    frac1 = 1 / (1 + (d / (2 * RLP))^2)
    frac2 = d^2 / (d^2 + 4 * RLP^2)
    frac3 = RLP^2 / ξd^2
    root = pre * frac3 * frac1 + n^2 * frac2 * (frac2 - 4/3)
    ΦLPa = n * frac1 - sqrt(root)
    ΦLPb = n * frac1 + sqrt(root)
    return ΦLPa, ΦLPb
end

function isdestructive(n, ξd, R, d)
    ΦLPa, ΦLPb = ΦLP(n, ξd, R, d)
    return !(ΦLPb - n <= 1/2)
end

"""
    is_in_lobe(Φ, ΦLPa, ΦLPb)

Checks whether the value `Φ` lies within the interval defined by `ΦLPa` and `ΦLPb` (inclusive).

# Arguments
- `Φ`: The value to check.
- `ΦLPa`: The lower bound of the interval.
- `ΦLPb`: The upper bound of the interval.

# Returns
- `Bool`: `true` if `Φ` is within `[ΦLPa, ΦLPb]`, otherwise `false`.
"""
function is_in_lobe(Φ, ΦLPa, ΦLPb)
    return (ΦLPa <= Φ) && (Φ <= ΦLPb)
end

