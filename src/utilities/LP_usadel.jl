include("AbrikosovSolver.jl")

function pairbreaking(Φ, n, Δ0, ξd, R, d)
    RLP = R + d/2
    Λ = ξd^2 * Δ0 / (1.76 * π * RLP^2) * (4 * (Φ - n)^2 + d^2 / RLP^2 * (Φ^2 + (n^2)/3))
    return Λ
end

# function ps_sign(x)
#     return real(x) < 0 ? -1.0 + imag(x)*im*1.0 : 1.0 + imag(x)*im*1.0
#     #return -1.0 + 1e-5im
# end

function ΔD(Λ, Δ0, ω)
    Δd = ΔΛ(real(Λ), real(Δ0))
    if !(Δd > 0)
        # Δd = ΔΛ(real(Δ0/2 - imag(ω)), real(Δ0))
        Δd = ΔΛ(real(Δ0/2 - 1e-5), real(Δ0))
    end
    return Δd
end

function Ω(Λ, Δ0, ω)
    return (ΔD(Λ, Δ0, ω)^(2/3) -  Λ^(2/3))^(3/2)
end

# function usimple(Δ0, Λ, ω)
#     return (ω/Ω(Λ, Δ0, ω)) * sqrt(complex(1 - (ω/Ω)^2))
# end

# function uUsadel(Δ0, Λ, ω)
#     ω = ifelse(real(ω) == 0, imag(ω) + imag(ω)*im, ω)
#     #Δd = Δ0 * (1 - π/4 * Λ/Δ0 - π^2/32 * (Λ/Δ0)^2 - π^3/96 * (Λ/Δ0)^3)
#     Δd = ΔD(Λ, Δ0, ω)
#     pep = complex(-Δd^6 + 3 * Δd^4 * (Λ^2 + ω^2) + (Λ^2 + ω^2)^3 - 3 * Δd^2 * (Λ^4 - 16 * Λ^2 * ω^2 + ω^4) + 
#     6 * Δd * Λ * sqrt(complex(-3*(Δd^2 - Λ^2)^3 * ω^2 + 9 * (Δd^4 + 7 * Δd^2 * Λ^2 + Λ^4) * ω^4 + 9 * (-Δd^2 + Λ^2) * ω^6 + 3 * ω^8)))^(1/3)
#     nun = ω^2 - Δd^2 + Λ^2
#     rai = sqrt(complex(ω^2 - 2 * nun / 3 + nun^2 / (3 * pep) + pep / 3))
#     usaIn = 1/(2 * Δd) *(ω + ps_sign(ω) * rai - ps_sign(ω) * sqrt(complex(2 * ω^2 - 4 * nun/3 - nun^2 /(3 * pep) - pep/ 3 - ps_sign(ω) * 2 * (Δd^2 + Λ^2) * ω / rai)))
#     return  usaIn
# end

# Removed ps_sign so uUsadel is analitical in the full complex plane. 
# WARNING: it is defined only for NEGATIVE frequencies.
function uUsadel(Δ0, Λ, ω)
    #ω = ifelse(real(ω) == 0, imag(ω) + imag(ω)*im, ω)
    #Δd = Δ0 * (1 - π/4 * Λ/Δ0 - π^2/32 * (Λ/Δ0)^2 - π^3/96 * (Λ/Δ0)^3)
    Δd = ΔD(Λ, Δ0, ω)
    pep = complex(-Δd^6 + 3 * Δd^4 * (Λ^2 + ω^2) + (Λ^2 + ω^2)^3 - 3 * Δd^2 * (Λ^4 - 16 * Λ^2 * ω^2 + ω^4) + 
    6 * Δd * Λ * sqrt(complex(-3*(Δd^2 - Λ^2)^3 * ω^2 + 9 * (Δd^4 + 7 * Δd^2 * Λ^2 + Λ^4) * ω^4 + 9 * (-Δd^2 + Λ^2) * ω^6 + 3 * ω^8)))^(1/3)
    nun = ω^2 - Δd^2 + Λ^2
    rai = sqrt(complex(ω^2 - 2 * nun / 3 + nun^2 / (3 * pep) + pep / 3))
    usaIn = 1/(2 * Δd) *(ω - rai + sqrt(complex(2 * ω^2 - 4 * nun/3 - nun^2 /(3 * pep) - pep/ 3 + 2 * (Δd^2 + Λ^2) * ω / rai)))
    return  usaIn
end

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

function is_in_lobe(Φ, ΦLPa, ΦLPb)
    return (ΦLPa <= Φ) && (Φ <= ΦLPb)
end

