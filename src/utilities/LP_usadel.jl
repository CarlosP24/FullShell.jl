function pairbreaking(Φ, n, Δ0, ξd, R, d)
    RLP = R + d/2
    Λ = ξd^2 * Δ0 / (1.76 * π * RLP^2) * (4 * (Φ - n)^2 + d^2 / RLP^2 * (Φ^2 + (n^2)/3))
    return Λ
end

function uUsadel(Δ0, Λ, ω)
    Δd = Δ0 * (1 - π/4 * Λ/Δ0 - π^2/32 * (Λ/Δ0)^2 - π^3/96 * (Λ/Δ0)^3)
    pep = complex(-Δd^6 + 3 * Δd^4 * (Λ^2 + ω^2) + (Λ^2 + ω^2)^3 - 3 * Δd^2 * (Λ^4 - 16 * Λ^2 * ω^2 + ω^4) + 
    6 * Δd * Λ * sqrt(complex(-3*(Δd^2 - Λ^2)^3 * ω^2 + 9 * (Δd^4 + 7 * Δd^2 * Λ^2 + Λ^4) * ω^4 + 9 * (-Δd^2 + Λ^2) * ω^6 + 3 * ω^8)))^(1/3)
    nun = ω^2 - Δd^2 + Λ^2
    rai = sqrt(complex(ω^2 - 2 * nun / 3 + nun^2 / (3 * pep) + pep / 3))
    usa = 1/(2 * Δd) *(ω + sign(ω) * rai - sign(ω) * sqrt(complex(2 * ω^2 - 4 * nun/3 - nun^2 /(3 * pep) - pep/ 3 - sign(ω) * 2 * (Δd^2 + Λ^2) * ω / rai)))
    return usa
end
