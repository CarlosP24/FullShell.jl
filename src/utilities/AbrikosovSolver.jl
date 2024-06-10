function P(z)
    if z <= 1
        p = π/4 * z 
    else
        p = log(z + sqrt(z^2 - 1)) + z/2 * atan(1/sqrt(z^2 - 1)) - sqrt(z^2 - 1)/(2*z)
    end
    return p
end

f(Δd, Δ0 = 0.23; Λ = 0) = log(Δd/Δ0) + P(Λ/Δd)

function ΔΛ(Λ, Δ0)
    g(Δd) = f(Δd, Δ0; Λ)
    return get(find_zeros(g, 0, Δ0), 1, 0)
end
