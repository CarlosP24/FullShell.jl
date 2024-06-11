using FullShell

base_model = (;
    R = 70,
    w = 0,
    d = 0,
    Δ0 = 0.23,
    ξd = 70,
    a0 = 5,
    preα = 0,
    g = 0,
    α = 0,
    μ = 0,
    τΓ = 1,
    Vexponent = 2,
    Vmin = 0,
    Vmax = 0,
    ishollow = true,
)

model = (;  base_model...,
    R = 30,
    w = 30,
    d = 10,
    Vmin = -30,
    g = 10,
    τΓ = 40,
    μ = 2,
    preα = 46.66,
    ishollow = false,
)

hSM, hSC, params = build_cyl(; model..., )

function ps_sign(x)
    return real(x) < 0 ? -1.0 + imag(x)*im*1.0 : 1.0 + imag(x)*im*1.0
end

function uUsadel(Δ0, Λ, ω)
    ω = ifelse(real(ω) == 0, imag(ω) + imag(ω)*im, ω)
    #Δd = Δ0 * (1 - π/4 * Λ/Δ0 - π^2/32 * (Λ/Δ0)^2 - π^3/96 * (Λ/Δ0)^3)
    Δd = ΔΛ(real(Λ), real(Δ0))
    if !(Δd > 0)
        Δd = ΔΛ(real(Δ0/2 - imag(ω)), real(Δ0))
    end
    pep = complex(-Δd^6 + 3 * Δd^4 * (Λ^2 + ω^2) + (Λ^2 + ω^2)^3 - 3 * Δd^2 * (Λ^4 - 16 * Λ^2 * ω^2 + ω^4) + 
    6 * Δd * Λ * sqrt(complex(-3*(Δd^2 - Λ^2)^3 * ω^2 + 9 * (Δd^4 + 7 * Δd^2 * Λ^2 + Λ^4) * ω^4 + 9 * (-Δd^2 + Λ^2) * ω^6 + 3 * ω^8)))^(1/3)
    nun = ω^2 - Δd^2 + Λ^2
    rai = sqrt(complex(ω^2 - 2 * nun / 3 + nun^2 / (3 * pep) + pep / 3))
    usaIn = 1/(2 * Δd) *(ω + ps_sign(ω) * rai - ps_sign(ω) * sqrt(complex(2 * ω^2 - 4 * nun/3 - nun^2 /(3 * pep) - pep/ 3 - ps_sign(ω) * 2 * (Δd^2 + Λ^2) * ω / rai)))
    return  usaIn
end

uOut(ω, Λ) =  ps_sign(ω) * ω / sqrt(complex(ω^2 + Λ^2))
function uIn(Δ0, Λ, ω)
    ω = ifelse(real(ω) == 0, imag(ω) + imag(ω)*im, ω)
    Δd = ΔΛ(real(Λ), real(Δ0))
    pep = complex(-Δd^6 + 3 * Δd^4 * (Λ^2 + ω^2) + (Λ^2 + ω^2)^3 - 3 * Δd^2 * (Λ^4 - 16 * Λ^2 * ω^2 + ω^4) + 
    6 * Δd * Λ * sqrt(complex(-3*(Δd^2 - Λ^2)^3 * ω^2 + 9 * (Δd^4 + 7 * Δd^2 * Λ^2 + Λ^4) * ω^4 + 9 * (-Δd^2 + Λ^2) * ω^6 + 3 * ω^8)))^(1/3)
    nun = ω^2 - Δd^2 + Λ^2
    rai = sqrt(complex(ω^2 - 2 * nun / 3 + nun^2 / (3 * pep) + pep / 3))
    usaIn = 1/(2 * Δd) *(ω + ps_sign(ω) * rai - ps_sign(ω) * sqrt(complex(2 * ω^2 - 4 * nun/3 - nun^2 /(3 * pep) - pep/ 3 - ps_sign(ω) * 2 * (Δd^2 + Λ^2) * ω / rai)))
    return usaIn
end

ΣS3DUsadel(Δ0, Λ, ω;) = -(uUsadel(Δ0, Λ, ω) * σ0τ0 - σ0τx) / sqrt(complex(1-uUsadel(Δ0, Λ, ω)^2))


using CairoMakie

Λ(Φ) = FullShell.pairbreaking(Φ, round(Int, Φ), params.Δ0, params.ξd, params.R, params.d)
h(Φ) = real(Λ(Φ)) - real(params.Δ0/2)
zs = find_zeros(h, 0, 3.5)

Φrng = range(0, 3.5, length = 50)
ωrng = range(-0.26, 0.26, length = 51) .+ 1e-4im;
Σs = map(Iterators.product(Φrng, ωrng)) do (Φ, ω)
    return ΣS3DUsadel(params.Δ0, Λ(Φ), ω)
end

us = map(Iterators.product(Φrng, ωrng)) do (Φ, ω)
    return uOut(ω, Λ(Φ))
end

uIns = map(Iterators.product(Φrng, ωrng)) do (Φ, ω)
    return uIn(params.Δ0, Λ(Φ), ω)
end

fig = Figure()
ax = Axis(fig[1, 1]; xlabel = L"\Phi / \Phi_0", ylabel = L"\omega")
heatmap!(ax, Φrng, real.(ωrng), imag.(tr.(inv.(Σs))); colormap = :thermal)
vlines!(ax, zs; color = :white, linestyle = :dash)
fig

##
fig = Figure()
ax = Axis(fig[1, 1]; xlabel = L"\Phi / \Phi_0", ylabel = L"\omega")
hmap = heatmap!(ax, Φrng, real.(ωrng), imag.(us); colormap = :thermal, )
Colorbar(fig[1, 2], hmap)
ax = Axis(fig[2, 1]; xlabel = L"\Phi / \Phi_0", ylabel = L"\omega")
hmap = heatmap!(ax, Φrng, real.(ωrng), imag.(uIns); colormap = :thermal)
vlines!(ax, zs; color = :white, linestyle = :dash)
Colorbar(fig[2, 2], hmap)
fig