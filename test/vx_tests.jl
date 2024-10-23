using FullShell
using Distributed
using CairoMakie

addprocs(32)

@everywhere begin
    using FullShell
    using Quantica
    using Parameters
    using ProgressMeter
end

##
function pldos(ρ, Φs, ωs, Zs; τ = 1, φ = 0)
    pts = Iterators.product(Φs, ωs, Zs)
    LDOS = @showprogress pmap(pts) do pt
        Φ, ω, Z = pt 
        ld = try 
            ρ(ω; ω = ω, Φ = Φ, Z = Z, τ = τ, phase = φ)
        catch
            0.0
        end
        return ld
    end
    LDOSarray = reshape(LDOS, size(pts)...)
    LDOS = sum.(LDOSarray)
    return Dict([ Zs[i] => LDOS[:, :, i] for i in 1:size(LDOS, 3)])
end

function pldos_x(g, ωs, Zs; Φ = 0.51,)
    pts = Iterators.product(ωs, Zs)
    LDOS = @showprogress pmap(pts) do pt
        ω, Z = pt 
        ld = try
            ldos(g(ω; ω, Φ, Z))[]
        catch
            0.0
        end
        return ld
    end
    LDOS = reshape(LDOS, size(pts)...)
    return Dict([Zs[i] => LDOS[:, i] for i in 1:size(LDOS, 2)])
end

function pfdo(ρ, xrng)
    fdo = @showprogress map(xrng) do x
        return ρ[cells = (x)]
    end
    return reshape(fdo, size(xrng)...)
end
##
wire = (;
    R = 70,
    w = 20,
    d = 10,
    μ = 0,
    α = 50,
    preα = 0,
    Vmax = -0.5,
    Δ0 = 0.23,
    ξd = 70,
    τΓ = 1.1,
    g = 10,
    ς = 1,                      
    Lstep = 10,                       
    Vshift = 0.5
)

hSM, hSC, params = build_cyl(; wire..., nforced = 1)
g = hSC |> greenfunction(GS.Schur(boundary = 0))


## Finite 
L = 1000
hf = hSC |> supercell(region = r -> 0 < r[1] < L)
gf = hf |> greenfunction()

##
Φrng = range(-5, 1.499, length = 100)
ωrng = range(-.26, 0, length = 101) .+ 1e-3im
Zs = 0

LDOS = pldos(ldos(gf[cells = (-1)]), Φrng, ωrng, Zs)

##
ω = 0.0 + 1e-3im
Φ = 0.51
Z = 0
xrng = 0:500
ρ = ldos(gf(ω; Φ, Z))
Ψ2 = pfdo(ρ)
##
ωrng = range(-.26, 0, length = 101) .+ 1e-3im
Zs = 0

xLDOS = pldos_x(gf, ωrng, Zs; Φ = 0.51)
##
sp = spectrum(hf(; Z = 0, Φ = 0.52, omegamap = ω -> (; ω)))

##
fig = Figure()
ax = Axis(fig[1, 1]; xlabel = L"\Phi/\Phi_0", ylabel = L"\omega")
fωrng = vcat(real.(ωrng), -real.(ωrng) |> reverse)
sLDOS = sum(values(LDOS))
fLDOS = abs.(hcat(sLDOS, reverse(sLDOS, dims = 2)))
heatmap!(ax, Φrng, fωrng, fLDOS; colormap = :thermal, colorrange = (5, 1e2))
fig

##
fig = Figure()
ax = Axis(fig[1, 1]; xlabel = "x", ylabel = "ω")
fωrng = vcat(real.(ωrng), -real.(ωrng) |> reverse)
sLDOS = mapreduce(permutedims, vcat, sum(values(xLDOS))) |> transpose
fLDOS = abs.(hcat(sLDOS, reverse(sLDOS, dims = 2)))
heatmap!(ax, (0:first(size(sLDOS))).*5, fωrng, fLDOS; colormap = :thermal, colorrange = (1e-4, 1e-1))
#hlines!(ax, [-0.0095]; color = :black)
fig