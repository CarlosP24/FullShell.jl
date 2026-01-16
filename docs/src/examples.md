# Examples

This page provides practical examples of using FullShell.jl for different scenarios.

## Basic Example: Building a Hamiltonian

The simplest way to build a full-shell Hamiltonian is to specify parameters directly:

```julia
using FullShell

# Build with default parameters
hSM, hSC, params = build_cyl()

# Build with custom parameters
hSM, hSC, params = build_cyl(
    R = 70,      # Cylinder radius (nm)
    w = 10,      # Semiconductor width (nm)
    d = 10,      # Superconductor width (nm)
    Δ0 = 0.23,   # Superconducting gap (meV)
    μ = 0.75,    # Chemical potential (meV)
    α = 0.2,     # Rashba SOC (meV·nm)
    Φ = 1.0,     # Magnetic flux (Φ/Φ₀)
)
```

## Using Parameter Structures

For more complex workflows, it's convenient to use the `Params` structure:

```julia
using FullShell

# Create parameter structure
p = Params(
    R = 100,
    w = 15,
    d = 10,
    Δ0 = 0.25,
    ξd = 70,
    μ = 0.5,
    α = 0.3,
    g = 12,
    Φ = 0.5,
)

# Build Hamiltonians
hSM, hSC, params = build_cyl(p)

# Access parameters
println("Radius: ", params.R, " nm")
println("Chemical potential: ", params.μ, " meV")
```

## Parametric Hamiltonians

The returned Hamiltonians are parametric and can be evaluated at different parameter values:

```julia
using FullShell, Quantica

hSM, hSC, params = build_cyl(R=70, μ=0.5, α=0.2, Φ=1.0)

# Evaluate at specific parameter values
h_μ1 = hSC(μ=0.5, Φ=1.0, α=0.2, ω=0)
h_μ2 = hSC(μ=0.8, Φ=1.0, α=0.2, ω=0)

# Compute band structure
using Plots
bands_μ1 = bandstructure(h_μ1, range(-π, π, 100))
bands_μ2 = bandstructure(h_μ2, range(-π, π, 100))

plot(bands_μ1, label="μ = 0.5 meV")
plot!(bands_μ2, label="μ = 0.8 meV")
```

## Scanning Parameters

Compute properties as a function of parameters:

```julia
using FullShell, Quantica

hSM, hSC, params = build_cyl(R=70, α=0.2)

# Scan chemical potential
μ_values = range(0, 2, length=50)
gap_closing = zeros(length(μ_values))

for (i, μ) in enumerate(μ_values)
    h = hSC(μ=μ, Φ=1.0, α=0.2, ω=0)
    bands = bandstructure(h, range(-π, π, 100))
    # Extract minimum gap (implementation depends on Quantica.jl API)
    gap_closing[i] = minimum(abs.(energies(bands)))
end

using Plots
plot(μ_values, gap_closing, xlabel="μ (meV)", ylabel="Min gap (meV)")
```

## Different Superconductor Models

You can choose between different models for the superconducting shell:

```julia
using FullShell

# Usadel (diffusive) superconductor (default)
hSM, hSC_usadel, _ = build_cyl(R=70, Δ0=0.25, shell="Usadel")

# Ballistic superconductor
hSM, hSC_ballistic, _ = build_cyl(R=70, Δ0=0.25, shell="Ballistic")

# Simplified gap suppression model
hSM, hSC_simple, _ = build_cyl(R=70, Δ0=0.25, shell="Simple")
```

## Electrostatic Dome Profile

Add a dome-shaped electrostatic potential:

```julia
using FullShell

# Create a dome profile with V(ρ) = Vmax + (Vmin - Vmax)(ρ/R)^n
hSM, hSC, params = build_cyl(
    R = 70,
    Vmax = 1.0,      # Maximum potential at inner radius (meV)
    Vmin = 0.0,      # Minimum potential at outer radius (meV)
    Vexponent = 2,   # Exponent controlling profile shape
)
```

## Multi-mJ Model

For explicit angular momentum channels:

```julia
using FullShell

# Build multi-mJ Hamiltonian
hSM, hSC, params = build_cyl_mm(
    R = 70,
    num_mJ = 5,      # Number of angular momentum channels
    B = 0.1,         # Magnetic field (Tesla)
    α = 0.2,
    μ = 0.5,
)

# Enable hopping between mJ channels
hSM, hSC, params = build_cyl_mm(
    R = 70,
    num_mJ = 5,
    hops0 = true,        # Enable mJ hopping
    range_hop_m = 15,    # Hopping range in units of a0
    B = 0.1,
)
```

## Chemical Potential Modulation

Apply position-dependent chemical potential:

```julia
using FullShell

# Build base Hamiltonian
hSM, hSC, params = build_cyl(R=70, α=0.2, L=500)

# Define step function
μ_profile(x) = 0.5 * (1 + tanh((x - 250) / 50))

# Apply modulation
h_modulated = mux(hSC, μ_profile, params.L)

# Or use the built-in step function
h_step = mu_step(hSC, L=500, Lstep=250, ς=50, μshift=0.5)
```

## Flux Conversion

Convert between magnetic field and flux:

```julia
using FullShell

params = Params(R=70, d=10)

# Get conversion functions
Φ_from_B = get_Φ(params)
B_from_Φ = get_B(params)

# Convert
B = 0.1  # Tesla
Φ = Φ_from_B(B)  # in units of Φ₀
println("B = $B T corresponds to Φ = $Φ Φ₀")

# Convert back
B_check = B_from_Φ(Φ)
println("Φ = $Φ Φ₀ corresponds to B = $B_check T")
```

## Little-Parks Oscillations

Calculate Little-Parks lobe boundaries:

```julia
using FullShell

params = Params(R=70, d=10, ξd=70)

# Get lobe boundaries for winding number n=1
Φ_min, Φ_max = LP_lobe(1, params.ξd, params.R, params.d)

println("Little-Parks lobe for n=1: [$Φ_min, $Φ_max] Φ₀")

# Check if a flux is in the lobe
Φ = 1.2
in_lobe = is_in_lobe(Φ, Φ_min, Φ_max)
println("Φ = $Φ Φ₀ is ", in_lobe ? "inside" : "outside", " the lobe")
```

## Effective Gap and Integration Path

Get the effective superconducting gap and valid imaginary frequency:

```julia
using FullShell

params = Params(R=70, d=10, Δ0=0.25, ξd=70)

# Get effective gap as function of flux
Ω_func = get_Ω(params)
Φ = 1.0
Ω_eff = Ω_func(Φ)
println("Effective gap at Φ = $Φ Φ₀: Ω = $Ω_eff meV")

# Get imaginary frequency threshold
itip_func = get_itip(params)
iω_max = itip_func(Φ)
println("Maximum safe imaginary frequency: iω = $iω_max meV")
```

## Angled Magnetic Field

Apply a magnetic field at an angle to the wire axis:

```julia
using FullShell

# θ = 0: field parallel to wire axis
hSM, hSC, _ = build_cyl(R=70, Φ=1.0, θ=0.0)

# θ = π/4: field at 45 degrees
hSM, hSC, _ = build_cyl(R=70, Φ=1.0, θ=π/4)

# θ = π/2: field perpendicular to wire axis
hSM, hSC, _ = build_cyl(R=70, Φ=1.0, θ=π/2)
```

## Bandwidth Calculation

Calculate the bandwidth of the semiconductor:

```julia
using FullShell

params = Params(R=70, μ=0.5, m0=0.023, a0=5)
bw = bandwidth(params)
println("Bandwidth: $bw meV")
```
