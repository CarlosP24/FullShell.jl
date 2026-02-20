"""
    FullShell

A Julia package for building full-shell hybrid semiconductor-superconductor Hamiltonians.

FullShell provides tools to construct and analyze tight-binding Hamiltonians for cylindrical
nanowire heterostructures consisting of a semiconductor core covered by a superconducting shell.
The package supports both standard and multi-mJ models, with various physical effects including:

- Spin-orbit coupling (Rashba and Dresselhaus)
- Zeeman coupling
- Magnetic field effects (orbital and Zeeman)
- Proximity-induced superconductivity (Usadel and ballistic regimes)
- Electrostatic dome profiles

The package is built on top of [Quantica.jl](https://github.com/pablosanjose/Quantica.jl).

# Exports
- `build_cyl`, `build_cyl_mm`: Main Hamiltonian constructors
- `Params`, `Params_mm`: Parameter structures
- `pairbreaking`, `LP_lobe`, `ΔΛ`: Superconducting proximity functions
- `get_itip`, `bandwidth`, `get_Φ`, `get_B`, `get_Ω`: Helper functions
- `mux`, `mu_step`: Hamiltonian modifiers
- Pauli matrices and operators: `σ0τx`, `σ0τy`, `σ0τz`, etc.

# Examples
```julia
using FullShell

# Define model parameters (default units: meV and nm)
model = (
    R = 70,      # Cylinder radius
    w = 10,      # Semiconductor width
    d = 10,      # Superconductor width
    Δ0 = 0.23,   # Superconducting gap
    μ = 0.75,    # Chemical potential
    α = 0.2,     # Rashba spin-orbit coupling
)

# Build Hamiltonians
hSM, hSC, parameters = build_cyl(; model...)
```
"""
module FullShell

using Quantica 
using Parameters
using Roots
using StaticArrays

export build_cyl, build_cyl_mm, ΣS3DUsadel, pairbreaking, LP_lobe, ΔΛ, get_itip, bandwidth, get_Φ, get_B, get_Ω,build_harmonic_deformations

export ΔD, is_in_lobe

export mux, mu_step

export Params, Params_mm

export σ0τx, σ0τy, σ0τz, σ0τ0, σzτ0, σzτz, σyτy, σyτz, σyτ0, σxτz, σxτ0, σ0, σx, σy, σz, c_up, c_down, Ω

include("utilities/LP_usadel.jl")
include("utilities/pauli_products.jl") 
include("HamiltonianBuilder.jl")
include("HamiltonianBuilder_MM.jl")
include("modifiers/mux.jl")

end