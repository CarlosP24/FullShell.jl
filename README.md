# FullShell.jl

[![DOI](https://zenodo.org/badge/794025449.svg)](https://doi.org/10.5281/zenodo.11450677)
[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://CarlosP24.github.io/FullShell.jl/dev/)

A Julia package for building full-shell hybrid semiconductor-superconductor Hamiltonians.

Uses [Quantica.jl](https://github.com/pablosanjose/Quantica.jl). 

## Documentation

📚 **[Read the full documentation](https://CarlosP24.github.io/FullShell.jl/dev/)**

The documentation includes:
- Comprehensive API reference with detailed docstrings
- Extensive usage examples
- Physical model descriptions
- Parameter guides

## Installation

### From GitHub (Current)
This package is not yet registered in Julia. To install it, run:
```julia
using Pkg
Pkg.add("https://github.com/CarlosP24/FullShell.jl.git")
```

### After Registration
Once registered, you will be able to install it simply with:
```julia
using Pkg
Pkg.add("FullShell")
```

## Quick Start

Default parameters are in meV and nm.
````
julia> using FullShell

julia> model = (;
        R = 70,
        w = 0,
        d = 0,
        preα = 0,
        Δ0 = 0.23,
        ξd = 70,
        τΓ = 1,
        g = 0,
        a0 = 5,
        μ = 0.75,
        α = 0,
        )
(R = 70, w = 0, d = 0, preα = 0, Δ0 = 0.23, ξd = 70, τΓ = 1, g = 0, a0 = 5, μ = 0.75, α = 0)

julia> hSM, hSC, parameters = build_cyl(; model...)
(ParametricHamiltonian{Float64,2,1}: Parametric Hamiltonian on a 1D Lattice in 2D space
  Bloch harmonics  : 3
  Harmonic size    : 1 × 1
  Orbitals         : [4]
  Element type     : 4 × 4 blocks (ComplexF64)
  Onsites          : 1
  Hoppings         : 2
  Coordination     : 2.0
  Parameters       : [:Vmax, :Vmin, :Z, :preα, :Φ, :α, :μ], ParametricHamiltonian{Float64,2,1}: Parametric Hamiltonian on a 1D Lattice in 2D space
  Bloch harmonics  : 3
  Harmonic size    : 1 × 1
  Orbitals         : [4]
  Element type     : 4 × 4 blocks (ComplexF64)
  Onsites          : 1
  Hoppings         : 2
  Coordination     : 2.0
  Parameters       : [:Vmax, :Vmin, :Z, :preα, :Φ, :α, :μ, :τΓ, :ω], FullShell.Params
  ħ2ome: Float64 76.1996
  μBΦ0: Float64 119.6941183
  m0: Float64 0.023
  g: Float64 0.0
  P: Float64 919.7
  Δg: Float64 417.0
  Δs: Float64 390.0
  preα: Float64 0.0
  a0: Float64 5.0
  t: Float64 66.26052173913044
  echarge: Float64 1.0
  R: Float64 70.0
  w: Float64 0.0
  d: Float64 0.0
  Vmax: Float64 0.0
  Vmin: Float64 0.0
  Vexponent: Float64 2.0
  Δ0: ComplexF64
  ξd: Float64 70.0
  α: Float64 0.0
  μ: Float64 0.75
  τΓ: Float64 1.0
  Φ: Float64 1.0
)
````

### Calculating the LDOS at the edge of the nanowire
First, we need to obtain the Greens function at the edge of the nanowire. 
As we work with a semi-infinite system, we need to use a Schur `Quantica.GreenSolver`(see [Quantica documentation](https://pablosanjose.github.io/Quantica.jl/dev/tutorial/greenfunctions/)).
```
julia> using Quantica
julia> g = hSC |> greenfunction(GS.Schur(boundary = 0))
GreenFunction{Float64,2,1}: Green function of a Hamiltonian{Float64,2,1}
  Solver          : AppliedSchurGreenSolver
  Contacts        : 0
  Contact solvers : ()
  Contact sizes   : ()
  ParametricHamiltonian{Float64,2,1}: Parametric Hamiltonian on a 1D Lattice in 2D space
    Bloch harmonics  : 3
    Harmonic size    : 1 × 1
    Orbitals         : [4]
    Element type     : 4 × 4 blocks (ComplexF64)
    Onsites          : 1
    Hoppings         : 2
    Coordination     : 2.0  
    Parameters       : [:Vmax, :Vmin, :Z, :preα, :Φ, :α, :μ, :τΓ, :ω]

```
Second, we build the LDOS function at the edge of the nanowire (cells `-1`or `1`as the Shur boundary is at cell `0`), which depends on the energy and all parameters of our Greens function (and Hamiltonian):
````
julia> ρ = ldos(g[cells = (-1)])
LocalSpectralDensitySlice{Float64} : local density of states at a fixed location and arbitrary energy
  kernel   : LinearAlgebra.UniformScaling{Bool}(true)
````

We are now able to calculate the LDOS at the edge of the semi-infinite for whatever set of parameters that we want, having set all defaults when defining the `model`object. If not specified there, defaults are those defined by `FullShell.jl`:
````
julia> ω = 0.0 + 1e-4im
0.0 + 0.0001im
julia> ρ(ω; ω = ω, Φ = 1, Z = 0)
1-element OrbitalSliceVector{Vector{Float64}}:
 1.6677498905272993e-6
````

For more examples and detailed documentation, see the [full documentation](https://CarlosP24.github.io/FullShell.jl/dev/).

## Features

- **Full-shell geometry**: Model semiconductor nanowires with complete superconducting coverage
- **Multiple physical effects**:
  - Rashba and Dresselhaus spin-orbit coupling
  - Zeeman coupling to magnetic fields
  - Orbital magnetic field effects
  - Proximity-induced superconductivity (Usadel and ballistic regimes)
  - Electrostatic dome profiles
- **Two modeling approaches**:
  - Standard approach with continuous angular coordinate
  - Multi-mJ approach with discrete angular momentum channels
- **Comprehensive documentation** with examples and API reference

## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines on:
- Setting up the development environment
- Building documentation locally
- Code style and testing
- Submitting pull requests

## Citation

If you use FullShell.jl in your research, please cite:

```bibtex
@software{fullshell_jl,
  author = {Payá, Carlos},
  title = {FullShell.jl: Full-shell hybrid semiconductor-superconductor Hamiltonians},
  year = {2024},
  doi = {10.5281/zenodo.11450677},
  url = {https://github.com/CarlosP24/FullShell.jl}
}
```

## License

See the LICENSE file for details.
