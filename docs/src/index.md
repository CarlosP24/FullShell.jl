# FullShell.jl

[![DOI](https://zenodo.org/badge/794025449.svg)](https://doi.org/10.5281/zenodo.11450677)

A Julia package for building full-shell hybrid semiconductor-superconductor Hamiltonians.

## Overview

FullShell.jl provides tools to construct and analyze tight-binding Hamiltonians for cylindrical
nanowire heterostructures consisting of a semiconductor core covered by a superconducting shell.
These structures are of great interest for quantum computing and topological physics applications.

## Features

- **Full-shell geometry**: Model semiconductor nanowires with complete superconducting coverage
- **Physical effects**:
  - Rashba and Dresselhaus spin-orbit coupling
  - Zeeman coupling to magnetic fields
  - Orbital magnetic field effects in cylindrical geometry
  - Proximity-induced superconductivity (Usadel and ballistic regimes)
  - Electrostatic dome profiles
- **Two modeling approaches**:
  - Standard approach with continuous angular coordinate
  - Multi-mJ approach with discrete angular momentum channels
- **Built on [Quantica.jl](https://github.com/pablosanjose/Quantica.jl)**: Leverages powerful tight-binding framework

## Installation

FullShell.jl is registered in the Julia General Registry. Install it using the Julia package manager:

```julia
using Pkg
Pkg.add("FullShell")
```

Or in the Julia REPL package mode (press `]`):

```
pkg> add FullShell
```

To install the development version from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/CarlosP24/FullShell.jl.git")
```

## Quick Start

```julia
using FullShell

# Define model parameters (units: meV and nm)
model = (
    R = 70,      # Cylinder radius (nm)
    w = 10,      # Semiconductor width (nm)
    d = 10,      # Superconductor width (nm)
    Δ0 = 0.23,   # Superconducting gap (meV)
    μ = 0.75,    # Chemical potential (meV)
    α = 0.2,     # Rashba SOC (meV·nm)
    Φ = 1.0,     # Magnetic flux (in units of Φ₀)
)

# Build Hamiltonians
hSM, hSC, parameters = build_cyl(; model...)

# hSM: Semiconductor in normal state
# hSC: Semiconductor with proximity-induced superconductivity
```

## Physical Model

The Hamiltonian describes a cylindrical semiconductor nanowire with inner radius R and width w,
covered by a superconducting shell of width d. The model includes:

1. **Kinetic energy**: Tight-binding with effective mass m₀
2. **Electrostatic potential**: Optional dome-shaped profile V(ρ)
3. **Spin-orbit coupling**: Rashba (and optionally Dresselhaus) terms
4. **Zeeman coupling**: g-factor coupling to magnetic field
5. **Orbital effects**: Vector potential in cylindrical gauge
6. **Proximity effect**: Superconducting pairing induced via proximity effect

The superconducting proximity effect can be modeled using:
- **Usadel formalism**: For diffusive superconductors (default)
- **Ballistic formalism**: For clean superconductors
- **Simplified model**: With phenomenological gap suppression

## Documentation Contents

```@contents
Pages = ["examples.md", "api.md"]
Depth = 2
```

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

This package is licensed under the terms specified in the LICENSE file.

## Related Work

FullShell.jl is built on top of [Quantica.jl](https://github.com/pablosanjose/Quantica.jl),
a powerful framework for quantum tight-binding models.
