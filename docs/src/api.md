# API Reference

Complete API documentation for FullShell.jl.

## Contents

```@contents
Pages = ["api.md"]
Depth = 3
```

## Main Module

```@docs
FullShell
```

## Hamiltonian Builders

```@docs
build_cyl
build_cyl_mm
```

## Parameter Structures

```@docs
Params
Params_mm
```

## Superconducting Proximity Functions

```@docs
pairbreaking
ΔD
Ω
LP_lobe
is_in_lobe
```

## Helper Functions

### Conversion Functions

```@docs
get_Φ
get_B
get_Ω
get_itip
```

### Analysis Functions

```@docs
bandwidth
```

## Hamiltonian Modifiers

```@docs
mux
mu_step
build_harmonic_deformations
```

## Pauli Matrices and Operators

The package exports tensor products of Pauli matrices in spin (σ) and Nambu (τ) spaces.

### 4×4 Operators (Spin ⊗ Nambu)

- **`σ0τx`, `σ0τy`, `σ0τz`, `σ0τ0`**: Identity in spin space ⊗ Pauli matrices in Nambu space
- **`σzτ0`, `σzτz`**: Pauli-z in spin space ⊗ Identity/Pauli-z in Nambu space  
- **`σyτy`, `σyτz`, `σyτ0`**: Pauli-y in spin space ⊗ Pauli matrices in Nambu space
- **`σxτz`, `σxτ0`**: Pauli-x in spin space ⊗ Pauli-z/Identity in Nambu space

### 2×2 Pauli Matrices

- **`σ0`**: 2×2 identity matrix
- **`σx`**: Pauli-x matrix
- **`σy`**: Pauli-y matrix
- **`σz`**: Pauli-z matrix

### Projection Operators

- **`c_up`**: Projects onto spin-up sector (both electron and hole)
- **`c_down`**: Projects onto spin-down sector (both electron and hole)

### Physical Interpretation

| Operator | Physical meaning |
|----------|------------------|
| `σ0τz` | Particle-hole asymmetric term (normal state Hamiltonian) |
| `σ0τx` | Real part of singlet s-wave pairing |
| `σ0τy` | Imaginary part of singlet s-wave pairing |
| `σzτ0` | Zeeman coupling to magnetic field |
| `σyτz` | Rashba spin-orbit coupling (y-component) |
| `σxτz` | Dresselhaus or Rashba spin-orbit coupling (x-component) |

### Basis Convention

The basis is ordered as: |↑e⟩, |↓e⟩, |↑h⟩, |↓h⟩, where:
- ↑, ↓ are spin-up and spin-down states
- e, h are electron and hole (particle/antiparticle) states in Nambu formalism

## Index

```@index
```
