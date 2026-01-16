"""
Pauli matrices and Nambu-space operators for full-shell Hamiltonians.

This module defines tensor products of Pauli matrices in spin (σ) and Nambu (τ) spaces,
which form the basis for constructing BdG (Bogoliubov-de Gennes) Hamiltonians.

# Nambu-Spin Basis
The basis is ordered as: |↑e⟩, |↓e⟩, |↑h⟩, |↓h⟩, where:
- ↑,↓ are spin states
- e,h are electron/hole (particle/antiparticle) states

# Exported Operators
## Full 4×4 matrices (σ ⊗ τ)
- `σ0τx`, `σ0τy`, `σ0τz`, `σ0τ0`: Identity in spin ⊗ Pauli in Nambu
- `σzτ0`, `σzτz`: Pauli-z in spin ⊗ Identity/Pauli-z in Nambu
- `σyτy`, `σyτz`, `σyτ0`: Pauli-y in spin ⊗ Pauli-y/z/identity in Nambu
- `σxτz`, `σxτ0`: Pauli-x in spin ⊗ Pauli-z/identity in Nambu

## 2×2 Pauli matrices
- `σ0`, `σx`, `σy`, `σz`: Standard Pauli matrices

## Projection operators
- `c_up`: Projects onto spin-up sector (electrons and holes)
- `c_down`: Projects onto spin-down sector (electrons and holes)
"""

const σ0τx = @SMatrix[0 0 1 0; 0 0 0 1; 1 0 0 0; 0 1 0 0]
const σ0τy = @SMatrix[0 0 -im 0; 0 0 0 -im; im 0 0 0; 0 im 0 0]
const σ0τz = @SMatrix[1 0 0 0; 0 1 0 0; 0 0 -1 0; 0 0 0 -1]
const σ0τ0 = @SMatrix[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]
const σzτ0 = @SMatrix[1 0 0 0; 0 -1 0 0; 0 0 1 0; 0 0 0 -1]
const σzτz = @SMatrix[1 0 0 0; 0 -1 0 0; 0 0 -1 0; 0 0 0 1]
const σyτy = @SMatrix[0 0 0 -1; 0 0 1 0; 0 1 0 0; -1 0 0 0]
const σyτz = @SMatrix[0 -im 0 0; im 0 0 0; 0 0 0 im; 0 0 -im 0]
const σyτ0 = @SMatrix[0 -im 0 0; im 0 0 0; 0 0 0 -im; 0 0 im 0]
const σxτz = @SMatrix[0 1 0 0; 1 0 0 0; 0 0 0 -1; 0 0 -1 0]
const σxτ0 = @SMatrix[0 1 0 0; 1 0 0 0; 0 0 0 1; 0 0 1 0]
const σ0 = SA[1 0; 0 1]
const σx = SA[0 1; 1 0]
const σy = SA[0 -im; im 0]
const σz = SA[1 0; 0 -1]

const c_up = @SMatrix[1 0 0 0; 0 1 0 0; 0 0 0 0; 0 0 0 0] 
const c_down = @SMatrix[0 0 0 0; 0 0 0 0; 0 0 1 0; 0 0 0 1] 