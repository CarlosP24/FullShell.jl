# FullShell
[![Build Status](https://github.com/CarlosP24/FullShell.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/CarlosP24/FullShell.jl/actions/workflows/CI.yml?query=branch%3Amain)

Package to build full-shell hybrid semiconductor-superconductor hamiltonians.

Uses [Quantica.jl](https://github.com/pablosanjose/Quantica.jl). 

## Example
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


