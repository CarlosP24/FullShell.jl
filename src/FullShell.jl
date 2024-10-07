module FullShell

using Quantica 
using Parameters
using Roots

export build_cyl, build_cyl_mm, ΣS3DUsadel, pairbreaking, LP_lobe, ΔΛ, itip

export Params, Params_mm

export σ0τx, σ0τy, σ0τz, σ0τ0, σzτ0, σzτz, σyτy, σyτz, σyτ0, σxτz, σxτ0, σ0, σx, σy, σz, c_up, c_down, Ω

include("utilities/LP_usadel.jl")
include("utilities/pauli_products.jl") 
include("HamiltonianBuilder.jl")
include("HamiltonianBuilder_MM.jl")

end
