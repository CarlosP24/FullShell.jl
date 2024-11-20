"""
    mux(h::Quantica.AbstractHamiltonian1D, f::Function, L::Real)

Create a 0D supercell of a 1D Hamiltonian applying a function `f(x)` to the onsite terms.
The wire is assumed to be 1D in the r[1] direction.

# Arguments
- `h::Quantica.AbstractHamiltonian1D`: The 1D Hamiltonian to be modified.
- `f::Function`: A function of one positional variable normalized to L to be applied to the onsite terms.
- `L::Real`: The length of the 0D supercell.

# Throws
- `ArgumentError`: If `f` is not a function of one positional variable.

# Returns
- The modified Hamiltonian.
"""
function mux(h::Quantica.AbstractHamiltonian1D, f::Function, L::Real;)
    # Check if f is a function of one positional variable
    if length(methods(f).ms[1].sig.parameters) != 2
        throw(ArgumentError("f must be a function of one positional variable"))
    end
    br = h.h.lattice.bravais.matrix
    idir = findfirst(x -> x != 0, br)
    hf = h |> supercell(region = r -> 0 <= r[1] <= L)
    mod! = @onsite!((o, r) -> o + f(r[1]/L) * σ0τz)
    return hf |> mod!
end