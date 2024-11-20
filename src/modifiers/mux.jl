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
    hf = h |> supercell(region = r -> 0 <= r[1] <= L)
    mod! = @onsite!((o, r) -> o + f(r[1]/L) * σ0τz)
    return hf |> mod!
end


@with_kw struct params_shift @deftype Float64
    L = 500
    ς = 0.1
    μshift = 0
end

"""
    mu_step(h::Quantica.AbstractHamiltonian1D; kw...)
    mu_step(h::Quantica.AbstractHamiltonian1D, p::params_shift)

Apply a smooth step-like chemical potential shift to a 1D Hamiltonian.

# Arguments
- `h::Quantica.AbstractHamiltonian1D`: The 1D Hamiltonian to be modified.
- `kw...`: Keyword arguments to customize the parameters of the step function.
- `p::params_shift`: A struct containing the parameters for the step function.

# Keyword Arguments
- `L::Real`: The length of the 0D supercell (default: 500).
- `ς::Real`: The smoothness parameter of the step function (default: 0.1).
- `μshift::Real`: The magnitude of the chemical potential shift (default: 0).

# Returns
- The modified Hamiltonian.
"""

mu_step(h; kw...) = mu_step(h, params_shift(; kw...))
function mu_step(h::Quantica.AbstractHamiltonian1D, p::params_shift)
    @unpack L, ς, μshift = p
    step(x) = ifelse(ς == 0, sign(x),  0.5 * (1 + tanh(x/ς)))
    μx(x) = μshift * (1 - step(x - 1/2))
    return mux(h, μx, L)
end